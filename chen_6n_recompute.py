#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
6N prime-distribution recompute script  (chen_6n_recompute.py)
================================================================================
Purpose: fix the factor-truncation bug caused by sympy.factorint(N, limit=1000)
in the old exploratory code. That bug treated a prime factor of N above 1000
(or a product of them) as a single 'prime factor', systematically corrupting the alpha/omega statistics.

This script does a COMPLETE factorisation of every N via a segmented sieve, and a DETERMINISTIC
primality test of 6N-1, 6N+1 (interval sieve, no probabilistic step), computing per node:

  - whether 6N-1, 6N+1 are prime (-> left / right / twin / singleton)
  - the complete prime factorisation of N
  - omega_big(N) = number of DISTINCT prime factors > 3 of N   (used by the singular-series local factor)
  - alpha(N)     = sum of exponents of prime factors > 3 (old definition, kept for comparison)
  - E_old(N) = prod_{distinct p|N, p>3} p/(p-2)     (survival-only term; first-order approximation)
  - E_new(N) = prod_{distinct p|N, p>3} (p-1)/(p-3) (correct local factor with the deficit correction)

Output (CSV, mapping directly to the paper's tables):
  table1_shell_counts.csv      absolute shell counts (all nodes vs prime-bearing nodes)
  table2_axial_symmetry.csv    left/right wing symmetry ratio
  table3_conditional_prob.csv  twin conditional probability by omega_big + Wilson CI
  table5_enrichment_test.csv   core new analysis: observed twin rate vs singular-series factor prediction
  table4_mean_omega.csv        mean omega/alpha of twin vs singleton nodes

Requires: numpy. sympy is only used for the optional self-check.
Usage:
    python chen_6n_recompute.py            # default to S8 (minutes; validate the pipeline first)
    MAX_K=10 python chen_6n_recompute.py   # full S10 (long; run in background)
    VERIFY=1 python chen_6n_recompute.py   # run the factorisation self-check first

Author's note: this script does not 'prove' any conjecture. It only produces reproducible, correctly factored numerical data.
================================================================================
"""

import os
import sys
import time
import math
import csv

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required:  pip install numpy")

# ------------------------------- config -------------------------------------
MAX_K     = int(os.environ.get("MAX_K", 8))          # deepest shell S_{MAX_K}
SEGMENT   = int(os.environ.get("SEGMENT", 2_000_000))# segment size in N (memory/speed trade-off)
OUTDIR    = os.environ.get("OUTDIR", "./recompute_out")
DO_VERIFY = bool(int(os.environ.get("VERIFY", 0)))   # whether to self-check factorisation on the first 1e5 N
os.makedirs(OUTDIR, exist_ok=True)

N_MAX = (10**MAX_K) // 6                              # max N (so that 6N < 10^MAX_K)
PRIME_BOUND = int(math.isqrt(10**MAX_K)) + 1          # base-prime bound for the sieve = sqrt(10^MAX_K)

# --------------------------- base primes (once) -----------------------------
def primes_upto(n):
    """Simple sieve of Eratosthenes; returns all primes <= n (numpy int64 array)."""
    if n < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.nonzero(sieve)[0].astype(np.int64)

print(f"[setup] generating base primes (<= {PRIME_BOUND}) ...")
BASE_PRIMES = primes_upto(PRIME_BOUND)
print(f"[setup] base prime count: {len(BASE_PRIMES)};  N_MAX = {N_MAX:,};  MAX_K = {MAX_K}")

# ----------------------- shell assignment (by number of digits of 6N) -------
# 6N in [10^{K-1}, 10^K)  <=>  K = number of decimal digits
POW10 = np.array([10**i for i in range(0, MAX_K + 2)], dtype=np.int64)

def shell_of(values6n):
    """Return the shell K of each 6N value (K such that 10^{K-1} <= 6N < 10^K)."""
    return np.searchsorted(POW10, values6n, side='right')

# ============= segmented: complete factorisation + interval primality =======
def factor_segment(n_lo, n_hi):
    """
    Complete factorisation of N in [n_lo, n_hi). Returns:
      omega_big : number of distinct prime factors > 3   (int16)
      alpha     : sum of exponents of prime factors > 3  (int16)
      enrich    : prod p/(p-2)        (float64)
    """
    size = n_hi - n_lo
    rem       = np.arange(n_lo, n_hi, dtype=np.int64)   # remainder after dividing out found factors
    omega_big = np.zeros(size, dtype=np.int16)
    alpha     = np.zeros(size, dtype=np.int16)
    enrich    = np.ones(size, dtype=np.float64)   # old factor prod p/(p-2)
    enrich2   = np.ones(size, dtype=np.float64)   # correct factor prod (p-1)/(p-3)

    for p in BASE_PRIMES:
        if p * p > (n_hi - 1):       # beyond sqrt(max): the remainder can only be a single large prime
            break
        first = ((n_lo + p - 1) // p) * p     # first multiple of p in the segment
        if first >= n_hi:
            continue
        idx = np.arange(first - n_lo, size, p)
        if idx.size == 0:
            continue
        sub = rem[idx]
        cnt = np.zeros(idx.size, dtype=np.int16)
        m = (sub % p) == 0
        while m.any():               # divide out all powers of p, count the exponent
            sub[m] //= p
            cnt[m] += 1
            m = (sub % p) == 0
        rem[idx] = sub
        if p > 3:                    # 2 and 3 are 'substrate factors', excluded from omega/alpha/enrich
            omega_big[idx] += 1
            alpha[idx]     += cnt
            enrich[idx]    *= p / (p - 2.0)
            enrich2[idx]   *= (p - 1.0) / (p - 3.0)

    # a remainder rem>1 must be a single prime > sqrt(max) (hence > 3)
    leftover = rem > 1
    if leftover.any():
        omega_big[leftover] += 1
        alpha[leftover]     += 1
        lp = rem[leftover].astype(np.float64)
        enrich[leftover]    *= lp / (lp - 2.0)
        enrich2[leftover]   *= (lp - 1.0) / (lp - 3.0)
    return omega_big, alpha, enrich, enrich2


def primality_6n(n_lo, n_hi):
    """
    Deterministically test whether 6N-1, 6N+1 (N in [n_lo,n_hi)) are prime via an interval sieve.
    Returns two boolean arrays (is_prime_minus, is_prime_plus).
    """
    v_lo = 6 * n_lo - 1
    v_hi = 6 * (n_hi - 1) + 1
    span = v_hi - v_lo + 1
    comp = np.zeros(span, dtype=bool)       # comp[i] = (v_lo+i) is composite
    sq = int(math.isqrt(v_hi)) + 1
    for p in BASE_PRIMES:
        if p > sq:
            break
        start = max(p * p, ((v_lo + p - 1) // p) * p)
        if start > v_hi:
            continue
        comp[start - v_lo : span : p] = True
    N = np.arange(n_lo, n_hi, dtype=np.int64)
    m_minus = 6 * N - 1
    m_plus  = 6 * N + 1
    pm = (~comp[m_minus - v_lo]) & (m_minus > 1)
    pp = (~comp[m_plus  - v_lo]) & (m_plus  > 1)
    return pm, pp

# ------------------------------ self-check (optional) -----------------------
def verify():
    """Compare the segmented factorisation with sympy on the first 100000 N to ensure no truncation error."""
    try:
        from sympy import factorint, isprime
    except ImportError:
        print("[verify] sympy not installed, skipping self-check.")
        return
    print("[verify] comparing factorisation and primality of the first 100000 N against sympy ...")
    n_lo, n_hi = 1, 100001
    ob, al, en = factor_segment(n_lo, n_hi)
    pm, pp = primality_6n(n_lo, n_hi)
    bad = 0
    for i, N in enumerate(range(n_lo, n_hi)):
        f = factorint(N)
        ob_t = sum(1 for q in f if q > 3)
        al_t = sum(e for q, e in f.items() if q > 3)
        en_t = 1.0
        for q in f:
            if q > 3:
                en_t *= q / (q - 2.0)
        if ob[i] != ob_t or al[i] != al_t or abs(en[i] - en_t) > 1e-9 * max(1, en_t):
            bad += 1
            if bad <= 5:
                print(f"   mismatch N={N}: script(ob={ob[i]},al={al[i]},E={en[i]:.4f}) "
                      f"vs sympy(ob={ob_t},al={al_t},E={en_t:.4f})")
        if pm[i] != isprime(6*N-1) or pp[i] != isprime(6*N+1):
            bad += 1
            if bad <= 5:
                print(f"   primality mismatch N={N}")
    print(f"[verify] done. mismatches: {bad} (should be 0).")

# ------------------------- Wilson binomial confidence interval --------------
def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z*z/n
    c = p + z*z/(2*n)
    h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return ((c - h)/d, (c + h)/d)

# =============================== main =======================================
def main():
    if DO_VERIFY:
        verify()

    NS = MAX_K + 2
    OW_MAX = 20  # omega_big histogram upper bound

    # shell-level accumulators (indexed by K)
    s_total = np.zeros(NS, dtype=np.int64)   # all nodes
    s_left  = np.zeros(NS, dtype=np.int64)   # only 6N-1 prime
    s_right = np.zeros(NS, dtype=np.int64)   # only 6N+1 prime
    s_twin  = np.zeros(NS, dtype=np.int64)   # twin
    s_alpha_single = np.zeros(NS, dtype=np.float64)  # sum of alpha over singleton nodes
    s_alpha_twin   = np.zeros(NS, dtype=np.float64)
    s_omega_single = np.zeros(NS, dtype=np.float64)
    s_omega_twin   = np.zeros(NS, dtype=np.float64)

    # stratified by omega_big (global)
    w_nodes = np.zeros(OW_MAX + 1, dtype=np.int64)
    w_twin  = np.zeros(OW_MAX + 1, dtype=np.int64)

    # enrichment test: accumulate nodes / twins / sum(E) by (shell K, omega_big)
    et_nodes = np.zeros((NS, OW_MAX + 1), dtype=np.int64)
    et_twin  = np.zeros((NS, OW_MAX + 1), dtype=np.int64)
    et_sumE  = np.zeros((NS, OW_MAX + 1), dtype=np.float64)   # sum of the old factor
    et_sumE2 = np.zeros((NS, OW_MAX + 1), dtype=np.float64)   # sum of the correct factor

    t0 = time.time()
    n = 1
    seg_id = 0
    while n <= N_MAX:
        n_hi = min(n + SEGMENT, N_MAX + 1)
        ob, al, en, en2 = factor_segment(n, n_hi)
        pm, pp = primality_6n(n, n_hi)

        N = np.arange(n, n_hi, dtype=np.int64)
        K = shell_of(6 * N)
        twin   = pm & pp
        single = (pm ^ pp)            # exactly one wing prime
        any_pr = pm | pp

        # shell accumulation
        np.add.at(s_total, K, 1)
        np.add.at(s_left,  K[pm & ~pp], 1)
        np.add.at(s_right, K[pp & ~pm], 1)
        np.add.at(s_twin,  K[twin], 1)
        np.add.at(s_alpha_twin,   K[twin],   al[twin])
        np.add.at(s_omega_twin,   K[twin],   ob[twin])
        np.add.at(s_alpha_single, K[single], al[single])
        np.add.at(s_omega_single, K[single], ob[single])

        # omega_big stratification (base = prime-bearing nodes only)
        obc = np.clip(ob, 0, OW_MAX)
        np.add.at(w_nodes, obc[any_pr], 1)
        np.add.at(w_twin,  obc[twin], 1)

        # enrichment-test accumulation (base = all nodes; studying P(twin | N))
        flat = K * (OW_MAX + 1) + obc
        np.add.at(et_nodes.ravel(), flat, 1)
        np.add.at(et_twin.ravel(),  flat[twin], 1)
        np.add.at(et_sumE.ravel(),  flat, en)
        np.add.at(et_sumE2.ravel(), flat, en2)

        seg_id += 1
        if seg_id % 20 == 0 or n_hi > N_MAX:
            done = n_hi - 1
            rate = done / max(1e-9, time.time() - t0)
            print(f"[scan] N={done:,}/{N_MAX:,}  ({100*done/N_MAX:5.1f}%)  "
                  f"{rate:,.0f} N/s  cumulative twins={s_twin.sum():,}")
        n = n_hi

    elapsed = time.time() - t0
    print(f"[done] scan complete in {elapsed:.1f}s. Writing results ...")

    # ----------------------------- write tables ----------------------------
    # table 1: shell counts
    with open(f"{OUTDIR}/table1_shell_counts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shell_K", "total_nodes", "prime_bearing_nodes",
                    "single_primes", "twin_pairs", "growth_factor_twin"])
        prev_twin = None
        for K in range(1, MAX_K + 1):
            single = s_left[K] + s_right[K]
            twin = s_twin[K]
            pb = single + twin
            gf = (twin / prev_twin) if (prev_twin and prev_twin > 0) else ""
            gf = f"{gf:.2f}" if gf != "" else ""
            w.writerow([f"S{K}", s_total[K], pb, single, twin, gf])
            prev_twin = twin

    # table 2: left/right wing symmetry
    with open(f"{OUTDIR}/table2_axial_symmetry.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shell_K", "left_6N-1", "right_6N+1", "ratio_L/R", "abs_dev"])
        for K in range(1, MAX_K + 1):
            L = s_left[K] + s_twin[K]   # total left-wing primes = left-only + twin
            R = s_right[K] + s_twin[K]
            ratio = (L / R) if R else 0.0
            w.writerow([f"S{K}", L, R, f"{ratio:.5f}", f"{abs(ratio-1):.2e}"])

    # table 3: twin conditional probability by omega_big + Wilson CI
    with open(f"{OUTDIR}/table3_conditional_prob.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["omega_big", "prime_bearing_nodes", "twin_pairs",
                    "cond_prob_%", "wilson_lo_%", "wilson_hi_%"])
        for ww in range(0, OW_MAX + 1):
            nn, tt = int(w_nodes[ww]), int(w_twin[ww])
            if nn == 0:
                continue
            lo, hi = wilson(tt, nn)
            w.writerow([ww, nn, tt, f"{100*tt/nn:.3f}",
                        f"{100*lo:.3f}", f"{100*hi:.3f}"])

    # table 5: enrichment test (observed vs singular-series factor prediction)
    # within a shell, 6N is ~constant in magnitude -> P(twin|N) ~ baseline_K * E(N).
    # baseline_K is set from all nodes of the shell: baseline_K = (total twins) / (sum E).
    # so predicted twin rate (per stratum) = baseline_K * meanE(stratum), compared to observed.
    with open(f"{OUTDIR}/table5_enrichment_test.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shell_K", "omega_big", "nodes", "twins", "obs_rate_%",
                    "meanE_old", "obs/pred_old", "meanE_new", "obs/pred_new"])
        for K in range(1, MAX_K + 1):
            tot_twin = et_twin[K].sum()
            sE  = et_sumE[K].sum()
            sE2 = et_sumE2[K].sum()
            if tot_twin == 0 or sE <= 0:
                continue
            base_old = tot_twin / sE
            base_new = tot_twin / sE2
            for ww in range(0, OW_MAX + 1):
                nn = int(et_nodes[K, ww])
                if nn < 50:
                    continue
                tt = int(et_twin[K, ww])
                obs = tt / nn
                mEo = et_sumE[K, ww] / nn
                mEn = et_sumE2[K, ww] / nn
                po = base_old * mEo
                pn = base_new * mEn
                w.writerow([f"S{K}", ww, nn, tt, f"{100*obs:.3f}",
                            f"{mEo:.4f}", f"{obs/po:.3f}" if po>0 else "",
                            f"{mEn:.4f}", f"{obs/pn:.3f}" if pn>0 else ""])

    # table 4: mean omega / alpha of twin vs singleton
    with open(f"{OUTDIR}/table4_mean_omega.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shell_K",
                    "mean_omega_single", "mean_omega_twin",
                    "mean_alpha_single", "mean_alpha_twin"])
        for K in range(1, MAX_K + 1):
            ns = s_left[K] + s_right[K]
            nt = s_twin[K]
            mos = (s_omega_single[K] / ns) if ns else 0.0
            mot = (s_omega_twin[K] / nt) if nt else 0.0
            mas = (s_alpha_single[K] / ns) if ns else 0.0
            mat = (s_alpha_twin[K] / nt) if nt else 0.0
            w.writerow([f"S{K}", f"{mos:.3f}", f"{mot:.3f}",
                        f"{mas:.3f}", f"{mat:.3f}"])

    print(f"[ok] all results written to {OUTDIR}/  (table1..table5)")
    print("     focus on table5_enrichment_test.csv -- if obs/pred is stably near 1.0,")
    print("     it is the empirical verification of the singular-series factor structure.")


if __name__ == "__main__":
    main()
