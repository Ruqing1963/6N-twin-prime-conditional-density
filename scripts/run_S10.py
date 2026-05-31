#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
6N 素数分布重算脚本  (chen_6n_recompute.py)
================================================================================
目的：修复旧勘探代码中 sympy.factorint(N, limit=1000) 造成的【因子截断 bug】。
该 bug 会把 N 中大于 1000 的素因子（或它们的乘积）错误地当成单个"素因子"，
导致 alpha / omega 统计系统性错误，使第三章的逐因子分析失去依据。

本脚本用【分段筛法】对每个 N 做*完整*素因子分解，并对 6N-1、6N+1 做*确定性*
素性判定（区间筛，无概率成分），逐节点计算：

  - 6N-1, 6N+1 是否素数（→ 左翼 / 右翼 / 孪生 / 单体）
  - N 的完整素因子分解
  - omega_big(N) = N 的【相异】大于3素因子个数      ← 奇异级数局部因子用这个
  - alpha(N)     = 大于3素因子的【幂次总和】（旧定义，保留以便对照）
  - E(N) = ∏_{相异 p|N, p>3} p/(p-2)                ← Hardy–Littlewood 奇异级数的逐因子放大率

输出（CSV，直接对应论文表格）：
  table1_shell_counts.csv      壳层绝对计数（区分"全部节点"与"含素数节点"）
  table2_axial_symmetry.csv    左右翼对称比
  table3_conditional_prob.csv  按 omega_big 分层的孪生条件概率 + Wilson 置信区间
  table5_enrichment_test.csv   ★核心新分析：实测孪生率 vs 奇异级数逐因子预测
  table4_mean_omega.csv         孪生节点 vs 单体节点的平均 omega/alpha

依赖：numpy（必需）。sympy 仅用于可选自检。
用法：
    python chen_6n_recompute.py            # 默认跑到 S8（几分钟，先验证管线）
    MAX_K=10 python chen_6n_recompute.py   # 跑完整 S10（耗时较长，建议放后台）
    VERIFY=1 python chen_6n_recompute.py   # 先做分解自检再跑

作者注：本脚本不"证明"任何猜想。它只产出可复现的、分解正确的数值数据。
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
    sys.exit("需要 numpy：  pip install numpy")

# ------------------------------- 配置 ---------------------------------------
MAX_K     = int(os.environ.get("MAX_K", 10))  # 本脚本默认 S10（完整）          # 最深壳层 S_{MAX_K}
SEGMENT   = int(os.environ.get("SEGMENT", 2_000_000))# N 的分段大小（内存/速度权衡）
OUTDIR    = os.environ.get("OUTDIR", "./recompute_out_S10")
DO_VERIFY = bool(int(os.environ.get("VERIFY", 0)))   # 是否对前 1e5 个 N 做分解自检
os.makedirs(OUTDIR, exist_ok=True)

N_MAX = (10**MAX_K) // 6                              # 最大 N（使 6N < 10^MAX_K）
PRIME_BOUND = int(math.isqrt(10**MAX_K)) + 1          # 筛法基素数上界 = sqrt(10^MAX_K)

# --------------------------- 基素数（一次性） --------------------------------
def primes_upto(n):
    """简单埃氏筛，返回 <= n 的全部素数（numpy int64 数组）。"""
    if n < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.nonzero(sieve)[0].astype(np.int64)

print(f"[setup] 生成基素数 (<= {PRIME_BOUND}) ...")
BASE_PRIMES = primes_upto(PRIME_BOUND)
print(f"[setup] 基素数个数: {len(BASE_PRIMES)};  N_MAX = {N_MAX:,};  MAX_K = {MAX_K}")

# ----------------------- 壳层归属（按 6N 的位数） ----------------------------
# 6N ∈ [10^{K-1}, 10^K)  ⇔  K = 该数的十进制位数
POW10 = np.array([10**i for i in range(0, MAX_K + 2)], dtype=np.int64)

def shell_of(values6n):
    """返回每个 6N 值所属壳层 K（K 使 10^{K-1} <= 6N < 10^K）。"""
    return np.searchsorted(POW10, values6n, side='right')

# =================== 分段：完整因子分解 + 区间素性筛 =========================
def factor_segment(n_lo, n_hi):
    """
    对 N ∈ [n_lo, n_hi) 做完整分解，返回:
      omega_big : 相异(>3)素因子个数   (int16)
      alpha     : (>3)素因子幂次总和   (int16)
      enrich    : ∏ p/(p-2)            (float64)
    """
    size = n_hi - n_lo
    rem       = np.arange(n_lo, n_hi, dtype=np.int64)   # 逐步除尽后的剩余
    omega_big = np.zeros(size, dtype=np.int16)
    alpha     = np.zeros(size, dtype=np.int16)
    enrich    = np.ones(size, dtype=np.float64)

    for p in BASE_PRIMES:
        if p * p > (n_hi - 1):       # 超过 sqrt(max)：剩下的只能是单个大素数
            break
        first = ((n_lo + p - 1) // p) * p     # 段内第一个 p 的倍数
        if first >= n_hi:
            continue
        idx = np.arange(first - n_lo, size, p)
        if idx.size == 0:
            continue
        sub = rem[idx]
        cnt = np.zeros(idx.size, dtype=np.int16)
        m = (sub % p) == 0
        while m.any():               # 除尽 p 的所有幂次，统计指数
            sub[m] //= p
            cnt[m] += 1
            m = (sub % p) == 0
        rem[idx] = sub
        if p > 3:                    # 2、3 是"基质因子"，不计入 omega/alpha/enrich
            omega_big[idx] += 1
            alpha[idx]     += cnt
            enrich[idx]    *= p / (p - 2.0)

    # 剩余 rem>1 必为单个大于 sqrt(max) 的素数（>3）
    leftover = rem > 1
    if leftover.any():
        omega_big[leftover] += 1
        alpha[leftover]     += 1
        lp = rem[leftover].astype(np.float64)
        enrich[leftover]    *= lp / (lp - 2.0)
    return omega_big, alpha, enrich


def primality_6n(n_lo, n_hi):
    """
    用区间筛确定性判定 6N-1, 6N+1 (N∈[n_lo,n_hi)) 是否素数。
    返回两个布尔数组 (is_prime_minus, is_prime_plus)。
    """
    v_lo = 6 * n_lo - 1
    v_hi = 6 * (n_hi - 1) + 1
    span = v_hi - v_lo + 1
    comp = np.zeros(span, dtype=bool)       # comp[i] = (v_lo+i) 是合数
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

# ------------------------------ 自检（可选） --------------------------------
def verify():
    """对前 100000 个 N，比对分段分解与 sympy 完整分解，确保无截断错误。"""
    try:
        from sympy import factorint, isprime
    except ImportError:
        print("[verify] 未安装 sympy，跳过自检。")
        return
    print("[verify] 用 sympy 比对前 100000 个 N 的分解与素性 ...")
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
                print(f"   不一致 N={N}: 脚本(ob={ob[i]},al={al[i]},E={en[i]:.4f}) "
                      f"vs sympy(ob={ob_t},al={al_t},E={en_t:.4f})")
        if pm[i] != isprime(6*N-1) or pp[i] != isprime(6*N+1):
            bad += 1
            if bad <= 5:
                print(f"   素性不一致 N={N}")
    print(f"[verify] 完成。不一致数: {bad}（应为 0）。")

# ------------------------- Wilson 二项置信区间 ------------------------------
def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z*z/n
    c = p + z*z/(2*n)
    h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return ((c - h)/d, (c + h)/d)

# =============================== 主流程 =====================================
def main():
    if DO_VERIFY:
        verify()

    NS = MAX_K + 2
    OW_MAX = 20  # omega_big 直方图上界

    # 壳层级累加器（按 K 索引）
    s_total = np.zeros(NS, dtype=np.int64)   # 全部节点
    s_left  = np.zeros(NS, dtype=np.int64)   # 仅 6N-1 素
    s_right = np.zeros(NS, dtype=np.int64)   # 仅 6N+1 素
    s_twin  = np.zeros(NS, dtype=np.int64)   # 孪生
    s_alpha_single = np.zeros(NS, dtype=np.float64)  # 单体节点 alpha 之和
    s_alpha_twin   = np.zeros(NS, dtype=np.float64)
    s_omega_single = np.zeros(NS, dtype=np.float64)
    s_omega_twin   = np.zeros(NS, dtype=np.float64)

    # 按 omega_big 分层（全局）
    w_nodes = np.zeros(OW_MAX + 1, dtype=np.int64)
    w_twin  = np.zeros(OW_MAX + 1, dtype=np.int64)

    # ★ 富集检验：按 (壳层K, omega_big) 累加 节点数 / 孪生数 / E之和
    et_nodes = np.zeros((NS, OW_MAX + 1), dtype=np.int64)
    et_twin  = np.zeros((NS, OW_MAX + 1), dtype=np.int64)
    et_sumE  = np.zeros((NS, OW_MAX + 1), dtype=np.float64)

    t0 = time.time()
    n = 1
    seg_id = 0
    while n <= N_MAX:
        n_hi = min(n + SEGMENT, N_MAX + 1)
        ob, al, en = factor_segment(n, n_hi)
        pm, pp = primality_6n(n, n_hi)

        N = np.arange(n, n_hi, dtype=np.int64)
        K = shell_of(6 * N)
        twin   = pm & pp
        single = (pm ^ pp)            # 恰好一翼为素
        any_pr = pm | pp

        # 壳层累加
        np.add.at(s_total, K, 1)
        np.add.at(s_left,  K[pm & ~pp], 1)
        np.add.at(s_right, K[pp & ~pm], 1)
        np.add.at(s_twin,  K[twin], 1)
        np.add.at(s_alpha_twin,   K[twin],   al[twin])
        np.add.at(s_omega_twin,   K[twin],   ob[twin])
        np.add.at(s_alpha_single, K[single], al[single])
        np.add.at(s_omega_single, K[single], ob[single])

        # omega_big 分层（仅统计"含素数节点"作为孪生/非孪生的基底）
        obc = np.clip(ob, 0, OW_MAX)
        np.add.at(w_nodes, obc[any_pr], 1)
        np.add.at(w_twin,  obc[twin], 1)

        # 富集检验累加（基底 = 全部节点，研究 P(twin | N)）
        flat = K * (OW_MAX + 1) + obc
        np.add.at(et_nodes.ravel(), flat, 1)
        np.add.at(et_twin.ravel(),  flat[twin], 1)
        np.add.at(et_sumE.ravel(),  flat, en)

        seg_id += 1
        if seg_id % 20 == 0 or n_hi > N_MAX:
            done = n_hi - 1
            rate = done / max(1e-9, time.time() - t0)
            print(f"[scan] N={done:,}/{N_MAX:,}  ({100*done/N_MAX:5.1f}%)  "
                  f"{rate:,.0f} N/s  累计孪生={s_twin.sum():,}")
        n = n_hi

    elapsed = time.time() - t0
    print(f"[done] 扫描完成，用时 {elapsed:.1f}s。开始写出结果 ...")

    # ----------------------------- 写出表格 --------------------------------
    # 表1：壳层计数
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

    # 表2：左右翼对称
    with open(f"{OUTDIR}/table2_axial_symmetry.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shell_K", "left_6N-1", "right_6N+1", "ratio_L/R", "abs_dev"])
        for K in range(1, MAX_K + 1):
            L = s_left[K] + s_twin[K]   # 左翼素数总数 = 仅左 + 孪生
            R = s_right[K] + s_twin[K]
            ratio = (L / R) if R else 0.0
            w.writerow([f"S{K}", L, R, f"{ratio:.5f}", f"{abs(ratio-1):.2e}"])

    # 表3：按 omega_big 的孪生条件概率 + Wilson CI
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

    # ★表5：富集检验（实测 vs 奇异级数逐因子预测）
    # 在每个壳层内，6N 量级近似恒定 → P(twin|N) ≈ baseline_K * E(N)。
    # baseline_K 由该壳层全部节点定标： baseline_K = (总孪生数) / (ΣE)。
    # 于是 预测孪生率(分层) = baseline_K * meanE(分层)，与实测率对比。
    with open(f"{OUTDIR}/table5_enrichment_test.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shell_K", "omega_big", "nodes", "twins",
                    "obs_rate_%", "mean_E", "pred_rate_%", "obs/pred"])
        for K in range(1, MAX_K + 1):
            tot_twin = et_twin[K].sum()
            tot_sumE = et_sumE[K].sum()
            if tot_sumE <= 0 or tot_twin == 0:
                continue
            baseline = tot_twin / tot_sumE         # 每单位 E 的孪生率
            for ww in range(0, OW_MAX + 1):
                nn = int(et_nodes[K, ww])
                if nn < 50:                        # 样本太小不报告
                    continue
                tt = int(et_twin[K, ww])
                meanE = et_sumE[K, ww] / nn
                obs = tt / nn
                pred = baseline * meanE
                ratio = (obs / pred) if pred > 0 else 0.0
                w.writerow([f"S{K}", ww, nn, tt, f"{100*obs:.3f}",
                            f"{meanE:.4f}", f"{100*pred:.3f}", f"{ratio:.3f}"])

    # 表4：孪生 vs 单体的平均 omega / alpha
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

    print(f"[ok] 全部结果写入 {OUTDIR}/  (table1..table5)")
    print("     重点看 table5_enrichment_test.csv —— obs/pred 若稳定接近 1.0，")
    print("     即为奇异级数逐因子结构的实测验证（论文第三章的核心证据）。")


if __name__ == "__main__":
    main()
