# 6N Twin-Prime Conditional Density

Exhaustive scan of prime occurrence on the **6N ± 1** skeleton up to **10¹⁰**
(all 1.67 × 10⁹ centres, 2.74 × 10⁷ twin pairs), together with a
**factor-resolved local model** for the conditional density of twin primes:
the conditional twin density at a centre *N* is proportional to

> **E(N) = ∏_{q | N, q > 3} (q − 1)/(q − 3).**

This repository contains the data, code, figures, and the manuscript that
together support the result. The model is validated factor by factor against
~10⁸ nodes and reproduces the squarefree bias reported by Puszkarz (2018) to
within −0.2 % at the matching 10¹⁰ scale.

> **Scope.** This is experimental / computational number theory. The model is a
> validated **Hardy–Littlewood heuristic**, *not* an unconditional theorem.
> **Nothing here bears on the infinitude of twin primes.** Priority for the
> underlying phenomenon ("twins prefer factor-rich centres") belongs to
> W. Puszkarz, arXiv:1807.00406 (2018); the cumulative twin counts coincide with
> OEIS A007508. The contribution claimed is the *factor-resolved local model and
> its verification*, nothing more.

---

## Repository layout

```
.
├── README.md
├── LICENSE                  (MIT)
├── CITATION.cff
├── data/                    five result tables (S1–S10, complete recomputation)
│   ├── table1_shell_counts_S10.csv
│   ├── table2_axial_symmetry_S10.csv
│   ├── table3_conditional_prob_S10.csv
│   ├── table4_mean_omega_S10.csv
│   ├── table5_enrichment_test_S10.csv
│   └── table5_enrichment_test_S8_dualmodel.csv
├── code/
│   ├── chen_6n_recompute.py        main scan: complete factorisation + interval-sieve primality
│   ├── run_S8.py                   convenience launcher (S8, ~minutes)
│   ├── run_S10.py                  convenience launcher (S10, ~hours)
│   ├── puszkarz_S10_correct.py     Puszkarz-bias recovery at 10^10 (§3.4)
│   ├── fig_enrichment_test.py      regenerates Figure 1
│   └── fig_puszkarz_recovery.py    regenerates Figure 2
├── figures/                 fig_enrichment_test.{pdf,png}, fig_puszkarz_recovery.{pdf,png}
└── paper/                   Chen_6N_TwinPrimes.{tex,pdf} + figure PDFs
```

---

## Reproducing the results

Requirements: Python 3.8+, `numpy`. Optional: `sympy` (self-check),
`matplotlib` (figures). A LaTeX install rebuilds the paper.

```bash
pip install numpy sympy matplotlib

# 1. Main scan. Writes table1..table5 CSVs into recompute_out_S8/ (or _S10/).
python code/run_S8.py            # S1–S8, a few minutes  (good first check)
python code/run_S10.py           # S1–S10, a few hours    (full data)

# Optional integrity self-check (compares complete factorisation + primality
# against sympy on the first 100,000 centres; should report 0 discrepancies):
VERIFY=1 python code/run_S8.py

# 2. Puszkarz-bias recovery at the 10^10 scale (produces rows (a),(b),(c) of §3.4):
python code/puszkarz_S10_correct.py

# 3. Figures:
python code/fig_enrichment_test.py
python code/fig_puszkarz_recovery.py
```

### Note on an earlier bug (why "complete factorisation" matters)

An earlier exploratory script factored each *N* only up to a small bound
(`factorint(N, limit=1000)`), which silently misreports any *N* with a prime
factor > 1000 and corrupts the factor counts ω, α. **All code here uses complete
factorisation** via a segmented sieve, self-verified against `sympy` with zero
discrepancies. Any analysis stratified by the factorisation of *N* must use the
complete factorisation; results from the truncated version are not valid.

---

## Data dictionary

**table1_shell_counts** — one row per shell S_K = [10^{K−1}, 10^K).
`total_nodes` (all 6N centres in the shell), `prime_bearing_nodes` (≥ one of
6N±1 prime), `single_primes`, `twin_pairs`, `growth_factor_twin`
(T₂(S_K)/T₂(S_{K−1})).
*Validation:* cumulative `twin_pairs` = OEIS A007508 minus 1 (the excluded pair
(3,5)) at every K.

**table2_axial_symmetry** — `left_6N-1`, `right_6N+1` prime counts per shell,
their ratio `ratio_L/R`, and `abs_dev` = |ratio − 1|. Deviation decays as
O(N^{−1/2}) (central-limit scale).

**table3_conditional_prob** — twin conditional probability stratified by
`omega_big` = ω_{>3}(N) (distinct prime factors of N exceeding 3), with
`prime_bearing_nodes`, `twin_pairs`, `cond_prob_%`, and Wilson 95 % interval
`wilson_lo_%`, `wilson_hi_%`.

**table4_mean_omega** — per shell, mean ω and mean α (sum of exponents of primes
> 3) for single-prime vs twin centres. Twin centres carry systematically more
factors, consistent with the shielding proposition.

**table5_enrichment_test** — the core test. Two files are provided:

- `table5_enrichment_test_S8_dualmodel.csv` (used for **Figure 1** and the §3.3
  table): per (shell, ω) stratum, `nodes`, `twins`, `obs_rate_%`, and
  observed/predicted ratios under **both** local models — `obs/pred_old` for the
  naïve ∏ q/(q−2) (tilts 0.93 → 1.17), and `obs/pred_new` for the model
  ∏ (q−1)/(q−3) (sits at 1.000 ± 0.008 across strata spanning a 3.3× enrichment
  range). The decisive panel uses S8 (~10⁸ nodes, the largest strata with
  reliable statistics for every ω = 1..5).
- `table5_enrichment_test_S10.csv` (single-model, full 10¹⁰ scan): the naïve
  model only, retained for completeness from the full-scale run.

The §3.4 Puszkarz-bias recovery (rows (a)/(b)/(c), residual −0.2 %) is produced
separately by `code/puszkarz_S10_correct.py` at the full 10¹⁰ scale.

---

## Citing

See `CITATION.cff`. Please also cite the prior work whose phenomenon this models:
W. Puszkarz, *Statistical bias in the distribution of prime pairs and isolated
primes*, arXiv:1807.00406 (2018).

## License

MIT — see `LICENSE`.
