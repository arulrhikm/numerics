# error_analysis/ — the exact-error figure (paper Figure 6)

`empirical_extrapolation_error.ipynb` is the notebook behind Figure 6. It builds the exact and
the Trotterized time-evolution operators for a small Heisenberg chain by dense/sparse linear
algebra, applies the Richardson schedules found in the bounds study, and plots the resulting
operator-norm error. Outputs are saved in the notebook, so the figure can be read without
running anything.

This is a classical simulation of *algorithmic* error. No quantum hardware, no sampling, and no
Hadamard tests are involved.

## What Figure 6 shows (paper, Sec. "Empirical error")

Anisotropic Heisenberg chain with open boundaries on `L = 8` qubits,

    H = J * sum_{i=1}^{L-1} ( X_i X_{i+1} + Y_i Y_{i+1} + 0.9 Z_i Z_{i+1} ) + h * sum_{i=1}^{L} ( Z_i + X_i ) ,   J = 1, h = -1,

time evolution `exp(-iHt)` at `t = 2` approximated by the **second-order** Suzuki-Trotter formula
(`p = 2`), over `5 … 60` Trotter steps.

- **Left panel.** Operator-norm error `|| U_approx - exp(-iHt) ||` against the maximum number of
  Trotter steps in one coherent circuit run, for plain Trotter and for the `p = 2` Richardson
  schedules. Markers are colored by the sample overhead `‖b‖₁²` on a viridis `LogNorm`, the same
  scale the other figures use.
- **Right panel.** The same with a layer of random single-qubit `X`-rotations (angles drawn
  Gaussian, `std = 0.008`) modelling measurement noise. A red dotted line marks the noise rate
  corresponding to one standard deviation. Only the best-conditioned schedules are plotted.

The extrapolated operator is `U_R = sum_k b_k * [P_2(t / (q_k r))]^(q_k r)` for integer `r`,
plotted at `x = max_k q_k r`. The schedules are imported "blind" from the bounds study: they were
optimized for the plane-wave gate-depth bounds, never re-tuned for this Hamiltonian.

## Notebook layout

| Cells | What they do |
|---|---|
| 2–12 | Hamiltonian terms, second-order Trotter operator, `X`-rotation noise layer, exact `expm` |
| 15–23 | Operator-norm error helper; plain-Trotter error vs step count |
| 26 | Extrapolated error, clean (left panel) |
| 30–33 | Noise layer applied, extrapolated error (right panel) |
| 34 | **Figure 6**, both panels side by side → `combined_trotter_plots.pdf` |
| 35–36 | Single-panel variant used for the summary figure → `trotter_error_plot.pdf` |

## Schedules and their provenance

`fig6_schedules.json` holds the 13 schedules the figure plots, extracted verbatim from the
notebook, each with a `provenance` field. Checked against `richardson.py`:

- **All 13 coefficient vectors `b` are reproduced exactly from their `q` grid** (rtol 1e-6) by
  `get_wc_richardson_coefficients` / `get_richardson_coefficients`. Nothing in the figure depends
  on an unverifiable number.
- The five `wc` grids are the LKW closed form at `m = 2, 3, 4, 5, 7`.
- **The eight `opt` grids are not produced by the committed brute-force search.** They appear in no
  committed sidecar, and re-running the search at the paper's stated `q_k ∈ [1, 10]` and at the
  repository default `q_max = 15` returns different grids (e.g. `(7,10)`, `(5,7,12)`, `(4,5,9,15)`
  rather than `(5,8)`, `(2,4,6)`, `(1,3,4,5)`). They therefore come from an earlier search run that
  is not in this repository.

The paper's Fig. 6 caption says these are "the same that were found in the optimization of Figures
4 and 5". Until the original search output is recovered, `fig6_schedules.json` is the record that
makes the figure reproducible.
