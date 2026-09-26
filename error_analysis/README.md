# error_analysis: the exact-error figure (paper Figure 6)

`empirical_extrapolation_error.ipynb` makes Figure 6. It builds the exact and Trotterized
time-evolution operators for a small Heisenberg chain with dense and sparse linear algebra, applies
the Richardson schedules from the bounds study, and plots the operator-norm error. The outputs are
saved in the notebook, so you can see the figure without running anything.

The notebook simulates algorithmic error classically. It involves no quantum hardware, no
sampling and no Hadamard tests.

## What Figure 6 shows (paper, Sec. "Empirical error")

The model is an anisotropic Heisenberg chain with open boundaries on `L = 8` qubits,

    H = J * sum_{i=1}^{L-1} ( X_i X_{i+1} + Y_i Y_{i+1} + 0.9 Z_i Z_{i+1} ) + h * sum_{i=1}^{L} ( Z_i + X_i ) ,   J = 1, h = -1,

with `exp(-iHt)` at `t = 2` approximated by the second-order Suzuki-Trotter formula (`p = 2`),
using 5 to 60 Trotter steps.

- Left panel: the operator-norm error `|| U_approx - exp(-iHt) ||` against the maximum number of
  Trotter steps in one circuit, for plain Trotter and for the `p = 2` Richardson schedules. Markers
  are colored by the sample overhead `‖b‖₁²` on a viridis `LogNorm`, the same scale as the other
  figures.
- Right panel: the same, with a layer of random single-qubit `X` rotations (Gaussian angles,
  `std = 0.008`) standing in for measurement noise. A red dotted line marks the noise rate at one
  standard deviation. Only the best-conditioned schedules are plotted.

The extrapolated operator is `U_R = sum_k b_k * [P_2(t / (q_k r))]^(q_k r)` for integer `r`,
plotted at `x = max_k q_k r`. The schedules were taken as-is from the bounds study. They were
optimized for the plane-wave gate-depth bounds and not re-tuned for this Hamiltonian.

## Notebook layout

| Cells | What they do |
|---|---|
| 2 to 12 | Hamiltonian terms, second-order Trotter operator, `X`-rotation noise layer, exact `expm` |
| 15 to 23 | Operator-norm error helper; plain-Trotter error against step count |
| 26 | Extrapolated error without noise (left panel) |
| 30 to 33 | Noise layer applied, extrapolated error (right panel) |
| 34 | Figure 6, both panels side by side, saved as `combined_trotter_plots.pdf` |
| 35 and 36 | Single-panel version used in the summary figure, saved as `trotter_error_plot.pdf` |

## Where the schedules come from

`fig6_schedules.json` has the 13 schedules the figure plots, copied exactly from the notebook, each
with a `provenance` field. Checked against `richardson.py`:

- `get_wc_richardson_coefficients` and `get_richardson_coefficients` reproduce all 13 coefficient
  vectors `b` from their `q` grids (rtol 1e-6), so every number in the figure can be checked.
- The five `wc` grids are the LKW closed form at `m = 2, 3, 4, 5, 7`.
- The eight `opt` grids do not come from the committed brute-force search. They are in no committed
  sidecar, and rerunning the search with the paper's `q_k ∈ [1, 10]` or the repository default
  `q_max = 15` gives different grids, for example `(7,10)`, `(5,7,12)` and `(4,5,9,15)` in place of
  `(5,8)`, `(2,4,6)` and `(1,3,4,5)`. They must come from an earlier search run that is not in this
  repository.

The caption of Fig. 6 in the paper says these are "the same that were found in the optimization of
Figures 4 and 5". Until the output of that search turns up, `fig6_schedules.json` is what makes the
figure reproducible.
