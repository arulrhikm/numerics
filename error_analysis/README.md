# error_analysis/ — the exact-error figure (paper Figure 6)

**Status: the script that produced Figure 6 is not in the repository yet.** It is being
contributed by Samson Wang. This file records what that script must compute and where it
must write its output, so that it drops into the existing pipeline without further changes.

`Exact_extrapolation_error.ipynb` is the earlier **prototype**: an isotropic chain at `n = 4`,
`J = 1`, `h = 0.5`, with the extrapolation schedules pasted in by hand. It is kept for
reference only and does not reproduce the paper's figure.

## What Figure 6 shows (paper, Sec. "Empirical error")

Anisotropic Heisenberg chain with open boundaries on `L = 8` qubits,

    H = J * sum_{i=1}^{L-1} ( X_i X_{i+1} + Y_i Y_{i+1} + 0.9 Z_i Z_{i+1} ) - h * sum_{i=1}^{L} Z_i ,   J = h = 1,

time evolution `exp(-iHt)` approximated by the **second-order** Suzuki-Trotter formula
(`p = 2`).

- **Left panel.** Operator-norm error `|| U_approx - exp(-iHt) ||` against the maximum
  number of Trotter steps in one coherent circuit run, for plain Trotter and for every
  Richardson schedule the overhead search found at `p = 2`. Markers are colored by the
  sample overhead `‖b‖₁²` on the same viridis `LogNorm` the other figures use.
- **Right panel.** The same with a layer of random single-qubit `X`-rotations inserted to
  model measurement noise; a red dotted line marks the noise rate corresponding to one
  standard deviation. Only the best-conditioned schedules, `‖b‖₁ ≤ √10`, are plotted.

The paper does **not** state the noise rate or the seed; both must come with the code.

## Interface the script should follow

1. **Schedules come from the committed sidecar, not from hard-coded arrays.**

   ```python
   from plotting.common import load_sidecar_schedules
   schedules = load_sidecar_schedules("plots/overhead.params.json", p=2)
   # each entry: {"mode": "wc"|"opt", "eps": float, "q_grid": [int, ...],
   #              "b": [float, ...], "bnorm_sq": float}
   ```

   This is what the paper means by importing the schedules "blind" from the bounds
   study: the extrapolated operator is `U_R = sum_k b_k * [P_2(t / (q_k r))]^(q_k r)` for
   integer `r`, plotted at `x = max_k q_k r`.

2. **Output.** Write `error_analysis/exact-error.pdf` (and `.png`), two panels side by
   side, clean on the left and noisy on the right. `plotting/plot_combined_figure.py`
   lifts the panels out of that PDF by position, so keep that order.

3. **Entry point.** A single `python error_analysis/exact_error.py` from the repo root,
   numpy/scipy/matplotlib only, with the model parameters, `t`, the Trotter-step window,
   the noise rate and the seed as module-level constants or CLI flags. Once it exists,
   add it to the `--all` list in `make_plots.py`.
