# Numerics for *Resource-efficient quantum eigenvalue transform with commutator scaling*

Code and figures for

> A. R. Mazumder, J. D. Watson, S. Wang, *Resource-efficient quantum eigenvalue transform with
> commutator scaling*, [arXiv:2608.13862](https://arxiv.org/abs/2608.13862) (2026).

The paper's central primitive is Randomized Extrapolation Trotterization: Richardson
extrapolation in the Trotter step size, compiled by a randomized scheme. This repository holds
the reference implementation of the extrapolation search and every script behind the paper's
numerical figures. Each figure ships with a parameter sidecar recording the exact schedules
plotted, so the numbers in the paper can be traced to a search setting.

## Figures in the paper

| Paper figure | File in `plots/` | Produced by |
|---|---|---|
| Fig. 2, summary of numerical results | `summary.pdf` | `plotting/plot_summary.py` (from the two sidecars below) |
| Fig. 4, bounds on Trotter step number | `overhead_multi_cap.pdf` | `plotting/plot_overhead.py --brute-bnorm-sq-caps 10,100,1000` |
| Fig. 5, bounds on gate depth | `gate_depth.pdf` | `plotting/plot_gate_depth.py` |
| Fig. 6, exact error, 8-qubit Heisenberg chain | `error_analysis/exact-error.pdf` | **pending**, see `error_analysis/README.md` |
| poster panel (Fig. 5 right + Fig. 6) | `exact-error-poster.pdf`, `combined-bounds-and-error.pdf` | `plotting/plot_combined_figure.py` (recomposes the PDFs above) |

The tag `figures-arxiv-v1` marks the script versions that produced the figures in the arXiv
submission and the TQC 2026 poster.

## Quick start

```bash
python -m pip install -r requirements.txt
python make_plots.py                          # Figs. 2, 4 (single cap), 5 + sidecars
python make_plots.py --multi-cap-overhead     # Fig. 4 as printed (three ‖b‖₁² caps)
python make_plots.py --all                    # + sanity-check plot + combined figure
```

The search is deterministic, so a rerun with default settings reproduces the committed
sidecars byte for byte. Extra flags after `--` go to both search scripts:

```bash
python make_plots.py -- --eps-points 30 --q-max 12
python plotting/plot_overhead.py --help       # full option list
```

Run everything from the repository root; the scripts locate `richardson.py` from there.
On Windows use `python`, not `python3`.

## Layout

```
richardson.py            Core algorithm: Richardson coefficients (LKW well-conditioned and
                         optimized Vandermonde), grid search, λ_comm bound, step-count and
                         gate-depth formulas. Everything else imports this.
make_plots.py            One-shot reproducer (see Quick start).
plotting/
  common.py              Shared search runner, sidecar writer, and
                         load_sidecar_schedules() for downstream studies.
  common_cli.py          argparse helpers shared by the search scripts.
  plot_overhead.py       -> plots/overhead{,_cropped,_multi_cap}.{png,pdf} + sidecars
  plot_gate_depth.py     -> plots/gate_depth{,_cropped}.{png,pdf} + sidecars
                            (also gate_depth_multi_cap.params.* for the summary)
  plot_summary.py        -> plots/summary.{png,pdf}, assembled from the sidecars
  plot_sanity_checks.py  -> plots/sanity_check.{png,pdf}, checks on the search objective
  plot_combined_figure.py-> plots/combined-bounds-and-error.pdf, exact-error-poster.pdf
plots/                   Generated figures and their parameter sidecars (committed).
error_analysis/          Exact-error study (paper Fig. 6). README.md states the contract for
                         the script; the notebook is an early n = 4 prototype.
checks/
  verify_tighter_bounds.py        Numerical checks behind the diagonal-restriction bounds.
  check_bnorms_wellconditioned.py WC vs optimized Richardson b-norms on the LKW grid.
requirements.txt, LICENSE (MIT), CITATION.cff
```

## What each figure computes

**Overhead (Fig. 4).** For each Trotter order `p ∈ {1, 2, 4}` and each target precision `ε`
on a log grid, the search picks the Richardson schedule (depth `m`, integer refinement
factors `q_k`, coefficients `b_k`) minimizing the bound on the maximum Trotter step count.
Two schedules are compared against plain Trotter: the well-conditioned LKW grid (`wc`,
even orders only) and a brute-force optimized grid (`opt`). Marker color is the sample
overhead `‖b‖₁²`; brute-force candidates above a cap are rejected, one curve per cap.

**Gate depth (Fig. 5).** The same schedules converted to gate depth, steps × per-step cost
`C_p`, for a plane-wave dual-basis Hamiltonian with `n = 100` basis functions. The left
panel shows each order, with an analytic `p = 6` Trotter line; the right panel is the
envelope, best Trotter over `p ∈ {1, 2, 4, 6}` against best extrapolated over `p ∈ {1, 2, 4}`.

**Summary (Fig. 2).** Left, the gate-depth envelope with two "best extrapolated" curves
(caps 10 and 100); right, the `p = 2` overhead panel. Both panels are read from the committed
sidecars so they match the published figures exactly.

**Exact error (Fig. 6).** Operator-norm error of second-order Trotter and of the extrapolated
operators built from the `p = 2` schedules above, on an 8-qubit anisotropic Heisenberg chain,
with and without a random single-qubit `X`-rotation noise layer. The script is being added;
`error_analysis/README.md` specifies inputs, model and output.

## Parameter sidecars

Every search-driven figure writes `<fig>.params.md` and `<fig>.params.json` beside it. The
Markdown table lists, per order `p` and mode (`wc`, `opt`), one row per `ε`:

| ε | m | q_k | ‖b‖₁ | ‖b‖₁² | ‖b̃‖₁ |
| - | - | --- | ---- | ----- | ----- |

- `m`, `q_k`: Richardson depth and integer refinement factors (nodes `s_k = 1/q_k`).
- `‖b‖₁`, `‖b‖₁²`: coefficient 1-norm and the sample overhead it implies.
- `‖b̃‖₁`: the suppressed norm `Σ |b_k (q_min/q_k)^e|` entering the step bound.

The JSON carries the same data plus the full coefficient vectors and the resulting Trotter and
Richardson step counts. `plotting.common.load_sidecar_schedules(path, p)` returns the
schedules for one order as plain dicts; downstream studies should read them from there
rather than copying arrays.

## Search settings

| Setting | Default | Flag |
|---|---|---|
| Trotter orders searched | `1, 2, 4` | `--orders` |
| ε grid (log10 edges, points) | `-6`, `log10(0.9)`, `50` | `--eps-log-min`, `--eps-log-max`, `--eps-points` |
| `q_min, q_max` (depth `m` also runs `1 … q_max`) | `1, 15` | `--q-min`, `--q-max` |
| `‖b‖₁²` brute-force cap | `100` | `--brute-bnorm-sq-max` |
| System size `n` (gate depth only) | `100` | `--n-sys` |
| Output directory | `plots` | `--out-dir` |

## Adding code

Put a new study in its own directory (as `error_analysis/`), import `richardson` and
`plotting.common` rather than re-deriving formulas, read schedules through
`load_sidecar_schedules`, write outputs to `plots/` (or the study's directory) and register
the entry point in `make_plots.py` so `--all` reproduces everything. Keep the paper's figure
numbering in the table above.

## License

MIT, see `LICENSE`. Please cite the paper (`CITATION.cff`) if you use this code.
