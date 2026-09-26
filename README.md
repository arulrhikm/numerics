# Numerics for *Resource-efficient quantum eigenvalue transform with commutator scaling*

Code and figures for

> A. R. Mazumder, J. D. Watson, S. Wang, *Resource-efficient quantum eigenvalue transform with
> commutator scaling*, [arXiv:2608.13862](https://arxiv.org/abs/2608.13862) (2026).

The paper's central primitive is Randomized Extrapolation Trotterization: Richardson
extrapolation in the Trotter step size, compiled by a randomized scheme. This repository holds
the extrapolation search and the scripts behind the paper's numerical figures, nothing else.
Each figure ships with a parameter sidecar recording the exact schedules plotted, so every
number in the paper traces back to a search setting.

## Figures in the paper

| Paper figure | File in `plots/` | Produced by |
|---|---|---|
| Fig. 2(a), gate-depth envelope | `summary.pdf`, left panel | `plotting/plot_summary.py`, from the sidecars below |
| Fig. 2(b), empirical error | not regenerated here, see below | `error_analysis/empirical_extrapolation_error.ipynb` |
| Fig. 4, bounds on Trotter step number | `overhead_multi_cap.pdf` | `plotting/plot_overhead.py --brute-bnorm-sq-caps 10,100,1000` |
| Fig. 5, bounds on gate depth | `gate_depth.pdf` | `plotting/plot_gate_depth.py` |
| Fig. 6, exact error, 8-qubit Heisenberg chain | inline in the notebook | `error_analysis/empirical_extrapolation_error.ipynb` (run all cells) |

**`plots/summary.pdf` is not the published Figure 2.** Its left panel matches Fig. 2(a), but its
right panel is the `p = 2` step-bound panel, whereas the published Fig. 2(b) is the empirical-error
panel from the notebook. See *Reproducibility status* below.

The tag `figures-arxiv-v1` marks the script versions that produced the figures in the arXiv
submission.

## Quick start

```bash
python -m pip install -r requirements.txt
python make_plots.py            # Figs. 2, 4, 5 and their sidecars
```

Figure 6 is a notebook: open `error_analysis/empirical_extrapolation_error.ipynb` and run all
cells (a few minutes; needs `jupyter`, which is not in `requirements.txt`). Its outputs are
committed with the notebook, so the figure can be read without running it.

The brute-force search at the paper's settings (`q_max = 15`, three caps, 50 precisions)
takes tens of minutes and a few gigabytes of memory on a laptop. The search is deterministic:
a rerun with default settings reproduces the committed sidecars byte for byte. For a quick
check of the pipeline, forward a smaller grid to both search scripts (under a minute):

```bash
python make_plots.py -- --eps-points 4 --q-max 4 --out-dir _scratch
python plotting/plot_overhead.py --help       # full option list
```

Run everything from the repository root; the scripts locate `richardson.py` from there.
On Windows use `python`, not `python3`.

## Layout

```
richardson.py            Core algorithm: Richardson coefficients (LKW well-conditioned and
                         optimized Vandermonde), grid search, λ_comm bound, step-count and
                         gate-depth formulas. Everything else imports this.
make_plots.py            One-shot reproducer with the paper's settings.
plotting/
  common.py              Shared search runner, sidecar writer, and
                         load_sidecar_schedules() for downstream studies.
  common_cli.py          argparse helpers shared by the search scripts.
  plot_overhead.py       -> plots/overhead_multi_cap.{png,pdf} + sidecars   (Fig. 4)
  plot_gate_depth.py     -> plots/gate_depth.{png,pdf} + sidecars           (Fig. 5)
  plot_summary.py        -> plots/summary.{png,pdf}, assembled from sidecars (Fig. 2)
plots/                   The three figures and four parameter sidecars (committed).
error_analysis/          Exact-error study (Fig. 6): empirical_extrapolation_error.ipynb,
                         fig6_schedules.json (the schedules it plots, with provenance)
                         and README.md describing the model and the known gaps.
verify_reproduction.py   Compares a regenerated run against the committed sidecars.
requirements.txt, LICENSE (MIT), CITATION.cff
```

## What each figure computes

**Step bounds (Fig. 4).** For each Trotter order `p ∈ {1, 2, 4}` and each target precision
`ε` on a log grid, the search picks the Richardson schedule (depth `m`, integer refinement
factors `q_k`, coefficients `b_k`) minimizing the bound on the maximum Trotter step count.
Two schedules are compared against plain Trotter: the well-conditioned LKW grid (`wc`, even
orders only) and a brute-force optimized grid (`opt`). Marker color is the sample overhead
`‖b‖₁²`; brute-force candidates above a cap are rejected, one curve per cap in
`{10, 100, 1000}`.

**Gate depth (Fig. 5).** The same schedules converted to gate depth, steps × per-step cost
`C_p`, for a plane-wave dual-basis Hamiltonian with `n = 100` basis functions. The left panel
shows each order, with an analytic `p = 6` Trotter line; the right panel is the envelope,
best Trotter over `p ∈ {1, 2, 4, 6}` against best extrapolated over `p ∈ {1, 2, 4}`.

**Summary (Fig. 2).** Left, the gate-depth envelope with two "best extrapolated" curves
(caps 10 and 100); right, the `p = 2` step-bound panel at cap 100. Both panels are read from
the committed sidecars so they match the published figures exactly.

**Exact error (Fig. 6).** Operator-norm error of second-order Trotter and of the extrapolated
operators built from the `p = 2` schedules above, on an 8-qubit anisotropic Heisenberg chain,
with and without a random single-qubit `X`-rotation noise layer. This is a classical simulation
of algorithmic error. It lives in a notebook rather than the `make_plots.py` pipeline;
`error_analysis/README.md` gives the model, the cell map and the known gaps.

## Reproducibility status

Checked against [arXiv:2608.13862](https://arxiv.org/abs/2608.13862) on 2026-09-24.

| Paper figure | Reproduced by this repository? |
|---|---|
| Fig. 4 | Yes, `plotting/plot_overhead.py` output matches the published panels |
| Fig. 5 | Yes, `plotting/plot_gate_depth.py` output matches the published panels |
| Fig. 2(a) | Yes, left panel of `plot_summary.py` |
| Fig. 2(b) | Produced by the notebook, not by `plot_summary.py`; the committed `summary.pdf` pairs Fig. 2(a) with a different right panel |
| Fig. 6 | Yes from `error_analysis/fig6_schedules.json`; the search run that *selected* its brute-force grids is not in the repository |

Two open points, both recorded so a reader does not have to rediscover them:

1. **The brute-force grids behind Fig. 6 are not reproduced by the committed search.** The five `wc`
   schedules are the LKW closed form at `m = 2, 3, 4, 5, 7`, and every coefficient vector `b` in the
   figure is reproduced exactly from its grid by `richardson.py`. The eight `opt` grids, however,
   appear in no committed sidecar and are not returned by the search at either `q_max = 10` or
   `q_max = 15`, so they come from an earlier search run. They are recorded verbatim in
   `error_analysis/fig6_schedules.json`.
2. **`q_max` differs between paper and code.** The paper states the brute-force domain is
   `q_k ∈ [1, 10]`; every committed sidecar records `q_min, q_max = 1, 15` and contains grids that
   reach 15. Figures 4 and 5 as published match the committed `q_max = 15` output.

To check a fresh run against the committed sidecars:

```bash
python make_plots.py -- --out-dir _verify     # tens of minutes, a few GB
python verify_reproduction.py _verify         # or --settings-only for a quick look
```

## Parameter sidecars

Each search writes `<stem>.params.md` and `<stem>.params.json` into `plots/`. Two of the four
are the provenance record of a published figure; the other two are the input the summary
figure is assembled from, and have no figure of their own:

| Search stem | What it is |
|---|---|
| `overhead_multi_cap` | provenance for Fig. 4 |
| `gate_depth` | provenance for Fig. 5 |
| `overhead` | input to the right panel of Fig. 2 |
| `gate_depth_multi_cap` | input to the left panel of Fig. 2 |

The Markdown table lists, per order `p` and mode (`wc`, `opt`), one row per `ε`:

| ε | m | q_k | ‖b‖₁ | ‖b‖₁² | ‖b̃‖₁ |
| - | - | --- | ---- | ----- | ----- |

- `m`, `q_k`: Richardson depth and integer refinement factors (nodes `s_k = 1/q_k`).
- `‖b‖₁`, `‖b‖₁²`: coefficient 1-norm and the sample overhead it implies.
- `‖b̃‖₁`: the suppressed norm `Σ |b_k (q_min/q_k)^e|` entering the step bound.

The JSON carries the same data plus the full coefficient vectors and the resulting Trotter and
Richardson step counts. `plotting.common.load_sidecar_schedules(path, p)` returns the
schedules for one order as plain dicts; downstream studies should read them from there rather
than copying arrays.

## Search settings

| Setting | Default | Flag |
|---|---|---|
| Trotter orders searched | `1, 2, 4` | `--orders` |
| ε grid (log10 edges, points) | `-6`, `log10(0.9)`, `50` | `--eps-log-min`, `--eps-log-max`, `--eps-points` |
| `q_min, q_max` (depth `m` also runs `1 … q_max`) | `1, 15` | `--q-min`, `--q-max` |
| `‖b‖₁²` brute-force caps (Fig. 4) | `10, 100, 1000` | `--brute-bnorm-sq-caps` |
| System size `n` (gate depth only) | `100` | `--n-sys` |
| Output directory | `plots` | `--out-dir` |

## Adding code

Put a new study in its own directory (as `error_analysis/`), import `richardson` and
`plotting.common` rather than re-deriving formulas, read schedules through
`load_sidecar_schedules`, write outputs to `plots/` (or the study's directory) and register
the entry point in `make_plots.py` so one command reproduces every figure.

## License

MIT, see `LICENSE`. Please cite the paper (`CITATION.cff`) if you use this code.
