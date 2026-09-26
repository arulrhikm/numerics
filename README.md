# Numerics for "Resource-efficient quantum eigenvalue transform with commutator scaling"

Code and figures for

> A. R. Mazumder, J. D. Watson, S. Wang, *Resource-efficient quantum eigenvalue transform with
> commutator scaling*, [arXiv:2608.13862](https://arxiv.org/abs/2608.13862) (2026).

The paper's method, Randomized Extrapolation Trotterization, applies Richardson extrapolation in
the Trotter step size and compiles it with a randomized scheme. This repository has the
extrapolation search and the scripts that make the paper's numerical figures. Each figure comes
with a parameter sidecar listing the schedules it plots, so the numbers in the paper can be traced
to a search setting.

## Figures in the paper

| Paper figure | File in `plots/` | Produced by |
|---|---|---|
| Fig. 2(a), gate-depth envelope | `summary.pdf`, left panel | `plotting/plot_summary.py`, from the sidecars below |
| Fig. 2(b), empirical error | not regenerated here, see below | `error_analysis/empirical_extrapolation_error.ipynb` |
| Fig. 4, bounds on Trotter step number | `overhead_multi_cap.pdf` | `plotting/plot_overhead.py --brute-bnorm-sq-caps 10,100,1000` |
| Fig. 5, bounds on gate depth | `gate_depth.pdf` | `plotting/plot_gate_depth.py` |
| Fig. 6, exact error, 8-qubit Heisenberg chain | inline in the notebook | `error_analysis/empirical_extrapolation_error.ipynb` (run all cells) |

Note that `plots/summary.pdf` is not the published Figure 2. Its left panel is Fig. 2(a), but its
right panel is the `p = 2` step-bound panel. The published Fig. 2(b) is the empirical-error panel
from the notebook.

The tag `figures-arxiv-v1` marks the script versions used for the arXiv v1 figures. The figures in
`plots/` have since been corrected; see [Reproducibility status](#reproducibility-status).

## Quick start

```bash
python -m pip install -r requirements.txt
python make_plots.py            # Figs. 2, 4, 5 and their sidecars
```

Figure 6 comes from a notebook. Open `error_analysis/empirical_extrapolation_error.ipynb` and run
all cells. This takes a few minutes and needs `jupyter`, which is not in `requirements.txt`. The
outputs are saved in the notebook, so you can see the figure without running it.

The brute-force search at the paper's settings (`q_max = 10`, three caps, 50 precisions) takes a few
minutes on a laptop (v1's `q_max = 15` took tens of minutes and a few gigabytes). The search is deterministic, and a rerun with
the default settings reproduces the committed sidecars byte for byte. To test the pipeline quickly,
pass a smaller grid to both search scripts (under a minute):

```bash
python make_plots.py -- --eps-points 4 --q-max 4 --out-dir _scratch
python plotting/plot_overhead.py --help       # full option list
```

Run everything from the repository root, since the scripts look for `richardson.py` there. On
Windows use `python`, not `python3`.

## Layout

```
richardson.py            Richardson coefficients (LKW well-conditioned and optimized
                         Vandermonde), grid search, λ_comm bound, step-count and gate-depth
                         formulas. The other scripts import it.
make_plots.py            Regenerates Figs. 2, 4 and 5 with the paper's settings.
plotting/
  common.py              Search runner, sidecar writer, and load_sidecar_schedules().
  common_cli.py          argparse helpers for the search scripts.
  plot_overhead.py       -> plots/overhead_multi_cap.{png,pdf} + sidecars   (Fig. 4)
  plot_gate_depth.py     -> plots/gate_depth.{png,pdf} + sidecars           (Fig. 5)
  plot_summary.py        -> plots/summary.{png,pdf}, built from sidecars    (Fig. 2)
plots/                   The three figures and four parameter sidecars.
error_analysis/          Exact-error study (Fig. 6): the notebook, fig6_schedules.json
                         (the schedules it plots, with where each came from) and a README.
verify_reproduction.py   Compares a fresh run against the committed sidecars.
requirements.txt, LICENSE (MIT), CITATION.cff
```

## What each figure computes

Fig. 4 bounds the Trotter step count. For each order `p ∈ {1, 2, 4}` and each precision `ε` on a
log grid, the search picks the Richardson schedule (depth `m`, integer refinement factors `q_k`,
coefficients `b_k`) with the smallest bound on the maximum number of Trotter steps. It compares
plain Trotter with two kinds of schedule: the well-conditioned LKW grid (`wc`, even orders only)
and a brute-force optimized grid (`opt`). Marker color is the sample overhead `‖b‖₁²`. Brute-force
candidates above a cap are rejected, and there is one curve for each cap in `{10, 100, 1000}`.

Fig. 5 converts the same schedules to gate depth, meaning steps times the per-step cost `C_p`, for
a plane-wave dual-basis Hamiltonian with `n = 100` basis functions, using brute-force schedules
with `‖b‖₁² ≤ 10`. The left panel shows each
order, with an analytic `p = 6` Trotter line. The right panel compares the best Trotter formula
over `p ∈ {1, 2, 4, 6}` with the best extrapolated schedule over `p ∈ {1, 2, 4}`.

Fig. 2(a) is the gate-depth envelope with two "best extrapolated" curves, for caps 10 and 100.
`plot_summary.py` reads it from the committed sidecars.

Fig. 6 shows the operator-norm error of second-order Trotter and of the extrapolated operators
built from the `p = 2` schedules, on an 8-qubit anisotropic Heisenberg chain, with and without a
layer of random single-qubit `X` rotations as noise. It is a classical simulation of the
algorithmic error, done in a notebook outside the `make_plots.py` pipeline.
`error_analysis/README.md` describes the model and the notebook cells.

## Reproducibility status

As of 2026-09-26, the figures in `plots/` differ from the ones published in
[arXiv:2608.13862v1](https://arxiv.org/abs/2608.13862). Three bugs were fixed so that the code
matches what the paper says it computes. The tag `figures-arxiv-v1` keeps the exact code and
figures behind v1.

| Paper figure | This repository |
|---|---|
| Fig. 4 | `plotting/plot_overhead.py`, corrected (items 1, 2 below) |
| Fig. 5 | `plotting/plot_gate_depth.py`, corrected (items 1–3 below) |
| Fig. 2(a) | left panel of `plot_summary.py`, corrected (items 1, 2 below) |
| Fig. 2(b) | made by the notebook, not `plot_summary.py`; the committed `summary.pdf` has a different right panel |
| Fig. 6 | from `error_analysis/fig6_schedules.json`; the search run that chose its brute-force grids is not in the repository |

Fixed since v1. None of these moves the ε ≈ 10⁻² crossover in Figs. 2(a) and 5 (cap 10): the last
grid point where extrapolation wins is still ε = 1.02 × 10⁻².

1. **λ_comm constant at `p = 1` and `p = 4`.** `LEMMA57_GEOMETRIC_RATIO_BY_P` held the ratios at
   `p = 1, 4` already raised to `1 + 1/p` (paper Eq. 254), and the step formulas raised them again.
   So v1 used 2.2605 and 1.0559 where the paper states 1.5035 and 1.0445. The table now holds the
   bare Lemma 52 ratios `{1: 1.2262, 2: 1.0968, 4: 1.0354}`, equal to
   `richardson.lemma57_geometric_ratio(p)`. As a result the `p = 1` extrapolated curves are 1.50×
   lower than in v1, and the `p = 4` curves are 1.1% lower. `p = 2` is unchanged.
2. **Brute-force domain.** The paper states `q_k ∈ [1, 10]`, but v1 searched `[1, 15]`. The default
   `--q-max` is now 10. Relative to v1, this raises the brute-force curves by up to 1.68× at
   `p = 1` and by up to 10% at `p = 2, 4`. The LKW grids are not bounded by `q_max`; they reach
   `q = 83`.
3. **Fig. 5 sample-overhead cap.** The paper states that Fig. 5 uses schedules with `‖b‖₁² ≤ 10`,
   but v1 used 100. `plot_gate_depth.py` now defaults to cap 10. The Fig. 5 crossover moves from
   ε ≈ 1.79 × 10⁻² (v1) to 1.02 × 10⁻².

Remaining differences between the paper text and the code:

4. **Fig. 6 grids.** The committed search does not reproduce the brute-force grids in Fig. 6. The
   five `wc` schedules are the LKW closed form at `m = 2, 3, 4, 5, 7`, and `richardson.py`
   reproduces every coefficient vector `b` in the figure exactly from its grid. But the eight `opt`
   grids are in no committed sidecar, and the search does not return them at `q_max = 10` or
   `q_max = 15`, so they come from an earlier run. `error_analysis/fig6_schedules.json` records them
   as used.
5. **Fig. 2(a) caption.** It says "factor 10", but the panel draws two envelopes, for caps 10 and
   100.
6. **Gate depth.** The paper's gate-depth formula (Eqs. 239–240) carries the stage count as
   `Υ^{2+1/p}`, but `compute_steps_gate_depth` uses `C_p^{1+1/p}`, the step-count scaling, and
   Figs. 2(a) and 5 follow the code: the `p = 6` Trotter line is ≈ 80 at ε = 1. Multiplying every
   curve by its `C_p` (1, 2, 10, 50) would leave all crossovers unchanged but visibly move the
   curves.
7. **Prefactors.** `RICHARDSON_K_BY_P = {1: 2.7232, 2: 3.627, 4: 5.15485}` is fixed per order,
   whereas the paper describes `a(ε)` as refined for each target precision. The derivation of these
   values is not in this repository.

Also note that `--n-sys` has no effect on any plotted number: `n` cancels in every ratio, and
`compute_lambda_scale` ignores it. "`n = 100` basis functions" labels the setting but does not
change the curves.

To check a fresh run against the committed sidecars:

```bash
python make_plots.py -- --out-dir _verify     # a few minutes
python verify_reproduction.py _verify         # or --settings-only for a quick look
```

## Parameter sidecars

Each search writes `<stem>.params.md` and `<stem>.params.json` to `plots/`. Two of the four record
the settings behind a published figure. The other two are inputs to the summary figure and have
no figure of their own.

| Search stem | Used for |
|---|---|
| `overhead_multi_cap` | Fig. 4 |
| `gate_depth` | Fig. 5 |
| `overhead` | input to the right panel of Fig. 2 |
| `gate_depth_multi_cap` | input to the left panel of Fig. 2 |

The Markdown file has one table per order `p` and mode (`wc`, `opt`), with one row per `ε`:

| ε | m | q_k | ‖b‖₁ | ‖b‖₁² | ‖b̃‖₁ |
| - | - | --- | ---- | ----- | ----- |

- `m`, `q_k`: Richardson depth and integer refinement factors (nodes `s_k = 1/q_k`).
- `‖b‖₁`, `‖b‖₁²`: the 1-norm of the coefficients and the sample overhead it gives.
- `‖b̃‖₁`: the suppressed norm `Σ |b_k (q_min/q_k)^e|` that enters the step bound.

The JSON file has the same data plus the full coefficient vectors and the resulting Trotter and
Richardson step counts. `plotting.common.load_sidecar_schedules(path, p)` returns the schedules for
one order as plain dicts. New studies should load schedules this way instead of copying arrays.

## Search settings

| Setting | Default | Flag |
|---|---|---|
| Trotter orders searched | `1, 2, 4` | `--orders` |
| ε grid (log10 edges, points) | `-6`, `log10(0.9)`, `50` | `--eps-log-min`, `--eps-log-max`, `--eps-points` |
| `q_min, q_max` (depth `m` also runs `1 … q_max`) | `1, 10` | `--q-min`, `--q-max` |
| `‖b‖₁²` brute-force caps (Fig. 4) | `10, 100, 1000` | `--brute-bnorm-sq-caps` |
| `‖b‖₁²` cap (Fig. 5, `plot_gate_depth.py`) | `10` | `--brute-bnorm-sq-max` |
| System size `n` (gate depth only) | `100` | `--n-sys` |
| Output directory | `plots` | `--out-dir` |

## Adding code

Put a new study in its own directory, as `error_analysis/` does. Import `richardson` and
`plotting.common` rather than rewriting the formulas, and load schedules with
`load_sidecar_schedules`. Write outputs to `plots/` or the study's own directory, and add the entry
point to `make_plots.py` so that one command regenerates every figure.

## License

MIT, see `LICENSE`. If you use this code, please cite the paper (`CITATION.cff`).
