# Randomized Richardson Trotter — numerical comparison

Reference implementation of the randomized Richardson extrapolation of Trotter
formulas, together with the two figures comparing it against the standard
Trotter baseline:

- **Overhead.** Trotter vs Richardson step counts, one panel per Trotter
  order.
- **Gate depth.** Total gate depth (steps × per-step cost `C_p`), per order
  and as an envelope over `p ∈ {1, 2, 4, 6}`.

## Layout

```
.
├── richardson.py            Core algorithm: Richardson coefficients, grid
│                            search, λ-comm bound, step-count formulas.
├── plotting/
│   ├── common.py            Shared search runner + sidecar writer.
│   ├── common_cli.py        argparse helpers and `--orders` parser.
│   ├── plot_overhead.py     → plots/overhead{,_cropped}.png + sidecars
│   └── plot_gate_depth.py   → plots/gate_depth{,_cropped}.png + sidecars
├── make_plots.py            One-shot reproducer (runs both scripts).
├── plots/                   Generated PNGs + parameter sidecars (committed).
├── requirements.txt
└── README.md
```

## Reproducing the figures

```bash
python -m pip install -r requirements.txt
python make_plots.py
```

After it finishes, `plots/` contains, for each of the two figures:

| File                  | Contents                                                                |
| --------------------- | ----------------------------------------------------------------------- |
| `<fig>.png`           | PNG preview (overhead: panel titles `$p=1,2,4$` only; gate depth: none). |
| `<fig>.pdf`           | Vector PDF companion (preferred for the manuscript).                    |
| `<fig>_cropped.png`   | Caption-friendly overhead (no panel titles).                            |
| `<fig>_cropped.pdf`   | Caption-friendly PDF (overhead + gate depth for `\includegraphics`).    |
| `<fig>.params.md`     | Settings used + per-(p, mode, ε) optimal grids and ‖b‖ norms.           |
| `<fig>.params.json`   | Same data plus the full Richardson coefficient vectors `b_k`.           |

### Including figures in the manuscript

Copy or symlink the PDFs into your paper tree:

```latex
% Three-panel overhead ($p=1,2,4$ panel labels only; no suptitle)
\begin{figure*}
  \centering
  \includegraphics[width=\textwidth]{plots/overhead.pdf}
  \caption{...}
  \label{fig:overhead}
\end{figure*}

% Gate depth (no titles)
\begin{figure*}
  \centering
  \includegraphics[width=\textwidth]{plots/gate_depth.pdf}
  \caption{...}
  \label{fig:gate-depth}
\end{figure*}
```

Use `plots/overhead_cropped.pdf` / `plots/gate_depth_cropped.pdf` instead if
the LaTeX caption supplies all panel labels and you want zero on-figure titles.

For the multi-cap overhead comparison, regenerate with
`python make_plots.py --multi-cap-overhead` and include
`plots/overhead_multi_cap.pdf` (panel labels `$p=1,2,4$` only).

Extra CLI flags can be forwarded to both underlying scripts via `--`:

```bash
python make_plots.py -- --eps-points 30 --q-max 12
```

## Per-plot parameters

Every figure ships with a Markdown sidecar (`<fig>.params.md`) that records the
exact settings used and, for each Trotter order `p` and schedule (`wc` =
well-conditioned grid, `opt` = brute-force optimized grid), a table

| ε | m | q_k | ‖b‖₁ | ‖b‖₁² | ‖b̃‖₁ |
| - | - | --- | ---- | ----- | ----- |

with one row per ε on the search grid. The columns are:

- `m`     — Richardson depth chosen by the search.
- `q_k`   — integer refinement factors (Vandermonde nodes are `s_k = 1/q_k`).
- `‖b‖₁`  — 1-norm of the Richardson coefficients.
- `‖b‖₁²` — sample overhead (colormap variable in the overhead panels).
- `‖b̃‖₁` — `Σ |b_i (q_i/q_min)^e|` (`e=1` WC, `e=p` optimized) in the step bound.

The full Richardson coefficient vector `b_k` for every (p, mode, ε), together
with the resulting Trotter and Richardson step counts, is in the companion
`<fig>.params.json`.

## CLI defaults

Both plot scripts share the same search backend, so the optimal `(m, q_k, b_k)`
triples are identical across the two figures; only the y-axis differs.

| Setting                                          | Default                              | Flag                              |
| ------------------------------------------------ | ------------------------------------ | --------------------------------- |
| Output directory                                 | `plots`                              | `--out-dir`                       |
| Titled output filename                           | `overhead.png` / `gate_depth.png`    | `--output`                        |
| Skip the cropped companion                       | off                                  | `--no-cropped`                    |
| Trotter orders searched                          | `1, 2, 4`                            | `--orders`                        |
| ε grid lower edge (log10)                        | `-6`                                 | `--eps-log-min`                   |
| ε grid upper edge (log10)                        | `log10(0.9) ≈ -0.0458`               | `--eps-log-max`                   |
| Number of ε samples                              | `50`                                 | `--eps-points`                    |
| `q_min, q_max` (Richardson depth `m` also runs `1 … q_max`) | `1, 15`                              | `--q-min`, `--q-max`              |
| `‖b‖₁²` brute-force cap (reject above)           | `10`                                 | `--brute-bnorm-sq-max`            |
| System size `n` (gate-depth only)                | `100`                                | `--n-sys`                         |
| Hide `p = 1` Richardson points with `‖b‖₁² >` …  | `1000`                               | `--omit-p1-sample-overhead-above` |

Run `python plotting/plot_overhead.py --help` (or `plot_gate_depth.py --help`)
for the complete option list.

## Notes on the gate-depth figure

The left panel shows, for each `p ∈ {1, 2, 4}`, the standard Trotter curve
(solid) and the better of the two Richardson schedules (dashed). A closed-form
vanilla Trotter `p = 6` line is overlaid (`C₆ = 2 · 5² = 50`, slope `ε^(-1/6)`);
no Richardson `p = 6` curve is computed.

The right-panel "Best Trotter" envelope is the pointwise minimum over
`p ∈ {1, 2, 4, 6}` (it includes the analytic `p = 6` line). "Best extrapolated"
is the pointwise minimum over `p ∈ {1, 2, 4}` and the two schedules. The faint
colored traces repeat the per-order curves for reference.
