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
├── verify_tighter_bounds.py Checks the diagonal-restriction derivation
│                            (`notes/derivation-tighter-bounds.tex`): the
│                            `R_p` / `R_p^{1+1/p}` ratios, pinching contraction,
│                            `p_0(ε)`, and the row-(b) collapse condition.
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
- `‖b̃‖₁` — `Σ |b_i (q_min/q_i)^e|` (`e=1` WC, `e=p` optimized) in the step bound.

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
| `‖b‖₁²` brute-force cap (reject above)           | `100`                                | `--brute-bnorm-sq-max`            |
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

## Summary figure

`plotting/plot_summary.py` assembles a side-by-side, two-square-panel
`summary.{png,pdf}` from committed sidecars (so it matches the published panels
exactly — regenerate those first):

- **Left**: the gate-depth envelope — gray "Best Trotter" plus two "Best
  extrapolated" envelopes (‖b‖₁² search caps 10 and 100). Each point is a
  filled triangle colored by sample overhead `‖b‖₁²`; cap budgets are read off
  the shared colorbar (not open vs filled markers). Data from
  `plots/gate_depth_multi_cap.params.json`.
- **Right**: the `p = 2` overhead panel — Trotter baseline plus LKW (squares)
  and brute-force (triangles) schedules, also colored by `‖b‖₁²`.

Both panels are square so the figure can sit beside the exact-error plot.
`plot_gate_depth.py` always emits the `gate_depth_multi_cap` sidecar (no PNG)
alongside the single-cap gate-depth figure, and `make_plots.py` runs the summary
last, after the overhead and gate-depth sidecars are written.

## Combined bounds-and-error figure

`plotting/plot_combined_figure.py` builds one figure holding paper Figure 5's
right panel plus both panels of Figure 6:

    (a) bound on gate depth, best curve across orders   <- plots/gate_depth.pdf
    (b) exact operator-norm error, 8-qubit Heisenberg    <- images/exact-error.pdf
    (c) the same under measurement noise                 <- images/exact-error.pdf

It is a **recomposition, not a rerun**: the panels are lifted as vector art out
of the PDFs the manuscript already uses, so the curves cannot drift from the
paper. It does not read the sidecars.

    # manuscript figure -> plots/combined-bounds-and-error.pdf
    python plotting/plot_combined_figure.py

    # poster block 7 -> plots/exact-error-poster.pdf
    python plotting/plot_combined_figure.py --layout row --panels bc \
           --width 300 --out exact-error-poster

Two deliverables, different constraints. `poster.tex` block 7 uses the second
one; the first is for Overleaf. `--out` keeps them from overwriting each other.

**Required input, and the reason this script exists at all:** Figure 6 cannot be
regenerated from this repo. `error_analysis/Exact_extrapolation_error.ipynb`
runs at n = 4 where the paper is n = 8 and has no `savefig`, and the
random-Pauli-X noise model behind panel (c) exists nowhere in code — only as
prose in `paper.tex`, which never states the noise rate. So `images/exact-error.pdf`
has to be copied in from the Overleaf project; the script exits with instructions
if it is missing.

Measured type, since it is the thing that decides the layout. At `\textwidth`
(165mm), for the manuscript:

| panels  | layout | type on the page |
|---------|--------|------------------|
| a, b, c | grid   | 6.3-7.3pt, 144mm tall — **default** |
| a, c    | row    | 5.6-6.5pt |
| a, b, c | row    | 3.9-4.6pt — unreadable |

Three across does not fit a single-column page; `--layout grid` stacks (a) over
(b)(c) and pays in height instead. At the poster's 300mm slot:

| panels  | type on the board | height | slot is 125mm |
|---------|-------------------|--------|---------------|
| b, c    | 11.5-13.3pt | 129mm | **poster default** |
| a, b    | 10.9-12.6pt | 118mm | |
| a, b, c | 7.1-8.3pt   |  80mm | unusable |

None of these clear the 16pt bar the rest of the poster holds, and none can:
the panels are lifted art and cannot be re-typeset, so enlarging them only
overflows the board. The poster therefore takes the two exact-error panels —
the strongest result — and demotes the gate-depth bound to a caption line.

Panel letters are assigned consecutively over whatever is kept, so dropping a
middle panel gives "(a) (b)", not "(a) (c)".

The script prints the placed type size for every panel and flags a spread over
1.5pt, since the panels come from sources with different native font sizes
(0.103 vs 0.120 pt per mm of plot-frame height).
