"""Overhead figure: Trotter vs two Richardson schedules, one panel per order.

Each panel shows, for a given Trotter order ``p``:

* the standard Trotter baseline (gray solid line);
* Richardson step counts under the well-conditioned grid (square markers);
* Richardson step counts under the brute-force optimized grid (circle markers).

Marker color encodes the sample overhead ``‖b‖₁²``. Writes a titled PNG plus a
caption-friendly ``_cropped`` companion, and a Markdown + JSON sidecar with the
exact settings and per-(p, mode, ε) grids.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import richardson as rt
from plotting.common_cli import (
    add_shared_grid_args,
    add_shared_search_args,
    parse_orders,
    resolve_output_dir,
)
from plotting import common as cm

FS_AXIS = 15
FS_LEGEND = 12
FS_SUPTITLE = 17
FS_CBAR = 15
LEGEND_MARKERSIZE = 10
LEGEND_LINEWIDTH = 2.25
LEGEND_HANDLELEN = 2.75
CBAR_VMIN_FLOOR = 1.0
CBAR_VMAX = 1000.0
# Per-order y-axis ranges on log scale.
OVERHEAD_YLIM_BY_P = {1: (1.0, 1e6), 2: (1.0, 1e3), 4: (1.0, 5e1)}


def _style_ax(ax, title=None, ylabel="Number of steps"):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Precision  $\varepsilon$", fontsize=FS_AXIS)
    ax.set_ylabel(ylabel, fontsize=FS_AXIS)
    if title is not None:
        ax.set_title(title, fontsize=FS_AXIS, pad=8)
    ax.tick_params(axis="both", labelsize=FS_AXIS - 2)
    ax.grid(True, which="major", ls="-", alpha=0.25)
    ax.grid(True, which="minor", ls=":", alpha=0.12)


def _panel_title(p: int, b2_cap: float, c_pref: float, brute_permutations: bool) -> str:
    cap_s = str(int(round(b2_cap))) if abs(b2_cap - round(b2_cap)) < 1e-9 else f"{b2_cap:g}"
    line1 = rf"Order $p = {p}$ ($\|\mathbf{{b}}\|_1^2 \leq {cap_s}$)"
    line2 = rf"static prefactor ($c = {c_pref:.3f}$)"
    if brute_permutations:
        line2 = f"{line2}, permuted $q$"
    return f"{line1}\n{line2}"


def parse_args():
    parser = argparse.ArgumentParser(description="Overhead figure (one panel per Trotter order).")
    parser.add_argument("--out-dir", type=str, default="plots", help="Output directory.")
    parser.add_argument("--output", type=str, default="overhead.png", help="Titled output filename.")
    add_shared_grid_args(parser, default_q_max=10)
    parser.set_defaults(eps_log_max=np.log10(0.9))
    parser.add_argument("--orders", type=str, default="1,2,4", help="Comma-separated Trotter orders p.")
    add_shared_search_args(parser)
    parser.add_argument(
        "--omit-p1-sample-overhead-above",
        type=float,
        default=1000.0,
        help="Hide p=1 Richardson points with ||b||_1^2 above this value.",
    )
    parser.add_argument("--cmap", type=str, default="viridis", help="Colormap for ||b||_1^2.")
    parser.add_argument(
        "--no-cropped",
        action="store_true",
        help="Skip the caption-friendly (no-title) companion file.",
    )
    return parser.parse_args()


def _build_figure(
    *,
    orders: list[int],
    errors: np.ndarray,
    results_by_order: dict,
    args: argparse.Namespace,
    b2_cap: float,
    omit_p1: float,
    omit_titles: bool,
) -> plt.Figure:
    norm = mcolors.LogNorm(vmin=CBAR_VMIN_FLOOR, vmax=CBAR_VMAX, clip=False)
    cmap = plt.get_cmap(args.cmap)

    n_p = len(orders)
    fig_w = max(6.5 * n_p, 6.5)
    fig_h = 5.5 if omit_titles else 6.5
    fig, axes = plt.subplots(1, n_p, figsize=(fig_w, fig_h), constrained_layout=True)
    fig.set_constrained_layout_pads(h_pad=0.06, w_pad=0.04, hspace=0.02, wspace=0.02)
    if n_p == 1:
        axes = np.array([axes])

    last_mappable = None
    for idx, p in enumerate(orders):
        ax = axes[idx]
        r_p = results_by_order[p]
        trotter_baseline = r_p.get("wc", r_p["opt"])["trotter"]
        ax.plot(
            errors,
            trotter_baseline,
            "-",
            color="gray",
            linewidth=2,
            alpha=0.35,
            zorder=1,
        )

        for mode, marker in (("wc", "s"), ("opt", "o")):
            if mode not in r_p:
                continue
            steps = np.asarray(r_p[mode]["richardson"], dtype=float)
            so = np.asarray(r_p[mode]["sample_overhead"], dtype=float)
            visible = cm.visible_mask_for_sample_overhead(
                order=int(p), sample_overhead=so, omit_p1_sample_overhead_above=omit_p1
            )
            sc = ax.scatter(
                errors[visible],
                steps[visible],
                c=so[visible],
                cmap=cmap,
                norm=norm,
                s=70,
                marker=marker,
                linewidths=0,
                edgecolors="none",
                zorder=3 if mode == "wc" else 4,
            )
            if mode == "opt":
                last_mappable = sc

        if not omit_titles:
            c_pref = rt.richardson_b_over_eps_prefactor(p)
            _style_ax(
                ax,
                title=_panel_title(int(p), b2_cap, float(c_pref), bool(args.brute_permutations)),
            )
        else:
            _style_ax(ax)
        ylim = OVERHEAD_YLIM_BY_P.get(int(p))
        if ylim is not None:
            ax.set_ylim(*ylim)

        legend_handles = [
            Line2D([0], [0], color="gray", linewidth=LEGEND_LINEWIDTH, alpha=0.45, label="Trotter"),
            Line2D([0], [0], marker="s", linestyle="None", markersize=LEGEND_MARKERSIZE * 0.7,
                   markerfacecolor="#888", markeredgecolor="none", label="Richardson (well-cond.)"),
            Line2D([0], [0], marker="o", linestyle="None", markersize=LEGEND_MARKERSIZE * 0.7,
                   markerfacecolor="#444", markeredgecolor="none", label="Richardson (optimized)"),
        ]
        ax.legend(
            handles=legend_handles,
            fontsize=FS_LEGEND,
            loc="upper left",
            framealpha=0.9,
            handlelength=LEGEND_HANDLELEN,
            handletextpad=0.65,
            borderpad=0.45,
            labelspacing=0.55,
        )

    if not omit_titles:
        fig.suptitle("Overhead", fontsize=FS_SUPTITLE, fontweight="600", y=1.03)

    cbar = fig.colorbar(
        last_mappable, ax=axes, orientation="vertical", fraction=0.02, pad=0.02, aspect=30
    )
    cbar.set_label(r"$\|\mathbf{b}\|_1^2$  (Richardson coefficients)", fontsize=FS_CBAR)
    cbar.ax.tick_params(labelsize=FS_AXIS - 2)
    return fig


def main() -> None:
    args = parse_args()
    rt.set_lambda_scale_mode(args.lambda_mode)
    plt.rcParams["pdf.fonttype"] = 42

    output_dir = resolve_output_dir(_ROOT, args.out_dir)
    orders = parse_orders(args.orders)
    errors = np.logspace(args.eps_log_max, args.eps_log_min, num=int(args.eps_points))
    b2_cap = float(args.brute_bnorm_sq_max)
    omit_p1 = float(args.omit_p1_sample_overhead_above)

    print(
        f"Overhead: eps in [{10**args.eps_log_min:g}, {10**args.eps_log_max:g}] "
        f"({args.eps_points} pts), q in [{args.q_min},{args.q_max}], "
        f"||b||_1^2 cap={b2_cap:g}, orders={orders}, lambda-mode={args.lambda_mode}"
    )

    results_by_order = cm.compute_results(
        errors=errors,
        orders=orders,
        q_max=int(args.q_max),
        q_min=int(args.q_min),
        b2_cap=b2_cap,
        brute_permutations=bool(args.brute_permutations),
        brute_permutations_max_count=int(args.brute_permutations_max_count),
        step_mode="plane_wave",
        n_sys=1,
    )

    def save(omit_titles: bool, name: str) -> None:
        fig = _build_figure(
            orders=orders, errors=errors, results_by_order=results_by_order,
            args=args, b2_cap=b2_cap, omit_p1=omit_p1, omit_titles=omit_titles,
        )
        path = output_dir / name
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {path}")

    stem = Path(args.output)
    save(omit_titles=False, name=args.output)
    if not args.no_cropped:
        save(omit_titles=True, name=f"{stem.stem}_cropped{stem.suffix}")

    md_path, json_path = cm.write_params_sidecar(
        output_dir=output_dir,
        stem=stem.stem,
        title="Overhead figure parameters",
        settings={
            "script": "plotting/plot_overhead.py",
            "step_mode": "plane_wave",
            "λ-mode": args.lambda_mode,
            "orders (p)": ", ".join(str(p) for p in orders),
            "ε grid": (
                f"{args.eps_points} points, log10 in "
                f"[{args.eps_log_min:.6g}, {args.eps_log_max:.6g}]"
            ),
            "q_min, q_max": f"{args.q_min}, {args.q_max}",
            "‖b‖₁² search cap": f"{b2_cap:g}",
            "brute permutations": str(bool(args.brute_permutations)),
            "drop p=1 points with ‖b‖₁² >": f"{omit_p1:g}",
            "system size n (λ-comm)": "1 (overhead omits the explicit n prefactor)",
        },
        errors=errors,
        orders=orders,
        results_by_order=results_by_order,
    )
    print(f"  Saved {md_path}")
    print(f"  Saved {json_path}")


if __name__ == "__main__":
    main()
