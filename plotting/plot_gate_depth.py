"""Gate-depth figure: per-order curves (left) and order-envelope (right).

Same ``compute_min_samples`` schedules as the overhead figure, but with
``compute_steps_gate_depth`` (``n = 100`` by default) on the y-axis. The left
panel additionally overlays a closed-form vanilla Trotter ``p = 6`` line (no
Richardson ``p = 6`` search is performed).
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
FS_LEGEND = 13
LEGEND_LINEWIDTH = 2.25
LEGEND_HANDLELEN = 2.75
PLOT_LW = 2.25
ORDER_COLORS = {1: "#2563eb", 2: "#16a34a", 4: "#9333ea", 6: "#ea580c"}
TROTTER_EXTRA_ORDER = 6
Y_LABEL = "Gate depth"


def _analytic_trotter_gate_depth(errors: np.ndarray, p: int) -> np.ndarray:
    """Closed-form vanilla Trotter gate depth (matches ``rt.compute_steps_gate_depth``)."""
    c = float(rt.gate_overhead(p))
    return c ** (1 + 1 / p) * (errors ** (-1.0 / p)) * ((2 / (1 + p)) ** (1.0 / p))


def _style_ax(ax, title=None, ylabel=Y_LABEL):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Precision  $\varepsilon$", fontsize=FS_AXIS)
    ax.set_ylabel(ylabel, fontsize=FS_AXIS)
    if title is not None:
        ax.set_title(title, fontsize=FS_AXIS)
    ax.tick_params(axis="both", labelsize=FS_AXIS - 2)
    ax.grid(True, which="major", ls="-", alpha=0.25)
    ax.grid(True, which="minor", ls=":", alpha=0.12)


def _legend(ax, *, ncol: int = 1, loc: str = "upper right"):
    leg = ax.legend(
        fontsize=FS_LEGEND - 1,
        loc=loc,
        framealpha=0.92,
        ncol=ncol,
        handlelength=LEGEND_HANDLELEN * 0.92,
        handletextpad=0.5,
        borderpad=0.35,
        labelspacing=0.35,
        columnspacing=0.9,
    )
    for line in leg.get_lines():
        line.set_linewidth(LEGEND_LINEWIDTH)
        mk = line.get_marker()
        if mk is not None and str(mk).lower() not in ("none", ""):
            line.set_markersize(8.5)
            line.set_markeredgewidth(0.7)
    return leg


def parse_args():
    parser = argparse.ArgumentParser(description="Gate-depth figure (per-order + envelope).")
    parser.add_argument("--out-dir", type=str, default="plots", help="Output directory.")
    parser.add_argument("--output", type=str, default="gate_depth.png", help="Titled output filename.")
    add_shared_grid_args(parser)
    parser.set_defaults(eps_log_max=np.log10(0.9))
    parser.add_argument("--orders", type=str, default="1,2,4", help="Comma-separated Trotter orders p.")
    add_shared_search_args(parser)
    parser.add_argument(
        "--n-sys",
        type=int,
        default=100,
        help="System size n passed to compute_steps_gate_depth (default 100).",
    )
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
    omit_titles: bool,
) -> plt.Figure:
    fig_h = 5.2
    fig, axes = plt.subplots(1, 2, figsize=(14.5, fig_h), constrained_layout=True)
    ax_left, ax_right = axes

    for p in orders:
        col = ORDER_COLORS.get(p, "#334155")
        r_p = results_by_order[p]
        trot = r_p.get("wc", r_p["opt"])["trotter"]
        ro = r_p["opt"]["richardson"]
        rw = r_p["wc"]["richardson"] if "wc" in r_p else ro
        rich_best = np.minimum(rw, ro) if "wc" in r_p else ro
        ax_left.plot(
            errors, trot, "-", color=col, linewidth=PLOT_LW,
            label=rf"Trotter $p={p}$",
        )
        ax_left.plot(
            errors, rich_best, "--", color=col, linewidth=PLOT_LW,
            label=rf"Extrapolated $p={p}$",
        )

    p6 = TROTTER_EXTRA_ORDER
    col6 = ORDER_COLORS[p6]
    trot6 = _analytic_trotter_gate_depth(errors, p6)
    ax_left.plot(
        errors, trot6, "-", color=col6, linewidth=PLOT_LW,
        label=rf"Trotter $p={p6}$",
    )

    _style_ax(ax_left, ylabel=Y_LABEL)
    _legend(ax_left, ncol=2, loc="upper right")

    trotter_envs = [results_by_order[p].get("wc", results_by_order[p]["opt"])["trotter"] for p in orders] + [trot6]
    rich_mins = []
    for p in orders:
        r_p = results_by_order[p]
        ro = r_p["opt"]["richardson"]
        if "wc" in r_p:
            rich_mins.append(np.minimum(r_p["wc"]["richardson"], ro))
        else:
            rich_mins.append(ro)
    trotter_best = np.minimum.reduce(trotter_envs)
    rich_best_env = np.minimum.reduce(rich_mins)
    best_trot_p_tex = ",".join(str(x) for x in sorted(set(orders) | {TROTTER_EXTRA_ORDER}))

    for p in orders:
        col = ORDER_COLORS.get(p, "#334155")
        r_p = results_by_order[p]
        ax_right.plot(errors, r_p.get("wc", r_p["opt"])["trotter"], "-",
                      color=col, linewidth=1.4, alpha=0.2, zorder=1)
        ro = r_p["opt"]["richardson"]
        rw = r_p["wc"]["richardson"] if "wc" in r_p else ro
        rich_line = np.minimum(rw, ro) if "wc" in r_p else ro
        ax_right.plot(
            errors,
            rich_line,
            "--", color=col, linewidth=1.4, alpha=0.2, zorder=1,
        )
    ax_right.plot(errors, trot6, "-", color=col6, linewidth=1.4, alpha=0.2, zorder=1)

    ext_p_tex = ",".join(str(p) for p in orders)
    ax_right.plot(
        errors, trotter_best, "k-", linewidth=PLOT_LW + 0.35, zorder=4,
        label=rf"Best Trotter $p \in \{{{best_trot_p_tex}\}}$",
    )
    ax_right.plot(
        errors, rich_best_env, "k--",
        linewidth=PLOT_LW + 0.35,
        marker="s", markersize=5.5,
        markevery=max(1, len(errors) // 10),
        markeredgecolor="black", markeredgewidth=0.6, markerfacecolor="white",
        zorder=4, label=rf"Best extrapolated $p \in \{{{ext_p_tex}\}}$",
    )
    _style_ax(ax_right, ylabel=Y_LABEL)
    leg_r = _legend(ax_right, ncol=1, loc="upper right")
    for line in leg_r.get_lines():
        line.set_linewidth(PLOT_LW + 0.35)

    return fig


SUMMARY_FS = 18
SUMMARY_CBAR_VMIN = 1.0
SUMMARY_CBAR_VMAX = 100.0
SUMMARY_CBAR_LABEL = r"$\|\mathbf{b}\|_1^2$  (sample overhead factor)"


def _envelopes(orders, errors, results_by_order):
    """Best-Trotter and best-extrapolated envelopes, plus the sample overhead
    ``‖b‖₁²`` of the winning Richardson schedule at each ε."""
    trot6 = _analytic_trotter_gate_depth(errors, TROTTER_EXTRA_ORDER)
    trot_curves = [
        results_by_order[p].get("wc", results_by_order[p]["opt"])["trotter"]
        for p in orders
    ] + [trot6]
    trot_env = np.minimum.reduce(trot_curves)

    rich_env = np.full(len(errors), np.inf)
    so_env = np.full(len(errors), np.nan)
    for p in orders:
        r_p = results_by_order[p]
        for mode in ("wc", "opt"):
            if mode not in r_p:
                continue
            steps = np.asarray(r_p[mode]["richardson"], dtype=float)
            so = np.asarray(r_p[mode]["sample_overhead"], dtype=float)
            better = steps < rich_env
            rich_env = np.where(better, steps, rich_env)
            so_env = np.where(better, so, so_env)
    return trot_env, rich_env, so_env


def _build_summary_figure(
    *,
    orders: list[int],
    errors: np.ndarray,
    results_by_order: dict,
) -> plt.Figure:
    """Standalone square version of the right (envelope) panel: gray Trotter
    line + best-extrapolated triangles colored by sample overhead."""
    trot_env, rich_env, so_env = _envelopes(orders, errors, results_by_order)
    trot_p_tex = ",".join(str(x) for x in sorted(set(orders) | {TROTTER_EXTRA_ORDER}))
    ext_p_tex = ",".join(str(p) for p in orders)

    norm = mcolors.LogNorm(vmin=SUMMARY_CBAR_VMIN, vmax=SUMMARY_CBAR_VMAX, clip=False)
    cmap = plt.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(6.6, 6.6), constrained_layout=True)
    ax.plot(
        errors, trot_env, "-", color="gray", linewidth=PLOT_LW + 0.35, zorder=2,
    )
    sc = ax.scatter(
        errors, rich_env, c=so_env, cmap=cmap, norm=norm,
        marker="^", s=70, linewidths=0, edgecolors="none", zorder=3,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Precision  $\varepsilon$", fontsize=SUMMARY_FS)
    ax.set_ylabel(Y_LABEL, fontsize=SUMMARY_FS)
    ax.tick_params(axis="both", labelsize=SUMMARY_FS)
    ax.grid(True, which="major", ls="-", alpha=0.25)
    ax.grid(True, which="minor", ls=":", alpha=0.12)
    ax.set_box_aspect(1)

    handles = [
        Line2D(
            [0], [0], color="gray", linewidth=PLOT_LW + 0.35,
            label=rf"Best Trotter $p \in \{{{trot_p_tex}\}}$",
        ),
        Line2D(
            [0], [0], marker="^", linestyle="None", markersize=10,
            markerfacecolor=cmap(0.6), markeredgecolor="none",
            label=rf"Best extrapolated $p \in \{{{ext_p_tex}\}}$",
        ),
    ]
    ax.legend(handles=handles, fontsize=SUMMARY_FS, loc="upper right", framealpha=0.92)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(SUMMARY_CBAR_LABEL, fontsize=SUMMARY_FS)
    cbar.ax.tick_params(labelsize=SUMMARY_FS)
    return fig


def main() -> None:
    args = parse_args()
    rt.set_lambda_scale_mode(args.lambda_mode)
    plt.rcParams["pdf.fonttype"] = 42

    output_dir = resolve_output_dir(_ROOT, args.out_dir)
    orders = parse_orders(args.orders)
    errors = np.logspace(args.eps_log_max, args.eps_log_min, num=int(args.eps_points))
    b2_cap = float(args.brute_bnorm_sq_max)
    n_sys = int(args.n_sys)

    print(
        f"Gate depth: eps in [{10**args.eps_log_min:g}, {10**args.eps_log_max:g}] "
        f"({args.eps_points} pts), q in [{args.q_min},{args.q_max}], "
        f"||b||_1^2 cap={b2_cap:g}, n_sys={n_sys}, orders={orders}, "
        f"lambda-mode={args.lambda_mode}"
    )

    results_by_order = cm.compute_results(
        errors=errors,
        orders=orders,
        q_max=int(args.q_max),
        q_min=int(args.q_min),
        b2_cap=b2_cap,
        brute_permutations=bool(args.brute_permutations),
        brute_permutations_max_count=int(args.brute_permutations_max_count),
        step_mode="gate_depth",
        n_sys=n_sys,
    )

    def save(omit_titles: bool, name: str) -> None:
        fig = _build_figure(
            orders=orders, errors=errors, results_by_order=results_by_order,
            omit_titles=omit_titles,
        )
        for path in cm.save_figure(fig, output_dir / name):
            print(f"  Saved {path}")
        plt.close(fig)

    stem = Path(args.output)
    save(omit_titles=False, name=args.output)
    if not args.no_cropped:
        save(omit_titles=True, name=f"{stem.stem}_cropped{stem.suffix}")

    summary_fig = _build_summary_figure(
        orders=orders, errors=errors, results_by_order=results_by_order,
    )
    for path in cm.save_figure(summary_fig, output_dir / f"{stem.stem}_summary{stem.suffix}"):
        print(f"  Saved {path}")
    plt.close(summary_fig)

    p6 = TROTTER_EXTRA_ORDER
    c6 = int(rt.gate_overhead(p6))
    md_path, json_path = cm.write_params_sidecar(
        output_dir=output_dir,
        stem=stem.stem,
        title="Gate-depth figure parameters",
        settings={
            "script": "plotting/plot_gate_depth.py",
            "step_mode": "gate_depth",
            "λ-mode": args.lambda_mode,
            "orders (p) searched": ", ".join(str(p) for p in orders),
            "ε grid": (
                f"{args.eps_points} points, log10 in "
                f"[{args.eps_log_min:.6g}, {args.eps_log_max:.6g}]"
            ),
            "q_min, q_max": f"{args.q_min}, {args.q_max}",
            "‖b‖₁² search cap": f"{b2_cap:g}",
            "brute permutations": str(bool(args.brute_permutations)),
            "system size n (λ-comm)": str(n_sys),
        },
        errors=errors,
        orders=orders,
        results_by_order=results_by_order,
        extras={
            f"Trotter p = {p6} (analytic, no Richardson)": (
                f"closed form C_p^(1+1/p) · ε^(-1/p) · (2/(1+p))^(1/p) with C_{p6} = {c6}"
            ),
            "envelope (right panel)": (
                "Best Trotter taken over p ∈ "
                f"{{{', '.join(str(x) for x in sorted(set(orders) | {p6}))}}}; "
                "Best extrapolated taken over p ∈ "
                f"{{{', '.join(str(x) for x in orders)}}} and both schedules."
            ),
        },
    )
    print(f"  Saved {md_path}")
    print(f"  Saved {json_path}")


if __name__ == "__main__":
    main()
