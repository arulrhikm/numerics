"""Side-by-side square summary figure assembled from the two committed sidecars.

Left panel  : the gate-depth envelope -- gray "Best Trotter" plus two
              "Best extrapolated" envelopes (‖b‖₁² caps 10 and 100), each
              point a filled triangle colored by sample overhead ``‖b‖₁²``.
              Cap budgets are distinguished by the shared colorbar, not by
              open vs filled markers.
Right panel : the p = 2 overhead panel -- Trotter baseline plus LKW (squares)
              and brute-force (triangles) schedules, also colored by ``‖b‖₁²``.

Both panels share the sample-overhead colorbar on the right.

The data is read from ``plots/overhead.params.json`` and
``plots/gate_depth_multi_cap.params.json``; regenerate via ``make_plots.py``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, LogFormatterMathtext

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import richardson as rt
from plotting import common as cm
from plotting.common_cli import resolve_output_dir

FS = 20                    # poster scale: summary.pdf is placed at ~1.1x on an A0 board
PLOT_LW = 2.25
TROTTER_EXTRA_ORDER = 6
OVERHEAD_PANEL_P = 2
CBAR_VMIN = 1.0
CBAR_VMAX = 1000.0
CBAR_LABEL = r"$\|\mathbf{b}\|_1^2$  (sample overhead factor)"
# Sample-overhead budgets drawn as separate "Best extrapolated" envelopes in the
# left panel: the published 100-sample line plus the tighter 10-sample line.
# Both crossovers with the Trotter baseline land above ε = 1e-2.
SUMMARY_CAPS = [100.0, 10.0]
# ε decades shown on both panels (matches the sidecar grid, 10^-6 … 10^0).
# Every OTHER decade. At poster font size (FS=20) all seven decade labels
# collide into an unreadable smear; four are enough to read the axis.
SUMMARY_EPS_XTICKS = [10.0 ** e for e in range(-6, 1, 2)]


def _cap_key(cap: float) -> str:
    """Match the ``opt_by_cap`` JSON key written by ``write_params_sidecar``."""
    return f"{cap:g}"


def _style(ax, *, xlabel=True, ylabel=None, title=None):
    ax.set_xscale("log")
    ax.set_yscale("log")
    if xlabel:
        ax.set_xlabel(r"Precision  $\varepsilon$", fontsize=FS)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=FS)
    if title is not None:
        ax.set_title(title, fontsize=FS, pad=8)
    ax.tick_params(axis="both", labelsize=FS - 2)
    ax.set_xticks(SUMMARY_EPS_XTICKS)
    ax.xaxis.set_major_formatter(LogFormatterMathtext())
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10))))
    ax.grid(True, which="major", ls="-", alpha=0.25)
    ax.grid(True, which="minor", ls=":", alpha=0.12)
    ax.set_box_aspect(1)


def _gate_depth_envelopes(g: dict, caps):
    """Best-Trotter envelope plus, for each sample-overhead cap, the
    best-extrapolated envelope and the winning schedule's ``‖b‖₁²``.

    The extrapolated envelope is the pointwise minimum over ``p`` of the
    brute-force optimized ``opt_by_cap`` schedule (the well-conditioned grid
    never wins on gate depth, so it is not part of the envelope).
    """
    eps = np.array(g["epsilon"])
    orders = g["orders"]
    trot6 = rt.gate_overhead(TROTTER_EXTRA_ORDER) ** (1 + 1 / TROTTER_EXTRA_ORDER) * (
        eps ** (-1.0 / TROTTER_EXTRA_ORDER)
    ) * ((2 / (1 + TROTTER_EXTRA_ORDER)) ** (1.0 / TROTTER_EXTRA_ORDER))
    trot_curves = [trot6]
    for p in orders:
        rp = g["results"][str(p)]
        any_block = rp["wc"] if "wc" in rp else next(iter(rp["opt_by_cap"].values()))
        trot_curves.append(np.array(any_block["trotter_steps"]))
    trot_env = np.minimum.reduce(trot_curves)

    env_by_cap: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for cap in caps:
        key = _cap_key(cap)
        rich_env = np.full(len(eps), np.inf)
        so_env = np.full(len(eps), np.nan)
        for p in orders:
            block = g["results"][str(p)]["opt_by_cap"][key]
            steps = np.array(block["richardson_steps"])
            so = np.array(block["bnorm1_sq"])
            better = steps < rich_env
            rich_env = np.where(better, steps, rich_env)
            so_env = np.where(better, so, so_env)
        env_by_cap[cap] = (rich_env, so_env)
    return eps, orders, trot_env, env_by_cap


def _draw_gate_depth_panel(ax, g: dict):
    eps, orders, trot_env, env_by_cap = _gate_depth_envelopes(g, SUMMARY_CAPS)
    trot_p_tex = ",".join(str(x) for x in sorted(set(orders) | {TROTTER_EXTRA_ORDER}))
    ext_p_tex = ",".join(str(p) for p in orders)
    norm = mcolors.LogNorm(vmin=CBAR_VMIN, vmax=CBAR_VMAX, clip=False)
    cmap = plt.get_cmap("viridis")

    ax.plot(
        eps, trot_env, "-", color="gray", linewidth=PLOT_LW + 0.35, zorder=2,
        label=rf"Best Trotter $p \in \{{{trot_p_tex}\}}$",
    )

    rich_pts: list[np.ndarray] = []
    so_pts: list[np.ndarray] = []
    for cap in SUMMARY_CAPS:
        rich, so = env_by_cap[cap]
        rich_pts.append(rich)
        so_pts.append(so)
    rich_all = np.concatenate(rich_pts)
    so_all = np.concatenate(so_pts)
    eps_all = np.concatenate([eps] * len(SUMMARY_CAPS))
    mappable = ax.scatter(
        eps_all, rich_all, c=so_all, cmap=cmap, norm=norm, marker="^", s=70,
        linewidths=0, edgecolors="none", zorder=3,
    )

    _style(ax, ylabel="Gate depth")
    handles = [
        Line2D([0], [0], color="gray", linewidth=PLOT_LW + 0.35,
               label=rf"Best Trotter $p \in \{{{trot_p_tex}\}}$"),
        Line2D([0], [0], marker="^", linestyle="None", markersize=10,
               markerfacecolor="#888", markeredgecolor="none",
               label=rf"Best extrapolated $p \in \{{{ext_p_tex}\}}$"),
    ]
    # FS-6, not FS: at poster font size the legend grows wider than the panel
    # and its opaque box slides over the rotated y-axis label.
    ax.legend(handles=handles, fontsize=FS - 6, loc="upper right", framealpha=0.92)
    return mappable


def _draw_overhead_panel(ax, o: dict, p: int):
    eps = np.array(o["epsilon"])
    norm = mcolors.LogNorm(vmin=CBAR_VMIN, vmax=CBAR_VMAX, clip=False)
    cmap = plt.get_cmap("viridis")
    r_p = o["results"][str(p)]
    ax.plot(
        eps, r_p.get("wc", r_p["opt"])["trotter_steps"], "-",
        color="gray", linewidth=2, alpha=0.35, zorder=1,
    )
    mappable = None
    for mode, marker, z in (("wc", "s", 3), ("opt", "^", 4)):
        if mode not in r_p:
            continue
        sc = ax.scatter(
            eps, np.array(r_p[mode]["richardson_steps"]),
            c=np.array(r_p[mode]["bnorm1_sq"]), cmap=cmap, norm=norm,
            s=70, marker=marker, linewidths=0, edgecolors="none", zorder=z,
        )
        if mode == "opt":
            mappable = sc
    _style(ax, ylabel="Maximum # Trotter steps", title=rf"$p = {p}$")
    ax.set_ylim(1.0, 1e3)
    handles = [
        Line2D([0], [0], color="gray", linewidth=PLOT_LW, alpha=0.45, label="Trotter"),
        Line2D([0], [0], marker="s", linestyle="None", markersize=8,
               markerfacecolor="#888", markeredgecolor="none", label="LKW well conditioned"),
        Line2D([0], [0], marker="^", linestyle="None", markersize=8,
               markerfacecolor="#444", markeredgecolor="none", label="Brute force optimization"),
    ]
    ax.legend(handles=handles, fontsize=FS - 6, loc="upper left", framealpha=0.9)
    return mappable


def parse_args():
    parser = argparse.ArgumentParser(description="Side-by-side numerics summary figure.")
    parser.add_argument("--out-dir", type=str, default="plots", help="Output/sidecar directory.")
    parser.add_argument("--output", type=str, default="summary.png", help="Output filename.")
    parser.add_argument(
        "--overhead-panel-p", type=int, default=OVERHEAD_PANEL_P,
        help="Which overhead order to show in the right panel.",
    )
    # tolerate (and ignore) the shared search flags forwarded by make_plots.py
    args, _unknown = parser.parse_known_args()
    return args


def main() -> None:
    args = parse_args()
    plt.rcParams["pdf.fonttype"] = 42
    output_dir = resolve_output_dir(_ROOT, args.out_dir)

    g = json.loads((output_dir / "gate_depth_multi_cap.params.json").read_text(encoding="utf-8"))
    o = json.loads((output_dir / "overhead.params.json").read_text(encoding="utf-8"))

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    left_mappable = _draw_gate_depth_panel(ax_left, g)
    right_mappable = _draw_overhead_panel(ax_right, o, int(args.overhead_panel_p))
    mappable = right_mappable if right_mappable is not None else left_mappable
    cbar = fig.colorbar(
        mappable, ax=[ax_left, ax_right], fraction=0.046, pad=0.02,
    )
    cbar.set_label(CBAR_LABEL, fontsize=FS - 1)
    cbar.ax.tick_params(labelsize=FS - 2)

    for path in cm.save_figure(fig, output_dir / args.output):
        print(f"  Saved {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
