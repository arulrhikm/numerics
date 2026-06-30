"""Side-by-side square summary figure assembled from the two committed sidecars.

Left panel  : the gate-depth envelope (right panel of the gate-depth figure),
              drawn in black only -- "Best Trotter" and "Best extrapolated"
              envelopes, with the faint per-order colored traces removed.
Right panel : the p = 2 panel of the overhead figure -- the Trotter baseline
              plus the LKW (squares) and brute-force (triangles) Richardson
              schedules colored by sample overhead ``‖b‖₁²`` (viridis colorbar,
              same 1..1000 LogNorm as the overhead figure).

Both panels are rendered square so the figure sits beside the exact-error plot.

The data is read straight from ``plots/overhead.params.json`` and
``plots/gate_depth.params.json`` (so the composite matches the published panels
exactly); regenerate those first via ``make_plots.py``.
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

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import richardson as rt
from plotting import common as cm
from plotting.common_cli import resolve_output_dir

FS = 17
PLOT_LW = 2.25
TROTTER_EXTRA_ORDER = 6
OVERHEAD_PANEL_P = 2
CBAR_VMIN = 1.0
CBAR_VMAX = 1000.0
CBAR_LABEL = r"$\|\mathbf{b}\|_1^2$  (sample overhead factor)"


def _style(ax, *, xlabel=True, ylabel=None, title=None):
    ax.set_xscale("log")
    ax.set_yscale("log")
    if xlabel:
        ax.set_xlabel(r"Precision  $\varepsilon$", fontsize=FS)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=FS)
    if title is not None:
        ax.set_title(title, fontsize=FS, pad=8)
    ax.tick_params(axis="both", labelsize=FS)
    ax.grid(True, which="major", ls="-", alpha=0.25)
    ax.grid(True, which="minor", ls=":", alpha=0.12)
    ax.set_box_aspect(1)


def _gate_depth_envelopes(g: dict):
    eps = np.array(g["epsilon"])
    orders = g["orders"]
    trot6 = rt.gate_overhead(TROTTER_EXTRA_ORDER) ** (1 + 1 / TROTTER_EXTRA_ORDER) * (
        eps ** (-1.0 / TROTTER_EXTRA_ORDER)
    ) * ((2 / (1 + TROTTER_EXTRA_ORDER)) ** (1.0 / TROTTER_EXTRA_ORDER))
    trot_curves = [trot6]
    rich_env = np.full(len(eps), np.inf)
    for p in orders:
        rp = g["results"][str(p)]
        trot_curves.append(np.array(rp.get("wc", rp["opt"])["trotter_steps"]))
        ro = np.array(rp["opt"]["richardson_steps"])
        if "wc" in rp:
            ro = np.minimum(ro, np.array(rp["wc"]["richardson_steps"]))
        rich_env = np.minimum(rich_env, ro)
    return eps, orders, np.minimum.reduce(trot_curves), rich_env


def _draw_gate_depth_panel(ax, g: dict):
    eps, orders, trot_env, rich_env = _gate_depth_envelopes(g)
    trot_p_tex = ",".join(str(x) for x in sorted(set(orders) | {TROTTER_EXTRA_ORDER}))
    ext_p_tex = ",".join(str(p) for p in orders)
    ax.plot(
        eps, trot_env, "k-", linewidth=PLOT_LW + 0.35, zorder=4,
        label=rf"Best Trotter $p \in \{{{trot_p_tex}\}}$",
    )
    ax.plot(
        eps, rich_env, "k--", linewidth=PLOT_LW + 0.35, marker="s", markersize=5.5,
        markevery=max(1, len(eps) // 10), markeredgecolor="black", markeredgewidth=0.6,
        markerfacecolor="white", zorder=4,
        label=rf"Best extrapolated $p \in \{{{ext_p_tex}\}}$",
    )
    _style(ax, ylabel="Gate depth")
    ax.legend(fontsize=FS, loc="upper right", framealpha=0.92)


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
    ax.legend(handles=handles, fontsize=FS - 3, loc="upper left", framealpha=0.9)
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

    g = json.loads((output_dir / "gate_depth.params.json").read_text(encoding="utf-8"))
    o = json.loads((output_dir / "overhead.params.json").read_text(encoding="utf-8"))

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13, 6.2), constrained_layout=True)
    _draw_gate_depth_panel(ax_left, g)
    mappable = _draw_overhead_panel(ax_right, o, int(args.overhead_panel_p))
    cbar = fig.colorbar(mappable, ax=ax_right, fraction=0.046, pad=0.02)
    cbar.set_label(CBAR_LABEL, fontsize=FS - 1)
    cbar.ax.tick_params(labelsize=FS - 2)

    for path in cm.save_figure(fig, output_dir / args.output):
        print(f"  Saved {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
