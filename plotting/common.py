"""Shared search-and-step computation for the overhead and gate-depth figures.

Runs ``compute_min_samples`` (well-conditioned grid + brute-force optimized grid)
for each Trotter order ``p``, converts the schedules into step counts under the
requested ``step_mode``, and exposes a parameter-sidecar writer that dumps a
Markdown + JSON summary of every grid the search picked.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import numpy as np

import richardson as rt

MODE_LABEL = {
    "wc": "well-conditioned grid",
    "opt": "brute-force optimized grid",
}

StepMode = Literal["plane_wave", "gate_depth"]


def visible_mask_for_sample_overhead(
    *,
    order: int,
    sample_overhead: np.ndarray,
    omit_p1_sample_overhead_above: float,
) -> np.ndarray:
    """Hide ill-conditioned ``p = 1`` scatter points above the threshold."""
    mask = np.isfinite(sample_overhead)
    threshold = float(omit_p1_sample_overhead_above)
    if order == 1 and np.isfinite(threshold) and threshold > 0.0:
        mask &= sample_overhead <= threshold
    return mask


def _steps_for_error(
    *, step_mode: StepMode, error: float, m: int, m_expr: float, p: int,
    A: float, n_sys: int,
) -> tuple[float, float]:
    if step_mode == "plane_wave":
        return rt.compute_steps_plane_wave(error, m, m_expr, p, A=A)
    return rt.compute_steps_gate_depth(error, m, m_expr, p, A=A, n=n_sys)


def compute_results(
    *,
    errors: np.ndarray,
    orders: list[int],
    q_max: int,
    q_min: int,
    b2_cap: float,
    brute_permutations: bool,
    brute_permutations_max_count: int,
    step_mode: StepMode,
    n_sys: int,
    A: float = 1.0,
) -> dict:
    """Return ``{p: {"wc"|"opt": {...}}}`` with step counts and grid diagnostics."""
    results_by_order: dict = {}
    schedules = [
        ("wc", dict(well_conditioned_formula=True)),
        ("opt", dict(well_conditioned_formula=False)),
    ]

    for p in orders:
        results_by_order[p] = {}
        for mode_name, mode_flags in schedules:
            if mode_name == "wc" and not rt.wc_extrapolation_valid_for_p(p):
                print(f"  p={p}, mode={mode_name}: skipped (WC requires even p)")
                continue
            print(f"  p={p}, mode={mode_name} ({step_mode}) ...")
            best_m, best_base, best_q = rt.compute_min_samples(
                errors,
                p=p,
                m_max=q_max,
                q_min=q_min,
                q_max=q_max,
                brute_permutations=brute_permutations,
                brute_permutations_max_count=brute_permutations_max_count,
                brute_force_b_norm1_sq_max=b2_cap,
                A=A,
                **mode_flags,
            )

            trotter_steps: list[float] = []
            richardson_steps: list[float] = []
            bnorm1: list[float] = []
            btilde1: list[float] = []
            b_coeffs: list[list[float]] = []
            q_grids: list[list[int]] = []
            use_wc = mode_name == "wc"
            for eps, m, m_expr, q_grid in zip(errors, best_m, best_base, best_q):
                t, r = _steps_for_error(
                    step_mode=step_mode, error=eps, m=m, m_expr=m_expr, p=p,
                    A=A, n_sys=n_sys,
                )
                trotter_steps.append(t)
                richardson_steps.append(r)
                s_list = [1.0 / q for q in q_grid]
                b, _ = (
                    rt.get_wc_richardson_coefficients(s_list, m)
                    if use_wc
                    else rt.get_richardson_coefficients(s_list, p, m)
                )
                bnorm1.append(rt.b_norm1(b))
                btilde1.append(rt.b_suppressed_norm(b, q_grid, m, p))
                b_coeffs.append([float(x) for x in b])
                q_grids.append([int(q) for q in q_grid])

            bnorm1_arr = np.array(bnorm1, dtype=float)
            results_by_order[p][mode_name] = dict(
                trotter=np.array(trotter_steps, dtype=float),
                richardson=np.array(richardson_steps, dtype=float),
                sample_overhead=bnorm1_arr**2,
                m_expr=np.array(best_base, dtype=float),
                m_per_error=np.array(best_m, dtype=int),
                bnorm1=bnorm1_arr,
                btilde1=np.array(btilde1, dtype=float),
                b_coeffs=b_coeffs,
                q_grids=q_grids,
            )
    return results_by_order


def write_params_sidecar(
    *,
    output_dir: Path,
    stem: str,
    title: str,
    settings: dict,
    errors: np.ndarray,
    orders: list[int],
    results_by_order: dict,
    extras: dict | None = None,
) -> tuple[Path, Path]:
    """Write ``<stem>.params.md`` + ``<stem>.params.json`` next to the figure.

    The Markdown is a human summary; the JSON has the full numerical data
    (Richardson coefficients ``b_k`` and per-curve step counts).
    """
    md_path = output_dir / f"{stem}.params.md"
    json_path = output_dir / f"{stem}.params.json"
    extras = extras or {}

    out: list[str] = [
        f"# {title}",
        "",
        f"Companion data for `plots/{stem}.png` (and `plots/{stem}_cropped.png`).",
        "",
        "## Settings",
        "",
        "| Parameter | Value |",
        "| --- | --- |",
    ]
    out.extend(f"| {k} | {v} |" for k, v in settings.items())

    if extras:
        out += ["", "## Additional curves", ""]
        out.extend(f"- **{k}**: {v}" for k, v in extras.items())

    out += [
        "",
        "## Optimal Richardson grids per (p, ε)",
        "",
        "For each Trotter order `p`, the search returned the schedule below "
        "(Richardson depth `m`, integer refinement grid `q_k`, and the "
        "coefficient norms `‖b‖₁`, `‖b‖₁²`, and the suppressed norm `‖b̃‖₁` "
        "that enters the step bound).",
        "",
    ]

    json_data: dict = {
        "title": title,
        "settings": settings,
        "extras": extras,
        "epsilon": [float(e) for e in errors],
        "orders": list(orders),
        "results": {},
    }

    for p in orders:
        json_data["results"][str(p)] = {}
        for mode in ("wc", "opt"):
            if mode not in results_by_order[p]:
                continue
            r = results_by_order[p][mode]
            out += [
                f"### p = {p}, {MODE_LABEL[mode]} (`mode = {mode}`)",
                "",
                "| ε | m | q_k | ‖b‖₁ | ‖b‖₁² | ‖b̃‖₁ |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
            for i, eps in enumerate(errors):
                q_str = "[" + ", ".join(str(q) for q in r["q_grids"][i]) + "]"
                out.append(
                    f"| {eps:.3e} | {int(r['m_per_error'][i])} | {q_str} | "
                    f"{r['bnorm1'][i]:.4g} | {r['sample_overhead'][i]:.4g} | "
                    f"{r['btilde1'][i]:.4g} |"
                )
            out.append("")

            json_data["results"][str(p)][mode] = {
                "m": [int(x) for x in r["m_per_error"]],
                "q_grids": r["q_grids"],
                "b_coeffs": r["b_coeffs"],
                "bnorm1": [float(x) for x in r["bnorm1"]],
                "bnorm1_sq": [float(x) for x in r["sample_overhead"]],
                "btilde1": [float(x) for x in r["btilde1"]],
                "trotter_steps": [float(x) for x in r["trotter"]],
                "richardson_steps": [float(x) for x in r["richardson"]],
            }

    out += [
        f"Full Richardson coefficients `b_k` and per-curve step counts are in "
        f"`{json_path.name}`.",
        "",
    ]
    md_path.write_text("\n".join(out), encoding="utf-8")
    json_path.write_text(json.dumps(json_data, indent=2), encoding="utf-8")
    return md_path, json_path
