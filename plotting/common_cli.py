"""Shared CLI helpers for the plotting scripts."""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_Q_MAX = 15


def add_shared_grid_args(parser: argparse.ArgumentParser, default_q_max: int | None = None) -> None:
    """Add the ε / q grid arguments shared by both plot scripts."""
    if default_q_max is None:
        default_q_max = DEFAULT_Q_MAX
    parser.add_argument("--q-max", type=int, default=default_q_max)
    parser.add_argument("--q-min", type=int, default=1)
    parser.add_argument("--brute-bnorm-sq-max", type=float, default=100.0)
    parser.add_argument("--eps-log-min", type=float, default=-6.0)
    parser.add_argument("--eps-log-max", type=float, default=0.0)
    parser.add_argument("--eps-points", type=int, default=50)


def add_shared_search_args(parser: argparse.ArgumentParser) -> None:
    """Add the brute-force / λ-mode arguments shared by both plot scripts."""
    parser.add_argument("--brute-permutations", action="store_true")
    parser.add_argument("--brute-permutations-max-count", type=int, default=2_000_000)
    parser.add_argument(
        "--lambda-mode",
        default="lemma57_fixed",
        choices=["lemma57_fixed", "legacy"],
        help="λ-comm model to use.",
    )


def resolve_output_dir(root_dir: Path, out_dir_arg: str) -> Path:
    """Return ``out_dir_arg`` resolved against ``root_dir`` (and mkdir -p)."""
    out = Path(out_dir_arg)
    if not out.is_absolute():
        out = root_dir / out
    out.mkdir(parents=True, exist_ok=True)
    return out


def parse_orders(orders_arg: str) -> list[int]:
    """Parse ``"1,2,4"`` -> ``[1, 2, 4]`` (rejects empty)."""
    orders = [int(x.strip()) for x in orders_arg.split(",") if x.strip()]
    if not orders:
        raise SystemExit("--orders must list at least one integer p.")
    return orders


def parse_brute_bnorm_sq_caps(caps_arg: str) -> list[float]:
    """Parse ``"10,100,1000"`` -> ``[10.0, 100.0, 1000.0]`` (rejects empty)."""
    caps = [float(x.strip()) for x in caps_arg.split(",") if x.strip()]
    if not caps:
        raise SystemExit("--brute-bnorm-sq-caps must list at least one cap.")
    return caps


def format_b2_cap(cap: float) -> str:
    """Human-readable cap for titles, legends, and sidecars."""
    if abs(cap - round(cap)) < 1e-9:
        return str(int(round(cap)))
    return f"{cap:g}"
