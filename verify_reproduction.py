"""Check that a regenerated run reproduces the committed parameter sidecars.

Usage::

    python make_plots.py -- --out-dir _verify     # ~tens of minutes, a few GB
    python verify_reproduction.py _verify

Compares, for each sidecar in ``plots/``, the ε grid and every searched
quantity (schedules, coefficients, norms, step counts) against the same file in
the given directory, and prints the first differing entry per block. The search
is deterministic, so a run at the same settings should report no differences.

``--settings-only`` skips the numerical comparison and just prints the settings
block of each pair, which is the quick way to see whether two runs used the
same ``q_max``, cap list and ε grid.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Settings keys carry λ and ‖·‖; the Windows console defaults to cp1252.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

_ROOT = Path(__file__).resolve().parent
_SIDECARS = [
    "overhead.params.json",
    "overhead_multi_cap.params.json",
    "gate_depth.params.json",
    "gate_depth_multi_cap.params.json",
]
_FIELDS = ("m", "q_grids", "b_coeffs", "bnorm1", "bnorm1_sq", "btilde1",
           "trotter_steps", "richardson_steps")


def _blocks(results: dict) -> dict[tuple[str, str, str], dict]:
    """Flatten results[p][mode] (and the per-cap blocks) to one dict per key."""
    out: dict[tuple[str, str, str], dict] = {}
    for p, by_mode in results.items():
        for mode, block in by_mode.items():
            if mode == "opt_by_cap":
                for cap, sub in block.items():
                    out[(p, mode, cap)] = sub
            else:
                out[(p, mode, "")] = block
    return out


def _compare(committed: Path, fresh: Path) -> int:
    a = json.loads(committed.read_text(encoding="utf-8"))
    b = json.loads(fresh.read_text(encoding="utf-8"))
    diffs = 0

    if a.get("settings") != b.get("settings"):
        print(f"  settings differ:")
        for k in sorted(set(a.get("settings", {})) | set(b.get("settings", {}))):
            va, vb = a.get("settings", {}).get(k), b.get("settings", {}).get(k)
            if va != vb:
                print(f"    {k}: committed={va!r} fresh={vb!r}")
        diffs += 1

    if a.get("epsilon") != b.get("epsilon"):
        print(f"  ε grid differs: {len(a.get('epsilon', []))} vs {len(b.get('epsilon', []))} points")
        return diffs + 1

    ba, bb = _blocks(a["results"]), _blocks(b["results"])
    if set(ba) != set(bb):
        print(f"  blocks differ: committed-only {sorted(set(ba) - set(bb))}, "
              f"fresh-only {sorted(set(bb) - set(ba))}")
        diffs += 1

    for key in sorted(set(ba) & set(bb)):
        for field in _FIELDS:
            va, vb = ba[key].get(field), bb[key].get(field)
            if va == vb:
                continue
            label = "/".join(x for x in key if x)
            if isinstance(va, list) and isinstance(vb, list) and len(va) == len(vb):
                i = next(i for i, (x, y) in enumerate(zip(va, vb)) if x != y)
                print(f"  {label} {field}: first difference at index {i}: "
                      f"committed={va[i]!r} fresh={vb[i]!r}")
            else:
                print(f"  {label} {field}: differs")
            diffs += 1
            break  # one report per block is enough to localize it
    return diffs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fresh_dir", help="directory holding a regenerated run (e.g. _verify)")
    ap.add_argument("--committed-dir", default="plots")
    ap.add_argument("--settings-only", action="store_true")
    args = ap.parse_args()

    fresh_dir = (_ROOT / args.fresh_dir).resolve()
    committed_dir = (_ROOT / args.committed_dir).resolve()
    total = 0
    for name in _SIDECARS:
        committed, fresh = committed_dir / name, fresh_dir / name
        print(f"\n{name}")
        if not fresh.exists():
            print(f"  missing in {fresh_dir} (was the run cut short?)")
            total += 1
            continue
        if args.settings_only:
            print(f"  committed: {json.loads(committed.read_text(encoding='utf-8'))['settings']}")
            print(f"  fresh    : {json.loads(fresh.read_text(encoding='utf-8'))['settings']}")
            continue
        n = _compare(committed, fresh)
        print("  identical" if n == 0 else f"  {n} differing block(s)")
        total += n

    if args.settings_only:
        return 0
    print(f"\n{'Reproduced: no differences.' if total == 0 else f'{total} difference(s) found.'}")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
