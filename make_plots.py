"""Regenerate the paper's figures (and their parameter sidecars) under ``plots/``.

Usage::

    python make_plots.py                      # Figs. 2, 4, 5 of arXiv:2608.13862
    python make_plots.py -- --eps-points 30   # forward extra flags to the two search scripts

Runs, in order:

    plotting/plot_overhead.py    -> plots/overhead_multi_cap.{png,pdf}   (Fig. 4)
                                    + overhead.params.* and overhead_multi_cap.params.*
    plotting/plot_gate_depth.py  -> plots/gate_depth.{png,pdf}           (Fig. 5)
                                    + gate_depth.params.* and gate_depth_multi_cap.params.*
    plotting/plot_summary.py     -> plots/summary.{png,pdf}              (Fig. 2)
                                    assembled from the sidecars above

Extra arguments after ``--`` are forwarded to the two search scripts only; the
summary takes none. The exact-error figure (Fig. 6) lives in ``error_analysis/``
and is added here once its script exists.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent

# (script, fixed arguments): the settings the paper's figures were made with.
_STEPS: list[tuple[Path, list[str]]] = [
    (_ROOT / "plotting" / "plot_overhead.py",
     ["--brute-bnorm-sq-caps", "10,100,1000", "--no-cropped"]),
    (_ROOT / "plotting" / "plot_gate_depth.py", ["--no-cropped"]),
    # reads the sidecars written above; must run last
    (_ROOT / "plotting" / "plot_summary.py", []),
]


def _extra_args(argv: list[str]) -> list[str]:
    if "--" in argv:
        idx = argv.index("--")
        return argv[idx + 1 :]
    return []


def main() -> int:
    extra = _extra_args(sys.argv[1:])
    for script, fixed in _STEPS:
        print(f"\n=== Running {script.relative_to(_ROOT)} ===")
        cmd = [sys.executable, str(script), *fixed]
        if script.name != "plot_summary.py":
            cmd.extend(extra)
        rc = subprocess.call(cmd, cwd=_ROOT)
        if rc != 0:
            return rc
    print("\nAll figures regenerated in plots/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
