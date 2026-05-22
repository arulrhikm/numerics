"""Regenerate every figure (and parameter sidecar) under ``plots/`` in one command.

Usage::

    python make_plots.py            # default settings
    python make_plots.py -- --eps-points 30   # forward extra flags to BOTH scripts

Extra arguments after ``--`` are forwarded to each underlying plotting script.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
_SCRIPTS = [
    _ROOT / "plotting" / "plot_overhead.py",
    _ROOT / "plotting" / "plot_gate_depth.py",
]


def _extra_args(argv: list[str]) -> list[str]:
    if "--" in argv:
        idx = argv.index("--")
        return argv[idx + 1 :]
    return []


def main() -> int:
    extra = _extra_args(sys.argv[1:])
    for script in _SCRIPTS:
        print(f"\n=== Running {script.relative_to(_ROOT)} ===")
        cmd = [sys.executable, str(script), *extra]
        rc = subprocess.call(cmd, cwd=_ROOT)
        if rc != 0:
            return rc
    print("\nAll plots regenerated in plots/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
