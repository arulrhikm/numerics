"""Regenerate every figure (and parameter sidecar) under ``plots/`` in one command.

Usage::

    python make_plots.py                        # overhead, gate depth, summary
    python make_plots.py -- --eps-points 30     # forward extra flags to the search scripts
    python make_plots.py --multi-cap-overhead   # also write overhead_multi_cap.png
    python make_plots.py --all                  # + sanity checks + combined figure
                                                #   (combined needs error_analysis/exact-error.pdf)

Extra arguments after ``--`` are forwarded to each search-driven plotting
script (overhead and gate depth); the summary and the ``--all`` extras take no
search flags.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
_SCRIPTS = [
    _ROOT / "plotting" / "plot_overhead.py",
    _ROOT / "plotting" / "plot_gate_depth.py",
    # reads the two sidecars written above; must run last
    _ROOT / "plotting" / "plot_summary.py",
]
# Extra steps run only with ``--all``. They take no search flags, so the
# forwarded ``--`` arguments are not passed to them.
_EXTRA_SCRIPTS = [
    _ROOT / "plotting" / "plot_sanity_checks.py",
    _ROOT / "plotting" / "plot_combined_figure.py",
]
_EXACT_ERROR = _ROOT / "error_analysis" / "exact-error.pdf"
_FLAGS = ("--multi-cap-overhead", "--all")


def _extra_args(argv: list[str]) -> list[str]:
    if "--" in argv:
        idx = argv.index("--")
        return argv[idx + 1 :]
    return []


def main() -> int:
    argv = [a for a in sys.argv[1:] if a not in _FLAGS]
    multi_cap = "--multi-cap-overhead" in sys.argv[1:]
    run_all = "--all" in sys.argv[1:]
    extra = _extra_args(argv)

    for script in _SCRIPTS:
        print(f"\n=== Running {script.relative_to(_ROOT)} ===")
        cmd = [sys.executable, str(script)]
        if script.name != "plot_summary.py":
            cmd.extend(extra)
        if script.name == "plot_overhead.py" and multi_cap:
            cmd.extend(["--brute-bnorm-sq-caps", "10,100,1000"])
        rc = subprocess.call(cmd, cwd=_ROOT)
        if rc != 0:
            return rc

    if run_all:
        for script in _EXTRA_SCRIPTS:
            if script.name == "plot_combined_figure.py" and not _EXACT_ERROR.is_file():
                print(
                    f"\n=== Skipping {script.relative_to(_ROOT)}: "
                    f"{_EXACT_ERROR.relative_to(_ROOT)} not present "
                    "(see error_analysis/README.md) ==="
                )
                continue
            print(f"\n=== Running {script.relative_to(_ROOT)} ===")
            rc = subprocess.call([sys.executable, str(script)], cwd=_ROOT)
            if rc != 0:
                return rc

    print("\nAll plots regenerated in plots/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
