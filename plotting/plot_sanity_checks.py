"""Sanity checks for the extrapolation objective (Samson's list).

Objective actually optimized (section.tex Eq. 139):

    C(q) = (q_max/q_min) * lambda_scale^(1+1/p)
           * ( K * ||b^(p)||_1 / eps )^(1/(sigma(m-1)+p))

where ||b^(p)||_1 is the *suppressed* norm (Eq. 48 / line 50): the norm of the
vector with entries b_k * (q_1/q_k)^p and q_1 = q_min (the reference point, so
the ratio is <= 1 and the bound is *smaller* than the plain ||b||_1).

Checks implemented here:
  #1  suppressed norm used in objective, and matches the Eq.48 definition
  #3  forcing q_1 = 1 in the search is never better (worse or equal)
  #5  err -> 1: trivial schedule (m=1, q=[1], ||b||^2=1) is optimal; the gap to
      Trotter is then a pure constant factor

Stubbed (need Samson's exact line definition / image.png):
  #4  optimized line parallel to Trotter in the trivial regime
  #6  per-order lambda_comm bound line; all extrapolated data must lie below it

Run: python plotting/plot_sanity_checks.py
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import richardson as rt

ORDERS = (1, 2, 4)
Q_MAX = 15


def _tex_suppressed_norm(b, q, p, *, well_conditioned=False):
    """Eq. 48 / line 50 with q_1 = q_min:  Sum |b_k| (q_min/q_k)^e."""
    b = np.asarray(b, dtype=float)
    q = np.asarray(q, dtype=float)
    q_min = q.min()
    e = 1 if well_conditioned else int(p)
    return float(np.sum(np.abs(b * (q_min / q) ** e)))


def check1_suppressed_norm() -> None:
    print("#1  Suppressed norm - code vs section.tex Eq. 48 (q_1 = q_min)")
    print("    code  : b_suppressed_norm = Sum |b_k| (q_k/q_min)^p   [ratio >= 1]")
    print("    Eq.48 : ||b^(p)||_1        = Sum |b_k| (q_min/q_k)^p   [ratio <= 1]")
    for p in (2, 4):
        for m in (3, 5):
            _, q = rt.get_sample_points(m)
            b, _ = rt.get_richardson_coefficients([1.0 / x for x in q], p, m)
            plain = rt.b_norm1(b)
            code = rt.b_suppressed_norm(b, q, m, p)
            tex = _tex_suppressed_norm(b, q, p)
            flag = "  <-- code > plain (contradicts 'much smaller bounds')" if code > plain else ""
            print(f"    p={p} m={m}: plain={plain:8.3f}  code={code:9.3f}  Eq.48={tex:8.3f}{flag}")
    print()


def check3_force_q1_one() -> None:
    print("#3  Forcing q_1 = 1 (grid must contain 1) should never beat the free search")
    q_range = list(range(1, Q_MAX + 1))
    for p in ORDERS:
        for eps in (1e-2, 1e-4):
            best_free, best_forced = np.inf, np.inf
            for m in range(1, 7):
                for combo in itertools.combinations(q_range, m):
                    s = [1.0 / x for x in combo]
                    try:
                        b, _ = rt.get_richardson_coefficients(s, p, m)
                    except np.linalg.LinAlgError:
                        continue
                    if rt.b_norm1(b) ** 2 > 100.0:
                        continue
                    base, _, _ = rt._score_grid(
                        p=p, m=m, q_grid=list(combo), error=eps, A=1.0, n=1, b_list=b,
                    )
                    best_free = min(best_free, base)
                    if 1 in combo:
                        best_forced = min(best_forced, base)
            ok = "OK (worse/equal)" if best_forced >= best_free - 1e-12 else "VIOLATED"
            print(f"    p={p} eps={eps:.0e}: free={best_free:.4g}  q1=1 forced={best_forced:.4g}  -> {ok}")
    print()


def check5_trivial_at_eps1() -> None:
    print("#5  err -> 1 (largest eps): trivial schedule optimal, gap = constant factor")
    g = json.loads((_ROOT / "plots" / "gate_depth.params.json").read_text(encoding="utf-8"))
    eps = g["epsilon"]
    i = int(np.argmax(eps))
    for p in ORDERS:
        opt = g["results"][str(p)]["opt"]
        m, q, so = opt["m"][i], opt["q_grids"][i], opt["bnorm1_sq"][i]
        ratio = g["results"][str(p)]["opt"]["trotter_steps"][i] / opt["richardson_steps"][i]
        trivial = "OK" if (m == 1 and q == [1] and abs(so - 1) < 1e-9) else "NOT trivial"
        print(f"    p={p} eps={eps[i]:.3g}: m={m} q={q} ||b||^2={so:.3g} -> {trivial}; trot/rich={ratio:.4f}")
    print()


def check4_parallel_to_trotter() -> None:
    print("#4  [stub] optimized line parallel to Trotter in trivial regime - "
          "needs confirmation of intended regime/slope from Samson")
    print()


def check6_lambda_comm_line() -> None:
    print("#6  [stub] per-order lambda_comm*4 bound line; all extrapolated data below it")
    print("    Need Samson's exact line definition (image.png didn't come through).")
    print("    Best guess from tex: line = (q_max/q_min) * a * lambda_scale with a=4.")
    print()


def main() -> None:
    check1_suppressed_norm()
    check3_force_q1_one()
    check5_trivial_at_eps1()
    check4_parallel_to_trotter()
    check6_lambda_comm_line()


if __name__ == "__main__":
    main()
