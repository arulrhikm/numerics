"""Diagnostic: WC vs optimized Richardson b-norms on the LKW grid.

WC extrapolation (Lemma 5 / LKW) is prescriptive: even error series starting at
s^2, Vandermonde rows [1, s_k^2, s_k^4, ...], closed-form b via Lagrange
weights. Same b for every even-p application (p=2 and p=4 share one matrix).

Optimized extrapolation is application-specific: Vandermonde rows use the
leading error powers of order-p Suzuki formulas (p, p+sigma, p+2sigma, ...).

Run: python check_bnorms_wellconditioned.py
"""

from __future__ import annotations

import numpy as np

import richardson as rt


M_VALUES = list(range(1, 11))
EVEN_P = (2, 4)
OPT_P = (1, 2, 4)


def _b_via_matrix_solve(s_list, m: int) -> np.ndarray:
    V = rt.setup_wc_vandermonde_matrix(s_list, m)
    e1 = np.zeros(m)
    e1[0] = 1.0
    return np.linalg.solve(V, e1)


def main() -> None:
    print("WC (prescriptive even-power Vandermonde, closed-form b)")
    print("  rows: [1, s_k^2, s_k^4, ..., s_k^(2(m-1))]")
    print("  q_k = ceil((sqrt(8) m / pi) / sin(pi(2k-1)/(8m)))")
    print(f"  m in {M_VALUES}\n")

    header = f"{'m':>2}  {'q_k':52s}  ||b||_1 (WC)  cond(V)  max|closed-solve|"
    print(header)
    print("-" * len(header))

    for m in M_VALUES:
        s_list, q_grid = rt.get_sample_points(m)
        q_str = "[" + ", ".join(str(q) for q in q_grid) + "]"
        b_closed, V = rt.get_wc_richardson_coefficients(s_list, m)
        b_solve = _b_via_matrix_solve(s_list, m)
        diff = float(np.max(np.abs(b_closed - b_solve)))
        cond = float(np.linalg.cond(V))
        print(
            f"{m:>2}  {q_str:52s}  {rt.b_norm1(b_closed):12.4f}  "
            f"{cond:9.3e}  {diff:12.2e}"
        )

    print()
    print("Samson example m=5, q=[58, 20, 12, 9, 7]: WC b is identical for p=2 and p=4")
    s_list, q_grid = rt.get_sample_points(5)
    assert q_grid == [58, 20, 12, 9, 7], q_grid
    b_wc, _ = rt.get_wc_richardson_coefficients(s_list, 5)
    for p in EVEN_P:
        b_opt, V_opt = rt.get_richardson_coefficients(s_list, p, 5)
        wc_powers = [2 * j for j in range(5)]
        opt_powers = [0] + [p + rt.sigma_parity(p) * j for j in range(4)]
        same_v = np.allclose(
            rt.setup_wc_vandermonde_matrix(s_list, 5),
            rt.setup_vandermonde_matrix(s_list, p, 5),
        ) if p == 2 else False
        print(
            f"  p={p}: ||b||_1 WC={rt.b_norm1(b_wc):.4f}  "
            f"opt={rt.b_norm1(b_opt):.4f}  "
            f"WC==opt Vandermonde? {same_v if p == 2 else 'n/a (p=4 opt uses [0,4,6,8,10])'}"
        )
        print(f"         WC row powers:  {wc_powers}")
        print(f"         opt row powers: {opt_powers}")

    print()
    print("Optimized extrapolation ||b||_1 (p-specific Vandermonde, matrix solve)")
    print(f"{'m':>2}  {'q_k':52s}", end="")
    for p in OPT_P:
        print(f"  ||b||_1 (p={p})", end="")
    print()
    print("-" * (2 + 2 + 52 + len(OPT_P) * 14))

    for m in M_VALUES:
        s_list, q_grid = rt.get_sample_points(m)
        q_str = "[" + ", ".join(str(q) for q in q_grid) + "]"
        row = f"{m:>2}  {q_str:52s}"
        for p in OPT_P:
            b, _ = rt.get_richardson_coefficients(s_list, p, m)
            row += f"  {rt.b_norm1(b):12.4f}"
        print(row)

    print()
    print("Note: p=1 must not use WC (odd error series). Plots skip WC for p=1.")


if __name__ == "__main__":
    main()
