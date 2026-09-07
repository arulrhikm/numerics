"""Numerical checks for the diagonal-restriction derivation.

Verifies the computable claims of the diagonal-restriction derivation
(`../notes/derivation/derivation-tighter-bounds.tex` in the working tree):

  1. `ratio_recall_sup`   — the supremum `R_p` of Eq. (lcomm-tight-recall) and the
     gate-depth ratio `R_p^(1+1/p)` of Eq. (ratio-recall), reconciled against the
     numbers quoted in the .tex/paper and the table in `richardson.py`. Surfaces
     the `richardson.py` LAMBDA-scale table bug (p=1,4 store the raised ratio
     where the bare sup is wanted).
  2. `pinching_demo`      — builds a small Hermitian `H = sum_g H_g`, pinches nested
     commutators onto its (block-)eigenbasis, and checks the contraction
     `||D(C)|| <= ||C||` (Lemma pinching-contraction), the vanishing of a purely
     off-diagonal commutator under `D`, and the strict `alpha^{D,(j)} < alpha^(j)`
     (Lemma acomm-D-smaller).
  3. `p0_table`           — the truncation order `p_0(eps) = ceil(log(2n/eps))` of
     the finitely-many-orders section, for `n = 100` and `eps` down to `1e-6`.
  4. `collapse_check`     — the row-(b) collapse condition (eq:collapse-condition).

numpy only. Run (from the repo root): python checks/verify_tighter_bounds.py
(`python`, not `python3`, on Windows.)
"""

from __future__ import annotations

import itertools
import math

import numpy as np

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import richardson as rt

# Values quoted in the .tex prose / paper.tex after Lemma lcomm-tight-bound.
# These are the *ratio* R_p^(1+1/p), not the bare sup R_p.
TEX_RATIO_QUOTED = {1: 1.5035, 2: 1.1487, 4: 1.0445}


# --------------------------------------------------------------------------- #
# 1. The ratio sup R_p and the gate-depth ratio R_p^(1+1/p)
# --------------------------------------------------------------------------- #
def ratio_recall_sup(p: int) -> tuple[float, float]:
    """Return `(R_p, R_p^(1+1/p))` for Trotter order `p`.

    `R_p = sup_{0<k<=1/p} ( e/(p+1)^2 (1/k - p) )^k` is the bare geometric ratio
    of Lemma lcomm-tight-bound divided by `n` (the `lambda_scale`); the quantity
    that actually multiplies the gate depth is `R_p^(1+1/p)`.
    """
    bare = rt.lemma57_geometric_ratio(p)
    return bare, bare ** (1.0 + 1.0 / p)


def report_ratios() -> None:
    print("[1] Gate-depth ratio  lambda_comm^(1+1/p) / (alpha_comm^(p+1))^(1/p)\n")
    header = (
        f"{'p':>2}  {'R_p (bare sup)':>14}  {'R_p^(1+1/p) (ratio)':>20}  "
        f"{'tex/paper quotes':>16}  {'richardson.py table':>20}"
    )
    print(header)
    print("-" * len(header))
    for p in (1, 2, 4):
        bare, ratio = ratio_recall_sup(p)
        table = rt.LEMMA57_GEOMETRIC_RATIO_BY_P[p]
        print(
            f"{p:>2}  {bare:>14.4f}  {ratio:>20.4f}  "
            f"{TEX_RATIO_QUOTED[p]:>16.4f}  {table:>20.4f}"
        )
    print()
    print("  - The tex/paper numbers {1.5035, 1.1487, 1.0445} are the RATIO R_p^(1+1/p),")
    print("    but the identity was written with a bare sup on the RHS: the exponent")
    print("    (1+1/p) was dropped. The numbers are right; the formula needs the power.")
    print("  - richardson.py's LAMBDA-scale table is used as the BARE R_p (it is then")
    print("    raised to (1+1/p) in the step counts). It stores the bare value only at")
    print("    p=2; at p=1,4 it stores R_p^(1+1/p), double-counting the exponent and")
    print("    inflating the plotted extrapolated-Trotter curves. Correct bare values:")
    print(f"    p=1 -> {ratio_recall_sup(1)[0]:.4f} (table {rt.LEMMA57_GEOMETRIC_RATIO_BY_P[1]:.4f}),  "
          f"p=4 -> {ratio_recall_sup(4)[0]:.4f} (table {rt.LEMMA57_GEOMETRIC_RATIO_BY_P[4]:.4f}).")
    print()


# --------------------------------------------------------------------------- #
# 2. Pinching contraction on a small Hamiltonian
# --------------------------------------------------------------------------- #
I2 = np.eye(2, dtype=complex)
PAULI = {
    "I": I2,
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_string(label: str) -> np.ndarray:
    """Kronecker product of single-qubit Paulis, e.g. 'XIZ'."""
    out = np.array([[1]], dtype=complex)
    for ch in label:
        out = np.kron(out, PAULI[ch])
    return out


def random_hermitian(d: int, rng: np.random.Generator) -> np.ndarray:
    """A generic dense d x d Hermitian matrix."""
    M = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    return M + M.conj().T


def spec_norm(X: np.ndarray) -> float:
    """Operator (spectral) norm."""
    return float(np.linalg.norm(X, 2))


def pinch(X: np.ndarray, H: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Block-diagonal pinching D(X) = sum_n Pi_n X Pi_n onto H's eigenbasis.

    Eigenvalues within `tol` are treated as one degenerate level, so this is the
    block-diagonal form of Eq. (pinching-def), not the strict diagonal.
    """
    w, V = np.linalg.eigh(H)
    Xe = V.conj().T @ X @ V
    same_level = np.abs(w[:, None] - w[None, :]) < tol
    return V @ (Xe * same_level) @ V.conj().T


def nested_commutator(ops: list[np.ndarray]) -> np.ndarray:
    """[A1, [A2, ..., [A_{j-1}, Aj] ... ]] for ops = [A1, ..., Aj]."""
    result = ops[-1]
    for A in reversed(ops[:-1]):
        result = A @ result - result @ A
    return result


def alpha_comm(fragments: list[np.ndarray], j: int, *, H=None) -> float:
    """Order-j commutator factor; dephased (alpha^{D,(j)}) if H is given."""
    total = 0.0
    for combo in itertools.product(fragments, repeat=j):
        C = nested_commutator(list(combo))
        if H is not None:
            C = pinch(C, H)
        total += spec_norm(C)
    return total


def pinching_demo() -> None:
    print("[2] Pinching contraction on H = sum_g H_g\n")
    # Generic dense Hermitian fragments (Gamma = 3): no accidental symmetry, so
    # the spectrum is nondegenerate (D = strict diagonal) and the dephased
    # factors are strictly smaller but nonzero -- the honest large-Gamma picture.
    rng = np.random.default_rng(0)
    d = 8
    fragments = [random_hermitian(d, rng) for _ in range(3)]
    H = sum(fragments)

    # (a) contraction ||D(C)|| <= ||C|| on every nested commutator up to order 3.
    max_ratio = 0.0
    for j in (2, 3):
        for combo in itertools.product(fragments, repeat=j):
            C = nested_commutator(list(combo))
            nC = spec_norm(C)
            nDC = spec_norm(pinch(C, H))
            assert nDC <= nC + 1e-9, (j, nDC, nC)
            if nC > 1e-12:
                max_ratio = max(max_ratio, nDC / nC)
    print(f"  (a) ||D(C)|| <= ||C|| held for all nested commutators; "
          f"max observed ||D(C)||/||C|| = {max_ratio:.4f}  (<= 1).")

    # (b) a purely off-diagonal commutator pinches to zero: with H diagonal, its
    #     eigenbasis is the computational basis, and any zero-diagonal operator
    #     is annihilated by D.
    Hdiag = np.diag(np.arange(1.0, d + 1.0)).astype(complex)
    C_off = random_hermitian(d, rng)
    np.fill_diagonal(C_off, 0.0)  # purely off-diagonal, still Hermitian
    n_off = spec_norm(pinch(C_off, Hdiag))
    assert n_off < 1e-9, n_off
    print(f"  (b) purely off-diagonal C under diagonal H: ||D(C)|| = {n_off:.2e}  (= 0).")

    # (c) exact Gamma=2 fact: for a two-fragment H = A + B, [A,B] = [A,H], and D
    #     annihilates [H, .], so the whole order-2 dephased factor vanishes.
    A, B = fragments[0], fragments[1]
    d2 = spec_norm(pinch(A @ B - B @ A, A + B))
    assert d2 < 1e-9, d2
    print(f"  (c) two-fragment H=A+B: ||D([A,B])|| = {d2:.2e}  (= 0, since [A,B]=[A,H]).")

    # (d) assembled alpha^{D,(j)} < alpha^(j), strictly, for Gamma = 3.
    print(f"\n  {'j':>2}  {'alpha^(j)':>14}  {'alpha^{D,(j)}':>16}  {'ratio':>8}")
    print("  " + "-" * 44)
    for j in (2, 3):
        a = alpha_comm(fragments, j)
        aD = alpha_comm(fragments, j, H=H)
        assert 0.0 < aD < a - 1e-9, (j, aD, a)
        print(f"  {j:>2}  {a:>14.4f}  {aD:>16.4f}  {aD / a:>8.4f}")
    print("\n  => 0 < alpha^{D,(j)} < alpha^(j) strictly (Lemma acomm-D-smaller), so lambda^D < lambda.\n")


# --------------------------------------------------------------------------- #
# 3. Truncation order p_0(eps) = ceil(log(2n/eps))
# --------------------------------------------------------------------------- #
def p0(eps: float, n: int) -> int:
    """Truncation order ceil(log(2n/eps)) (natural log, per the k-local lemma)."""
    return math.ceil(math.log(2 * n / eps))


def p0_table(n: int = 100) -> None:
    print(f"[3] BCH truncation order  p_0(eps) = ceil(log(2n/eps))  at n = {n}\n")
    print(f"  {'eps':>8}  {'p_0':>4}")
    print("  " + "-" * 14)
    for k in range(1, 7):
        eps = 10.0 ** (-k)
        print(f"  {eps:>8.0e}  {p0(eps, n):>4d}")
    # extrapolation variant p_0(eps') with eps' = s_hat eps / (2||b||_1)
    b_norm1, s_hat, eps = 12.0, 1.0 / 15.0, 1e-6
    p0_ext = math.ceil(math.log(4 * n * b_norm1 / (s_hat * eps)))
    print(f"\n  extrapolation variant (||b||_1={b_norm1}, s_hat={s_hat:.3f}, eps={eps:.0e}): "
          f"p_0(eps') = {p0_ext}")
    print("  => a list of order tens, not infinitely many orders.\n")


# --------------------------------------------------------------------------- #
# 4. Row-(b) collapse condition
# --------------------------------------------------------------------------- #
def collapse_check(alpha_p1: float, Lambda: float, T: float, p: int) -> bool:
    """True iff (alpha^(p+1))^(1/p) T^(1/p) <= Lambda (eq:collapse-condition)."""
    return (alpha_p1 ** (1.0 / p)) * (T ** (1.0 / p)) <= Lambda


def report_collapse() -> None:
    print("[4] Row-(b) collapse condition  (alpha^(p+1))^(1/p) T^(1/p) <= Lambda\n")
    Lambda, T, p = 50.0, 10.0, 2
    print(f"  {'alpha^(p+1)':>12}  {'collapses to Lambda*T?':>22}  {'depth term':>18}")
    print("  " + "-" * 56)
    for alpha_p1 in (1.0, 25.0, 250.0, 2500.0):
        collapses = collapse_check(alpha_p1, Lambda, T, p)
        depth = "Lambda*T" if collapses else "(alpha^(p+1))^(1/p) T^(1+1/p)"
        print(f"  {alpha_p1:>12.1f}  {str(collapses):>22}  {depth:>28}")
    print(f"\n  (Lambda={Lambda}, T={T}, p={p}: collapse holds while alpha^(p+1) <= "
          f"{Lambda ** p / T:.0f}.)\n")


# --------------------------------------------------------------------------- #
# 5. Non-multiplicativity of the pinching, and the first-order shift
# --------------------------------------------------------------------------- #
def eigenprojectors(H: np.ndarray, tol: float = 1e-9) -> list[np.ndarray]:
    """Eigenprojectors Pi_n of H, one per distinct eigenvalue (degeneracies merged)."""
    w, V = np.linalg.eigh(H)
    projectors = []
    start = 0
    for k in range(1, len(w) + 1):
        if k == len(w) or w[k] - w[start] > tol:
            block = V[:, start:k]
            projectors.append(block @ block.conj().T)
            start = k
    return projectors


def nonmultiplicativity_demo() -> None:
    """Checks the claims of `conjecture-1-explained.md` section 8.

    (a) D(AB) != D(A)D(B), and the difference is exactly the k != n cross terms.
    (b) X purely off-diagonal => D(X) = 0 but D(X @ X) != 0: off-diagonal times
        off-diagonal returns to the diagonal. This is the operator-level shadow of
        the second-order Dyson term sum_k <m|E|k><k|E|n> at m = n, k != n.
    (c) V = i[H_1, H] is Hermitian with D(V) = 0, and the spectrum of H + eps*V
        matches that of H to O(eps^2) -- i.e. no first-order energy shift. The
        vanishing shift IS the statement D(V) = 0, since first-order perturbation
        theory reads dE_n = <n|V|n>, the diagonal part.

    Its own rng, so sections [1]-[4] print unchanged.
    """
    print("[5] Non-multiplicativity of the pinching  (conjecture-1-explained.md section 8)\n")
    rng = np.random.default_rng(1234)
    d = 8

    # (a) D(AB) vs D(A)D(B), against the explicit cross-term formula.
    A, B = random_hermitian(d, rng), random_hermitian(d, rng)
    H = A + B
    lhs = pinch(A @ B, H)
    rhs = pinch(A, H) @ pinch(B, H)
    gap = spec_norm(lhs - rhs)
    # D(AB) = sum_{n,k} Pi_n A Pi_k B Pi_n -- insert a resolution of the identity.
    projs = eigenprojectors(H)
    cross = sum(Pn @ A @ Pk @ B @ Pn for Pn in projs for Pk in projs)
    assert spec_norm(lhs - cross) < 1e-9, spec_norm(lhs - cross)
    assert gap > 1e-6, gap
    print(f"  (a) ||D(AB) - D(A)D(B)||                 = {gap:.4e}   (!= 0)")
    print(f"      ||D(AB) - sum_{{n,k}} Pi_n A Pi_k B Pi_n|| = {spec_norm(lhs - cross):.2e}   (= 0)")
    print("      => the difference is exactly the k != n cross terms.")

    # (b) off-diagonal squared returns to the diagonal.
    Hdiag = np.diag(np.arange(1.0, d + 1.0)).astype(complex)
    X = random_hermitian(d, rng)
    np.fill_diagonal(X, 0.0)  # purely off-diagonal, still Hermitian
    nDX, nDX2 = spec_norm(pinch(X, Hdiag)), spec_norm(pinch(X @ X, Hdiag))
    assert nDX < 1e-9, nDX
    assert nDX2 > 1e-6, nDX2
    print(f"\n  (b) X purely off-diagonal:  ||D(X)|| = {nDX:.2e}  (= 0),  "
          f"||D(X@X)|| = {nDX2:.4e}  (!= 0)")
    print("      => off-diagonal x off-diagonal feeds back onto the diagonal.")

    # (c) V = i[H_1, H] is Hermitian, D(V) = 0, and gives no first-order shift.
    #     Reseed until H's minimum gap comfortably exceeds the largest perturbation,
    #     so sorted-spectrum comparison cannot swap levels and fake a shift.
    while True:
        H1, H2 = random_hermitian(d, rng), random_hermitian(d, rng)
        H = H1 + H2
        w0 = np.linalg.eigvalsh(H)
        if np.min(np.diff(w0)) > 0.5:
            break
    V = 1j * (H1 @ H - H @ H1)
    assert spec_norm(V - V.conj().T) < 1e-9, "V must be Hermitian"
    assert spec_norm(pinch(V, H)) < 1e-9, spec_norm(pinch(V, H))
    print(f"\n  (c) V = i[H_1,H]:  Hermitian, ||D(V)|| = {spec_norm(pinch(V, H)):.2e}  (= 0)")
    print(f"\n      {'eps':>8}  {'max|dE|':>12}  {'drop vs prev':>13}")
    print("      " + "-" * 37)
    prev = None
    for eps in (1e-2, 1e-3, 1e-4):
        shift = float(np.max(np.abs(np.linalg.eigvalsh(H + eps * V) - w0)))
        drop = "" if prev is None else f"{prev / shift:>13.1f}x"
        if prev is not None:
            # quadratic scaling => ~100x per decade; linear would be ~10x.
            assert 30.0 < prev / shift < 300.0, (eps, prev / shift)
        print(f"      {eps:>8.0e}  {shift:>12.3e}  {drop:>13}")
        prev = shift
    print("\n      => ~100x per decade of eps: the shift is O(eps^2), so there is no")
    print("         first-order energy shift. That vanishing IS D(V) = 0, since")
    print("         first-order perturbation theory reads dE_n = <n|V|n>.\n")


def main() -> None:
    report_ratios()
    pinching_demo()
    p0_table()
    report_collapse()
    nonmultiplicativity_demo()
    print("All assertions passed.")


if __name__ == "__main__":
    main()
