"""The spin-harmonic axial rule of the ``altmann-herzig`` basis convention
(wayfinder ticket #21).

Altmann & Herzig write every two-dimensional irrep of an axial group in a
basis of two spherical or spin harmonics ``<|j m1>, s |j -m1>|`` (their
tables T n.6) and tabulate the Wigner matrices of the group's operations in
that basis (T n.7). This module holds the closed form of those matrices
on the generators of the rotation group ``R = D_N`` (or ``C_N``), as
``docs/spec/basis-conventions.md`` states it, exact in ``Q(zeta_M)`` with
``M = lcm(2N, 4)``:

* ``D(a) = diag(zeta_N^{-m1}, zeta_N^{m1})`` for the principal rotation
  ``a`` (``zeta_2N`` and an odd ``m1`` in halves for the projective irreps),
  where ``m1 = +k`` or ``-k`` is fixed by the chain rule below;
* ``D(b) = exp(-i pi m~) [[0, y], [conj(y), 0]]`` for the secondary rotation
  ``b``, with ``y = exp(-2 i m1 alpha)``, ``alpha`` the azimuth of the pole
  the reference gives the anchored secondary element (``0``, or ``pi/2``
  under the letter-preserving twist of ``Cnv`` with ``n = 2 (mod 4)``),
  and ``m~`` the residue of ``m1`` at the bottom of the chain.

Chain rule: let ``n0`` be the odd part of ``N`` (``2`` when ``N`` is a power
of two) and walk the levels ``n0, 2 n0, ..., N``; at the first level ``L``
where the residues of ``+k`` and ``-k`` modulo ``L`` (taken in
``(-L/2, L/2]``) differ, ``m1`` is the sign whose residue is positive; if
they agree on every level, ``m1 = +k``.

    uv run python tools/ah_axial_rule.py check   # the rule against the transcribed samples and kets
    uv run python tools/ah_axial_rule.py show N  # the generator images of every irrep of D_N

``check`` compares, entry by entry and exactly, (1) the rule's matrices with
the six transcribed axial sample tables of ticket #19
(``tests/reference/altmann-herzig/matrices/{D3,C4v,C6v,D3h,D2d,D6h}.toml``)
and (2) the Wigner matrices of the kets the book lists in T n.6 (``KETS``)
with the rule, which is how the rule was found. Nothing here is production
code.
"""

from __future__ import annotations

import math
import sys
import tomllib
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
import canonical_basis_prototype as cb

MAT_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "matrices"


# --------------------------------------------------------------------------
# The rule
# --------------------------------------------------------------------------


def residue(m2, L):
    """The representative of ``m2/2`` modulo ``L`` in ``(-L/2, L/2]``, in halves."""
    r = m2 % (2 * L)
    if r > L:
        r -= 2 * L
    return r


def chain(N):
    n0 = N
    while n0 % 2 == 0:
        n0 //= 2
    if n0 == 1:
        n0 = 2
    levels = []
    L = n0
    while L <= N:
        levels.append(L)
        L *= 2
    assert levels[-1] == N, (N, levels)
    return levels


def first_partner(N, k2):
    """``(m1, m~)`` in halves for the irrep with exponent ``k = k2/2`` on the principal element."""
    levels = chain(N)
    m2 = k2
    for L in levels:
        plus, minus = residue(k2, L), residue(-k2, L)
        if plus != minus:
            m2 = k2 if plus > 0 else -k2
            break
    return m2, residue(m2, levels[0])


def conductor(N):
    return math.lcm(2 * N, 4)


def zeta_pow(F, e):
    """``zeta_M^e`` as an element of ``F = Q(zeta_M)``."""
    return F.zeta(e % F.conductor)


def rule(N, k2, twisted, F=None):
    """The generator images ``(D(a), D(b))`` of the irrep with exponent ``k2/2``.

    ``k2`` is even for an ordinary irrep ``E_k`` (``k = k2/2``) and odd for a
    projective irrep ``E_{k2/2}``; ``twisted`` selects ``alpha = pi/2``.
    """
    F = F or cb.field(conductor(N))
    M = F.conductor
    assert M % (2 * N) == 0 and M % 4 == 0
    m2, mt2 = first_partner(N, k2)
    zero, one = F.zero(), F.rational(1)
    # zeta_2N^{-/+ m2} = zeta_M^{-/+ m2 * M/(2N)}
    Da = [
        [zeta_pow(F, -m2 * M // (2 * N)), zero],
        [zero, zeta_pow(F, m2 * M // (2 * N))],
    ]
    # exp(-i pi m~) = zeta_4^{-mt2} ; exp(-2 i m1 alpha) = zeta_4^{-m2} for alpha = pi/2
    t = zeta_pow(F, -mt2 * M // 4)
    y = zeta_pow(F, -m2 * M // 4) if twisted else one
    Db = [[zero, F.mul(t, y)], [F.mul(t, F.conj(y)), zero]]
    return Da, Db


def ket(N, j2, m2, s, bullet, twisted, F=None):
    """Wigner matrices of the generators in the basis ``<|j m1>, s |j -m1>|``.

    ``j2``, ``m2`` are ``2j``, ``2 m1``; ``bullet`` marks an ungerade spin
    harmonic (parity ``-1``); integral kets have parity ``(-1)^j``. Returns
    ``(D(a), D(b), parity)`` with ``D(b)`` the matrix of the *rotation* about
    the anchored pole (``alpha = 0`` or ``pi/2``).
    """
    F = F or cb.field(conductor(N))
    M = F.conductor
    zero = F.zero()
    Da = [
        [zeta_pow(F, -m2 * M // (2 * N)), zero],
        [zero, zeta_pow(F, m2 * M // (2 * N))],
    ]
    # R(alpha) |j m> = exp(-i pi j) exp(2 i m alpha) |j -m>
    ph = zeta_pow(F, -j2 * M // 4)
    y = zeta_pow(F, -m2 * M // 4) if twisted else F.rational(1)
    sgn = F.rational(s)
    Db = [[zero, F.mul(F.mul(ph, y), sgn)], [F.mul(F.mul(ph, F.conj(y)), sgn), zero]]
    if j2 % 2 == 0:
        parity = -1 if (j2 // 2) % 2 else 1
    else:
        parity = -1 if bullet else 1
    return Da, Db, parity


# --------------------------------------------------------------------------
# Samples: the transcribed generator matrices of ticket #19
# --------------------------------------------------------------------------

# group -> (N, twisted, book element -> role, book irrep -> k2)
# roles: ("c", p) = the principal element to the power p, "b" = the secondary
# element, "i" = the inversion (D(i) = +-1 by g/u)
SAMPLES = {
    "D3": (3, False, {"C3+": ("c", 1), "C2'1": "b"}, {"E": 2, "E1/2": 1}),
    "C4v": (4, False, {"C4+": ("c", 1), "sigma_v1": "b"}, {"E": 2, "E1/2": 1, "E3/2": 3}),
    "C6v": (
        6,
        True,
        {"C6+": ("c", 1), "sigma_v1": "b"},
        {"E1": 2, "E2": 4, "E1/2": 1, "E3/2": 3, "E5/2": 5},
    ),
    # D3h: the principal element is S3^-1 = i C6+, of order 6; C3+ is its square and
    # sigma_h its cube; E' and E'' are E_2 and E_1 of D6 (chi(S3) = -1 and +1)
    "D3h": (
        6,
        False,
        {"C3+": ("c", 2), "C2'1": "b", "sigma_h": ("c", 3)},
        {"E'": 4, "E''": 2, "E1/2": 1, "E3/2": 3, "E5/2": 5},
    ),
    # D2d: the principal element is S4^-1 = i C4+, the book's S4-
    "D2d": (4, False, {"S4-": ("c", 1), "C2'1": "b"}, {"E": 2, "E1/2": 1, "E3/2": 3}),
    "D6h": (
        6,
        False,
        {"C6+": ("c", 1), "C2'1": "b", "i": "i"},
        {
            "E1g": 2, "E2g": 4, "E1u": 2, "E2u": 4,
            "E1/2,g": 1, "E3/2,g": 3, "E5/2,g": 5, "E1/2,u": 1, "E3/2,u": 3, "E5/2,u": 5,
        },
    ),
}  # fmt: skip


def mpow(F, A, p):
    P = cb.midentity(F, len(A))
    for _ in range(p):
        P = cb.mmul(F, P, A)
    return P


def check_samples():
    bad = 0
    for group, (N, twisted, roles, kmap) in SAMPLES.items():
        data = tomllib.loads((MAT_DIR / f"{group}.toml").read_text())
        M = conductor(N)
        assert data["conductor"] == M, (group, data["conductor"], M)
        F = cb.field(M)
        problems = []
        for symbol, entry in data["irreps"].items():
            k2 = kmap[symbol]
            Da, Db = rule(N, k2, twisted, F)
            for g, A in zip(data["elements"], entry["matrices"]):
                want = [[cb.parse_value(F, x) for x in row] for row in A]
                role = roles[g]
                if role == "b":
                    got = Db
                elif role == "i":
                    sign = -1 if symbol.endswith("u") else 1
                    got = cb.mscale(F, cb.midentity(F, 2), F.rational(sign))
                else:
                    got = mpow(F, Da, role[1])
                if not cb.meq(got, want):
                    problems.append(f"{symbol} {g}")
        m1s = ", ".join(
            f"{s}: m1={Fr(first_partner(N, k2)[0], 2)}" for s, k2 in kmap.items()
        )
        print(
            f"{group} (N={N}{', twisted' if twisted else ''}): "
            f"{'ok' if not problems else 'MISMATCH'}  [{m1s}]"
        )
        for p in problems:
            print("  rule differs from the book on " + p)
        bad += bool(problems)
    return bad


# --------------------------------------------------------------------------
# Kets: the bases the book lists in T n.6 (first ket of the lowest j)
# --------------------------------------------------------------------------

# group -> (N, twisted, table, {irrep symbol: (k2, j2, m2, s, bullet)})
KETS = {
    "C4v": (4, False, "T 52.6, p. 490", {
        "E": (2, 2, 2, -1, False),
        "E1/2": (1, 1, 1, 1, False),
        "E3/2": (3, 3, -3, -1, False),
    }),
    "C6v": (6, True, "T 54.6, p. 498", {
        "E1": (2, 2, 2, -1, False),
        "E2": (4, 4, -4, -1, False),
        "E1/2": (1, 1, 1, 1, False),
        "E3/2": (3, 3, 3, 1, False),
        "E5/2": (5, 5, -5, 1, False),
    }),
}  # fmt: skip


def check_kets():
    bad = 0
    for group, (N, twisted, table, kets) in KETS.items():
        F = cb.field(conductor(N))
        problems = []
        for symbol, (k2, j2, m2, s, bullet) in kets.items():
            Da, Db = rule(N, k2, twisted, F)
            Ka, Kb, parity = ket(N, j2, m2, s, bullet, twisted, F)
            # Cnv: the secondary element is the plane i * b, so the ket's matrix carries its parity
            Kb = cb.mscale(F, Kb, F.rational(parity))
            if not (cb.meq(Da, Ka) and cb.meq(Db, Kb)):
                problems.append(symbol)
        print(f"{group} kets ({table}): {'ok' if not problems else 'MISMATCH'}")
        for p in problems:
            print("  ket and rule differ on " + p)
        bad += bool(problems)
    return bad


def show(N):
    F = cb.field(conductor(N))
    print(f"D_{N}: conductor {F.conductor}, chain {chain(N)}")
    for k2 in list(range(2, N, 2)) + list(range(1, N, 2)):
        m2, mt2 = first_partner(N, k2)
        Da, Db = rule(N, k2, False, F)
        name = f"E{k2 // 2}" if k2 % 2 == 0 else f"E{k2}/2"
        fmt = lambda A: [[cb.to_grammar(F, x) for x in row] for row in A]
        print(
            f"  {name:6s} m1 = {Fr(m2, 2)!s:5s} m~ = {Fr(mt2, 2)!s:4s} "
            f"D(a) = {fmt(Da)}  D(b) = {fmt(Db)}"
        )


def main(argv):
    if len(argv) >= 2 and argv[1] == "check":
        return 1 if check_samples() + check_kets() else 0
    if len(argv) >= 3 and argv[1] == "show":
        show(int(argv[2]))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
