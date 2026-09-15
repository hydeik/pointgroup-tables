"""The spin-harmonic axial rule of the ``altmann-herzig`` basis convention
(wayfinder tickets #21 and #22).

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
* ``D(b) = exp(-i pi m~) [[0, y], [conj(y), 0]]`` for the secondary element
  ``b``, with ``y = exp(-2 i m1 alpha)``, ``alpha`` the azimuth of the pole
  the reference gives the anchored secondary element (``0``, or ``pi/2``
  where the book's ``sigma_v1`` of ``Cnv`` lies in the ``yz`` plane, which
  is every ``n`` not divisible by four; the ordinary irreps of the odd
  ``Cnv`` groups take ``y = 1`` all the same, since the book prints for
  them the ``D_n`` matrices of ``C2'(k)``), and ``m~`` the residue of
  ``m1`` at the deciding level of the chain.

Chain rule: divide ``N`` by its smallest prime factor, then the quotient by
its smallest prime factor, and so on down to ``1``; the levels are the
quotients above ``1`` in ascending order (``2, 4, 8`` for ``N = 8``;
``3, 9`` for ``N = 9``; ``3, 6, 12`` for ``N = 12``; ``5, 10`` for
``N = 10``). At the first level ``L`` where the residues of ``+k`` and
``-k`` modulo ``L`` (taken in ``(-L/2, L/2]``) differ (they always differ
at ``N``), ``m1`` is the sign whose residue is positive, except that the
ordinary ``E_{N/4}`` with ``N/4`` odd and greater than one takes ``-k``.

    uv run python tools/ah_axial_rule.py check    # the rule against every transcribed table and ket
    uv run python tools/ah_axial_rule.py survey   # one line per printed ket: book vs rule
    uv run python tools/ah_axial_rule.py show N   # the generator images of every irrep of D_N

``check`` compares, entry by entry and exactly, (1) the rule's matrices with
the transcribed axial matrix tables
(``tests/reference/altmann-herzig/matrices/*.toml`` other than the polyhedral
ones) and (2) the Wigner matrices of every ket the book lists in T n.6
(``tests/reference/altmann-herzig/bases/*.toml``) with the rule, and (3) the
book's irrep label of each ket with the exponent the ket carries. Nothing
here is production code.
"""

from __future__ import annotations

import math
import re
import sys
import tomllib
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
import canonical_basis_prototype as cb

MAT_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "matrices"
BASES_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "bases"


# --------------------------------------------------------------------------
# The rule
# --------------------------------------------------------------------------


def residue(m2, L):
    """The representative of ``m2/2`` modulo ``L`` in ``(-L/2, L/2]``, in halves."""
    r = m2 % (2 * L)
    if r > L:
        r -= 2 * L
    return r


def smallest_prime_factor(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            return p
        p += 1
    return n


def chain(N):
    """The levels ``N, N/p1, N/(p1 p2), ...`` above ``1``, ascending (``p1 <= p2 <= ...``)."""
    levels = []
    L = N
    while L > 1:
        levels.append(L)
        L //= smallest_prime_factor(L)
    return levels[::-1]


def first_partner(N, k2):
    """``(m1, m~)`` in halves for the irrep with exponent ``k = k2/2`` on the principal element.

    ``m~`` is the residue of ``m1`` at the deciding level: the first level at
    which the residues of ``+k`` and ``-k`` differ (they always differ at
    ``N``). One exception to "the positive residue": the ordinary irrep
    ``E_{N/4}`` with ``N/4`` an odd number greater than one is decided at
    ``N`` only, and there the book takes ``m1 = -k`` (``D6d``, ``D10d``);
    with ``N/4`` even or ``N = 4`` it takes ``+k`` (``D8``, ``D8d``, ``D2d``).
    """
    for L in chain(N):
        plus, minus = residue(k2, L), residue(-k2, L)
        if plus != minus:
            m2 = k2 if plus > 0 else -k2
            if L == N and k2 % 2 == 0 and N % 4 == 0 and N > 4 and (N // 4) % 2 == 1:
                m2 = -k2
            return m2, residue(m2, L)
    raise ValueError((N, k2))


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
    ``(D(a), D(b), parity)`` with ``D(a)``, ``D(b)`` the matrices of the
    *rotations* ``C_N`` about ``z`` and ``C_2`` about the anchored pole
    (``alpha = 0`` or ``pi/2``); an improper generator is ``parity`` times
    the rotation's matrix.
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
# Families: which abstract generator the book's principal element and
# secondary element are, and how the book's irrep labels map to exponents
# --------------------------------------------------------------------------

GROUP_RE = re.compile(r"^([CD])(\d+)(v|h|d)?$")


def family(group):
    """``(kind, n, N)``: ``N`` is the order of the principal element of ``R = D_N``.

    kinds: ``Dn`` (proper), ``Cnv`` (``b`` is a plane), ``Dnh_even``
    (``D_n x Ci``), ``Dnh_odd`` (twisted ``D_2n``, ``a = S_n^-1 = i C_2n``),
    ``Dnd_even`` (twisted ``D_2n``, ``a = i C_2n``).
    """
    m = GROUP_RE.match(group)
    assert m, group
    letter, n, suffix = m.group(1), int(m.group(2)), m.group(3)
    if letter == "D" and suffix is None:
        return "Dn", n, n
    if letter == "C" and suffix == "v":
        return "Cnv", n, n
    if letter == "D" and suffix == "h":
        return ("Dnh_even", n, n) if n % 2 == 0 else ("Dnh_odd", n, 2 * n)
    if letter == "D" and suffix == "d":
        assert n % 2 == 0, group
        return "Dnd_even", n, 2 * n
    raise ValueError(group)


def pole_twisted(kind, n, setting=None):
    """``alpha = pi/2``: the book's ``sigma_v1`` of ``Cnv`` has pole ``(0, 1, 0)``
    unless ``4 | n`` (T n.1); setting B of ``C3v`` puts it at ``(1, 0, 0)``."""
    if kind != "Cnv":
        return False
    if setting == "B":
        return False
    return n % 4 != 0


def rule_twisted(kind, n, k2, setting=None):
    """The ``alpha`` the rule uses: the pole's azimuth, except that the ordinary
    irreps of ``Cnv`` with odd ``n`` take ``alpha = 0`` (the book prints for them
    the ``D_n`` matrices of ``C2'(k)`` under ``sigma_v(k) -> C2'(k)``, T n.7)."""
    if not pole_twisted(kind, n, setting):
        return False
    return n % 2 == 0 or k2 % 2 == 1


IRREP_RE = re.compile(r"^E(\d*)(/2)?(,g|,u|g|u|''|')?$")


def label_k2(kind, n, label):
    """The exponent ``k2`` (in halves) of the abstract irrep behind the book's label."""
    m = IRREP_RE.match(label)
    assert m, label
    num, half, tag = m.groups()
    k = int(num) if num else 1
    if kind in ("Dn", "Cnv", "Dnh_even"):
        return k if half else 2 * k
    if kind == "Dnh_odd":
        # E_k' and E_k'' are the irreps of D_2n with exponent k' = k or n - k,
        # the even one for a prime (chi(sigma_h) = +2), the odd one for two primes
        if half:
            return k
        want_even = tag == "'"
        return 2 * next(c for c in (k, n - k) if (c % 2 == 0) == want_even)
    if kind == "Dnd_even":
        # S_2n = a^(n+1) up to the parity of the harmonic: E_k of the book is E_k'
        # of D_2n with k' = k (n+1) mod 2n; the projective labels carry k directly
        if half:
            return k
        kk = (k * (n + 1)) % (2 * n)
        return 2 * min(kk, 2 * n - kk)
    raise ValueError(kind)


def label_parity(label):
    """``+1``/``-1`` for a ``g``/``u`` label, ``0`` otherwise."""
    if label.endswith("g"):
        return 1
    if label.endswith("u"):
        return -1
    return 0


def exponent_of(F, N, Da):
    """``k2`` in ``[0, N]`` and the signed ``m2`` of a diagonal ``D(a) = diag(zeta_2N^-m2, zeta_2N^m2)``."""
    M = F.conductor
    for e in range(M):
        if Da[0][0] == zeta_pow(F, e):
            break
    else:
        raise ValueError("D(a) is not a root of unity")
    assert (e * 2 * N) % M == 0, (e, N, M)
    m2 = (-e * 2 * N // M) % (2 * N)
    if m2 > N:
        m2 -= 2 * N
    return abs(m2), m2


# --------------------------------------------------------------------------
# Samples: the transcribed generator matrices (T n.7)
# --------------------------------------------------------------------------

# group -> (book element -> role); roles: ("c", p) = the principal element to
# the power p, "b" = the secondary element, "i" = the inversion (D(i) = +-1 by g/u)
SAMPLES = {
    "D3": {"C3+": ("c", 1), "C2'1": "b"},
    "D5": {"C5+": ("c", 1), "C2'1": "b"},
    "D8": {"C8+": ("c", 1), "C2'1": "b"},
    "D9": {"C9+": ("c", 1), "C2'1": "b"},
    "D10": {"C10+": ("c", 1), "C2'1": "b"},
    "C4v": {"C4+": ("c", 1), "sigma_v1": "b"},
    "C6v": {"C6+": ("c", 1), "sigma_v1": "b"},
    # D3h: the principal element is S3^-1 = i C6+, of order 6; C3+ is its square and
    # sigma_h its cube; E' and E'' are E_2 and E_1 of D6 (chi(S3) = -1 and +1)
    "D3h": {"C3+": ("c", 2), "C2'1": "b", "sigma_h": ("c", 3)},
    # D2d: the principal element is S4^-1 = i C4+, the book's S4-
    "D2d": {"S4-": ("c", 1), "C2'1": "b"},
    "D6h": {"C6+": ("c", 1), "C2'1": "b", "i": "i"},
}


def mpow(F, A, p):
    P = cb.midentity(F, len(A))
    for _ in range(p):
        P = cb.mmul(F, P, A)
    return P


def check_samples():
    bad = 0
    for group, roles in SAMPLES.items():
        data = tomllib.loads((MAT_DIR / f"{group}.toml").read_text())
        kind, n, N = family(group)
        twisted = pole_twisted(kind, n)
        M = conductor(N)
        assert data["conductor"] == M, (group, data["conductor"], M)
        F = cb.field(M)
        problems = []
        kmap = {}
        for symbol, entry in data["irreps"].items():
            k2 = label_k2(kind, n, symbol)
            kmap[symbol] = k2
            Da, Db = rule(N, k2, rule_twisted(kind, n, k2), F)
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
            f"{group} matrices (N={N}{', alpha=pi/2' if twisted else ''}): "
            f"{'ok' if not problems else 'MISMATCH'}  [{m1s}]"
        )
        for p in problems:
            print("  rule differs from the book on " + p)
        bad += bool(problems)
    return bad


# --------------------------------------------------------------------------
# Kets: the bases the book lists in T n.6
# --------------------------------------------------------------------------


def load_bases():
    out = []
    for path in sorted(BASES_DIR.glob("*.toml")):
        data = tomllib.loads(path.read_text())
        out.append((path.stem, data))
    return out


def ket_rows(name, data):
    """Yield ``(label, row, Ka, Kb, parity, k2_ket, m2_ket)`` for every printed ket.

    ``Ka`` and ``Kb`` are the matrices of the group's generators ``a`` and
    ``b`` (parity applied where the generator is improper).
    """
    kind, n, N = family(data["group"])
    twisted = pole_twisted(kind, n, data.get("setting"))
    F = cb.field(conductor(N))
    for label, entry in data["irreps"].items():
        for row in entry["kets"]:
            j, m1 = Fr(row["j"]), Fr(row["m1"])
            j2, m2 = int(2 * j), int(2 * m1)
            assert j2 % 2 == m2 % 2 and abs(m2) <= j2, (name, label, row)
            Ka, Kb, parity = ket(N, j2, m2, row["s"], row["bullet"], twisted, F)
            if kind == "Cnv":
                Kb = cb.mscale(F, Kb, F.rational(parity))
            elif kind in ("Dnh_odd", "Dnd_even"):
                Ka = cb.mscale(F, Ka, F.rational(parity))
            k2, m2k = exponent_of(F, N, Ka)
            yield F, kind, n, N, twisted, label, row, Ka, Kb, parity, k2, m2k


def check_kets():
    bad = 0
    for name, data in load_bases():
        problems = []
        count = 0
        for F, kind, n, N, twisted, label, row, Ka, Kb, parity, k2, m2k in ket_rows(
            name, data
        ):
            count += 1
            where = f"{label} <|{row['j']} {row['m1']}>, {'-' if row['s'] < 0 else ''}|{row['j']} -({row['m1']})>|{'*' if row['bullet'] else ''}"
            if kind == "Dnh_even" and parity != label_parity(label):
                problems.append(f"{where}: parity of the harmonic is not the label's g/u")
            want_k2 = label_k2(kind, n, label)
            if k2 != want_k2:
                problems.append(
                    f"{where}: the ket carries exponent {Fr(k2, 2)}, the label says {Fr(want_k2, 2)}"
                )
            Da, Db = rule(N, k2, rule_twisted(kind, n, k2, data.get("setting")), F)
            if not cb.meq(Da, Ka):
                problems.append(
                    f"{where}: D(a) differs (book m1 = {Fr(m2k, 2)}, rule m1 = {Fr(first_partner(N, k2)[0], 2)})"
                )
            elif not cb.meq(Db, Kb):
                problems.append(f"{where}: D(b) differs (sign or phase)")
        kind, n, N = family(data["group"])
        print(
            f"{name} kets ({data['table']}, N={N}{', pole at pi/2' if pole_twisted(kind, n, data.get('setting')) else ''}): "
            f"{count} kets, {'ok' if not problems else 'MISMATCH'}"
        )
        for p in problems:
            print("  " + p)
        bad += bool(problems)
    return bad


def survey():
    """One line per printed ket: what the book chose and what the rule gives."""
    print("group   N  label     j     m1    s  par  | rule m1  m~  | a  b")
    for name, data in load_bases():
        for F, kind, n, N, twisted, label, row, Ka, Kb, parity, k2, m2k in ket_rows(
            name, data
        ):
            Da, Db = rule(N, k2, rule_twisted(kind, n, k2, data.get("setting")), F)
            m2, mt2 = first_partner(N, k2)
            print(
                f"{name:6s} {N:3d}  {label:8s} {row['j']:>5s} {row['m1']:>5s} {row['s']:+2d}  {parity:+2d}"
                f"{'*' if row['bullet'] else ' '} | {str(Fr(m2, 2)):>6s} {str(Fr(mt2, 2)):>4s} | "
                f"{'ok' if cb.meq(Da, Ka) else 'XX'} {'ok' if cb.meq(Db, Kb) else 'XX'}"
            )


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
    if len(argv) >= 2 and argv[1] == "survey":
        survey()
        return 0
    if len(argv) >= 3 and argv[1] == "show":
        show(int(argv[2]))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
