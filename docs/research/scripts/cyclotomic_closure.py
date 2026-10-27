"""Cyclotomic closure of chain-adapted bases and Clebsch-Gordan coefficients.

Computational evidence for ``docs/research/cyclotomic-closure.md`` (wayfinder
research ticket #40). Run with ``uv run python3
docs/research/scripts/cyclotomic_closure.py [GROUP ...]``; without arguments
the default set of groups is surveyed. The script uses only the standard
library and the exact arithmetic of ``pointgroup_tables.field``.

What it computes
----------------

For every irrep pair ``(alpha, beta)`` of a group and every irrep ``gamma`` of
the product, the projection method gives the unnormalised Clebsch-Gordan
vector ``w = Q e_j`` where ``Q`` is the Hermitian projector onto the
``gamma``-component (first partner), restricted to the complement of the
already found copies; its designated component is real and equals its squared
norm ``N = w_j = <w|w>``. The normalised vector with that component real and
positive is ``w / sqrt(N)``, so every Clebsch-Gordan coefficient of the pair
lies in ``K(sqrt(N_1), sqrt(N_2), ...)``, ``K = Q(zeta_M)`` the document
field. The same holds for the chain-adapted vectors ``P e_j`` of a group chain
(``P`` the product of the isotypic projectors of the chain's subgroups).

For every squared norm ``N`` (a totally positive element of ``K``) the script
decides exactly:

- whether ``sqrt(N)`` lies in ``K``;
- otherwise whether ``K(sqrt(N))`` is abelian over ``Q`` (Galois-theoretic
  criterion, decided with square roots in ``K`` only), i.e. whether
  ``sqrt(N)`` lies in any cyclotomic field at all (Kronecker-Weber);
- if so, the smallest multiple ``M'`` of ``M`` with ``sqrt(N)`` in
  ``Q(zeta_M')`` (the conductor of ``K(sqrt(N))``), found by reducing
  membership in ``Q(zeta_M')`` to square roots in ``K`` through the
  quadratic subextensions of ``Q(zeta_M')/K``;
- whether ``N`` is a Hermitian norm ``mu * conj(mu)`` from ``K`` itself (a
  unit-modulus rephasing then puts the vector into ``K``), by a search over
  products of short cyclotomic integers.

Every positive answer is verified exactly in the field; a negative answer to
the first three questions is a proof (the sign enumeration of the square-root
search is exhaustive); a negative answer to the fourth only says that the
search did not find a factor.
"""

from __future__ import annotations

import cmath
import math
import random
import sys
import time
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cache
from itertools import combinations, product

from pointgroup_tables.document import Model, build_model
from pointgroup_tables.field import (
    CyclotomicField,
    Matrix,
    Value,
    cyclotomic_field,
    euler_phi,
)
from pointgroup_tables.irreps import Irrep
from pointgroup_tables.symbols import parse_group
from pointgroup_tables.values import pretty

DEFAULT_GROUPS = [
    "C3v",
    "D4h",
    "D5h",
    "D6h",
    "D7h",
    "D8h",
    "D9d",
    "T",
    "Td",
    "O",
    "Oh",
    "I",
    "Ih",
]
# the chains probed per rotation group: (subgroup name, how its generators are found)
CHAINS: dict[str, list[list[tuple[str, str]]]] = {
    "T": [[("C3", "axis:C3^1(1)")], [("D2", "d2"), ("C2", "c2")]],
    "O": [
        [("D4", "dihedral:C4^1(1)"), ("C4", "axis:C4^1(1)")],
        [("T", "tetrahedral"), ("C3", "axis:C3^1(1)")],
        [("D3", "dihedral:C3^1(1)"), ("C3", "axis:C3^1(1)")],
        [("T", "tetrahedral"), ("D2", "d2"), ("C2", "c2")],
    ],
    "I": [
        [("D5", "dihedral:C5^1(1)"), ("C5", "axis:C5^1(1)")],
        [("T", "tetrahedral"), ("C3", "axis:C3^1(1)")],
        # I > D3 > C3 is not canonical (H and I5/2 restrict to D3 with multiplicity
        # two, found by the rank check) and is left out of the default survey
        [("T", "tetrahedral"), ("D2", "d2"), ("C2", "c2")],
    ],
}
MAX_DEGREE = 100  # the largest phi(M') searched for sqrt(N) in Q(zeta_M')
TWO_PI = 2 * math.pi


# ---- exact square roots in Q(zeta_M) ----------------------------------------------


def embed(F: CyclotomicField, x: Value, k: int = 1) -> complex:
    """The complex embedding ``zeta -> exp(2 pi i k / M)`` of ``x``."""
    M = F.conductor
    return sum(
        (float(c) * cmath.exp(1j * TWO_PI * j * k / M) for j, c in enumerate(x) if c),
        0j,
    )


def galois(F: CyclotomicField, x: Value, k: int) -> Value:
    """The Galois conjugate ``sigma_k: zeta -> zeta^k`` of ``x`` (``gcd(k, M) = 1``)."""
    return F.total(F.scale(F.zeta(j * k), c) for j, c in enumerate(x) if c)


def units(M: int) -> list[int]:
    return [k for k in range(1, M + 1) if math.gcd(k, M) == 1]


def _denominator_lcm(x: Value) -> int:
    out = 1
    for c in x:
        out = math.lcm(out, c.denominator)
    return out


@cache
def _inverse_vandermonde(M: int) -> list[list[complex]]:
    """The inverse of ``V[k][j] = exp(2 pi i k j / M)`` over ``k`` in ``(Z/M)^*``,
    ``j < phi(M)``: the map from the embeddings of an element to its coefficients."""
    ks = units(M)
    n = len(ks)
    rows = [
        [cmath.exp(1j * TWO_PI * k * j / M) for j in range(n)]
        + [1.0 + 0j if r == c else 0j for c in range(n)]
        for r, k in enumerate(ks)
    ]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(rows[r][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        p = rows[col][col]
        rows[col] = [v / p for v in rows[col]]
        for r in range(n):
            if r != col:
                f = rows[r][col]
                if f:
                    rows[r] = [
                        a - f * b for a, b in zip(rows[r], rows[col], strict=True)
                    ]
    return [row[n:] for row in rows]


def real_sqrt(F: CyclotomicField, a: Value) -> Value | None:
    """The positive real square root of ``a`` in ``F``, or ``None`` when ``a`` has
    none in ``F``.

    ``a`` must be totally positive (real and positive in every embedding); a
    square root in ``F`` of such an element is real, so its embeddings are
    ``+-sqrt(sigma_k(a))`` with equal signs on ``k`` and ``-k``. Clearing
    denominators makes the root an algebraic integer, hence an integer vector
    over the power basis (``Z[zeta_M]`` is the maximal order), which is
    recovered from each sign pattern by the inverse Vandermonde matrix and
    verified exactly. No sign pattern giving an integer vector that squares to
    ``a`` is a proof that ``a`` is not a square in ``F``.
    """
    if not F.is_real(a):
        return None
    M = F.conductor
    if F.is_rational(a):
        r = a[0]
        if r < 0:
            return None
        p, q = math.isqrt(r.numerator), math.isqrt(r.denominator)
        if p * p == r.numerator and q * q == r.denominator:
            return F.rational(Fraction(p, q))
        # a rational non-square may still be a square in F (sqrt2 in Q(zeta_8))
    D = _denominator_lcm(a)
    A = F.scale(a, D * D)
    ks = units(M)
    n = len(ks)
    try:
        images = [embed(F, A, k) for k in ks]
    except OverflowError:
        return None  # beyond the floating-point reconstruction (huge denominators)
    if any(abs(z.imag) > 1e-6 or z.real <= 0 for z in images):
        return None
    roots = [math.sqrt(z.real) for z in images]
    position = {k: i for i, k in enumerate(ks)}
    pairs = [k for k in ks if k < M - k]  # one representative per {k, -k}
    Vinv = _inverse_vandermonde(M)
    for signs in product((1, -1), repeat=len(pairs) - 1):
        s = [0.0] * n
        for sign, k in zip((1,) + signs, pairs, strict=True):
            s[position[k]] = sign * roots[position[k]]
            s[position[M - k]] = s[position[k]]
        coefficients = []
        ok = True
        for row in Vinv:
            z = sum((v * t for v, t in zip(row, s, strict=True)), 0j)
            c = round(z.real)
            if abs(z.real - c) > 1e-6 or abs(z.imag) > 1e-6:
                ok = False
                break
            coefficients.append(Fraction(c, D))
        if not ok:
            continue
        x = F.reduce(coefficients)
        if F.mul(x, x) == a:
            return x
    return None


def is_abelian_over_q(K: CyclotomicField, a: Value) -> bool:
    """Whether ``K(sqrt(a))`` is abelian over ``Q`` for a totally positive ``a`` in
    ``K`` that is not a square in ``K``.

    ``K(sqrt a)/Q`` is Galois iff ``sigma(a)/a`` is a square ``b_sigma^2`` in
    ``K`` for every ``sigma`` in ``Gal(K/Q)``, and then abelian iff the lifts
    commute: ``sigma(b_rho) b_sigma = rho(b_sigma) b_rho`` for all pairs (the
    condition does not depend on the signs of the ``b``).
    """
    ks = units(K.conductor)
    b: dict[int, Value] = {}
    inverse = K.inverse(a)
    for k in ks:
        root = real_sqrt(K, K.mul(galois(K, a, k), inverse))
        if root is None:
            return False
        b[k] = root
    for k in ks:
        for m in ks:
            if k < m:
                left = K.mul(galois(K, b[m], k), b[k])
                right = K.mul(galois(K, b[k], m), b[m])
                if left != right:
                    return False
    return True


def _characters_of_order_two(H: list[int], modulus: int) -> list[dict[int, int]]:
    """The non-trivial characters ``H -> {+1, -1}`` of a finite abelian group of
    residues: coordinates of ``H / H^2`` over ``F_2`` are built greedily."""
    squares = {t * t % modulus for t in H}
    coords: dict[int, tuple[int, ...]] = {t: () for t in squares}
    group = set(squares)
    for t in H:
        if t in group:
            continue
        new = {t * e % modulus: coords[e] + (1,) for e in group}
        for e in list(group):
            coords[e] = coords[e] + (0,)
        coords.update(new)
        group |= set(new)
    rank = len(next(iter(coords.values())))
    characters = []
    for w in product((0, 1), repeat=rank):
        if any(w):
            characters.append(
                {
                    t: -1
                    if sum(x * y for x, y in zip(w, coords[t], strict=True)) % 2
                    else 1
                    for t in H
                }
            )
    return characters


def root_in_extension(
    K: CyclotomicField, a: Value, Fp: CyclotomicField
) -> Value | None:
    """``sqrt(a)`` as an element of ``Fp = Q(zeta_M')`` (``M | M'``), or ``None``.

    The quadratic extensions of ``K`` inside ``Fp`` are ``K(y_psi)`` for the
    characters ``psi`` of order two of ``H = Gal(Fp/K)``, with ``y_psi = sum_t
    psi(t) sigma_t(theta)`` for a generic real ``theta`` (``y_psi`` is real and
    ``y_psi^2`` lies in ``K``). ``sqrt(a)`` lies in ``Fp`` iff ``a`` or some
    ``a / y_psi^2`` is a square in ``K``.
    """
    M, Mp = K.conductor, Fp.conductor
    direct = real_sqrt(K, a)
    if direct is not None:
        return K.lift(direct, Fp)
    H = [t for t in units(Mp) if t % M == 1]
    lifted = K.lift(a, Fp)
    for psi in _characters_of_order_two(H, Mp):
        for m in range(1, 6):
            theta = Fp.total(
                Fp.scale(Fp.add(Fp.zeta(j), Fp.zeta(-j)), j) for j in range(1, m + 1)
            )
            y = Fp.total(Fp.scale(galois(Fp, theta, t), psi[t]) for t in H)
            if not Fp.is_zero(y):
                break
        else:
            raise ArithmeticError("no generic element for the Gauss sum")
        delta = Fp.descend(Fp.mul(y, y), K)
        if delta is None:
            raise ArithmeticError("y_psi^2 is not in K")
        c = real_sqrt(K, K.divide(a, delta))
        if c is not None:
            root = Fp.mul(y, K.lift(c, Fp))
            if Fp.mul(root, root) != lifted:
                raise ArithmeticError("the reconstructed root does not square to a")
            return root
    return None


@dataclass
class NormAnalysis:
    """What is known about one squared norm ``N`` of the field ``K``."""

    value: Value
    in_field: bool
    abelian: bool
    conductor: int | None  # smallest M' (multiple of M) with sqrt(N) in Q(zeta_M')
    root: Value | None  # sqrt(N) in Q(zeta_conductor)
    hermitian_factor: Value | None  # mu in K with mu conj(mu) = N, when found

    def describe(self, K: CyclotomicField) -> str:
        if self.in_field:
            where = f"sqrt in Q(zeta_{K.conductor})"
        elif not self.abelian:
            where = "NOT cyclotomic (K(sqrt N)/Q non-abelian)"
        elif self.conductor is None:
            where = f"cyclotomic, conductor > searched bound (phi <= {MAX_DEGREE})"
        else:
            where = f"sqrt in Q(zeta_{self.conductor}), not in Q(zeta_{K.conductor})"
        if self.hermitian_factor is None:
            norm = "no Hermitian-norm factor in K found"
        else:
            norm = f"= |mu|^2, mu = {pretty(K, self.hermitian_factor)}"
        return f"N = {pretty(K, self.value):<24} {where}; {norm}"


@cache
def analyse_norm(conductor: int, a: Value) -> NormAnalysis:
    K = cyclotomic_field(conductor)
    root = real_sqrt(K, a)
    if root is not None:
        return NormAnalysis(a, True, True, conductor, root, root)
    abelian = is_abelian_over_q(K, a)
    found: int | None = None
    found_root: Value | None = None
    if abelian:
        k = 2
        while euler_phi(k * conductor) <= MAX_DEGREE:
            Fp = cyclotomic_field(k * conductor)
            r = root_in_extension(K, a, Fp)
            if r is not None:
                found, found_root = k * conductor, r
                break
            k += 1
    return NormAnalysis(a, False, abelian, found, found_root, hermitian_factor(K, a))


def _short_sums(K: CyclotomicField, k_max: int) -> list[Value]:
    """The non-zero sums ``1 + zeta^j_2 + ... + zeta^j_k`` of ``k <= k_max`` distinct
    roots of unity; ``1 - zeta^j`` is among them for even ``M``, and odd ``M`` gets
    it explicitly."""
    M = K.conductor
    out: list[Value] = []
    if M % 2:
        out += [K.sub(K.one, K.zeta(j)) for j in range(1, M)]
    for k in range(2, k_max + 1):
        for exponents in combinations(range(1, M), k - 1):
            x = K.total([K.one] + [K.zeta(j) for j in exponents])
            if not K.is_zero(x):
                out.append(x)
    # sums with a doubled constant term, ``2 + zeta^j + ...`` (``phi^2 = 2 + zeta_5 +
    # zeta_5^-1`` is one; ``phi^2 + i phi^-2`` has norm 7), with fewer roots so that
    # the list stays short for large ``M``
    extra = k_max - 1 if M <= 24 else 2
    for k in range(1, extra + 1):
        for exponents in combinations(range(1, M), k):
            x = K.total([K.rational(2)] + [K.zeta(j) for j in exponents])
            if not K.is_zero(x):
                out.append(x)
    # when i is in K, every element is u + i v with u, v real: the norms u^2 + v^2
    # of short real sums (``phi^2 + i phi^-2`` of norm 7 in Q(zeta_20) is one)
    if M % 4 == 0:
        reals: list[Value] = []
        size = 2 if M <= 24 else 1
        for c in range(3):
            for k in range(size + 1):
                for exponents in combinations(range(1, M // 2), k):
                    reals.append(
                        K.total(
                            [K.rational(c)]
                            + [K.add(K.zeta(j), K.zeta(-j)) for j in exponents]
                        )
                    )
        i = K.zeta(M // 4)
        for a, u in enumerate(reals):
            for v in reals[a:]:
                x = K.add(u, K.mul(i, v))
                if not K.is_zero(x):
                    out.append(x)
    return out


def _is_prime(p: int) -> bool:
    if p < 2:
        return False
    for q in range(2, math.isqrt(p) + 1):
        if p % q == 0:
            return False
    return True


@cache
def _quadratic_characters(M: int, count: int) -> list[tuple[int, int]]:
    """``count`` pairs ``(p, r)`` with ``p`` a prime ``1 mod M`` and ``r`` a primitive
    ``M``-th root of unity mod ``p``: each is a prime ideal ``P`` of ``Q(zeta_M)`` of
    residue degree one (``zeta -> r``), carrying the quadratic character
    ``x -> (x mod P / p)``."""
    out: list[tuple[int, int]] = []
    p = 1000 - 1000 % M + 1  # large primes keep the candidates out of the ideals
    while len(out) < count:
        if _is_prime(p):
            for g in range(2, p):
                r = pow(g, (p - 1) // M, p)
                if all(
                    pow(r, M // q, p) != 1
                    for q in range(2, M + 1)
                    if M % q == 0 and _is_prime(q)
                ):
                    break
            else:
                raise ArithmeticError(f"no primitive root of unity mod {p}")
            for k in units(M):
                out.append((p, pow(r, k, p)))
        p += M
    return out[:count]


def _legendre(x: Value, p: int, r: int) -> int:
    """The quadratic character of ``x`` at the prime ideal ``(p, zeta - r)``: ``1``,
    ``-1``, or ``0`` when ``x`` lies in the ideal or has ``p`` in a denominator."""
    value = 0
    for j, c in enumerate(x):
        if c:
            if c.denominator % p == 0:
                return 0
            value = (value + c.numerator * pow(c.denominator, -1, p) * pow(r, j, p)) % p
    if value == 0:
        return 0
    return 1 if pow(value, (p - 1) // 2, p) == 1 else -1


@cache
@cache
def _candidate_table(
    M: int, k_max: int, rows: int
) -> tuple[list[Value], list[Value], list[tuple[int, int]], list[int]]:
    """The short sums of ``Q(zeta_M)`` of at most ``k_max`` roots of unity whose
    norms avoid the first ``rows`` prime ideals, their norms, those prime ideals,
    and the character bit vectors of the norms (bit ``i`` set when the ``i``-th
    character is ``-1``)."""
    K = cyclotomic_field(M)
    characters = _quadratic_characters(M, rows)
    candidates: list[Value] = []
    norms: list[Value] = []
    vectors: list[int] = []
    for nu in _short_sums(K, k_max):
        nn = K.norm_squared(nu)
        bits = 0
        for i, (p, r) in enumerate(characters):
            sign = _legendre(nn, p, r)
            if sign == 0:
                break
            if sign < 0:
                bits |= 1 << i
        else:
            candidates.append(nu)
            norms.append(nn)
            vectors.append(bits)
    return candidates, norms, characters, vectors


def _solve_f2(vectors: list[int], order: list[int], target: int) -> int | None:
    """A combination (bit set over ``vectors``) whose XOR is ``target``, by Gaussian
    elimination over ``F_2`` in the given processing order; ``None`` if none."""
    pivots: dict[int, tuple[int, int]] = {}  # leading bit -> (vector, combination)
    for index in order:
        v, combo = vectors[index], 1 << index
        while v:
            lead = v.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = (v, combo)
                break
            pv, pc = pivots[lead]
            v, combo = v ^ pv, combo ^ pc
    v, combo = target, 0
    while v:
        lead = v.bit_length() - 1
        if lead not in pivots:
            return None
        pv, pc = pivots[lead]
        v, combo = v ^ pv, combo ^ pc
    return combo


def hermitian_factor(K: CyclotomicField, a: Value) -> Value | None:
    """A ``mu`` in ``K`` with ``mu conj(mu) = a``, or ``None`` when none was found.

    ``mu`` is sought as ``c * prod nu_i`` with ``nu_i`` short sums of roots of unity
    and ``c`` in ``K^+``: the classes of ``a`` and of the norms ``|nu|^2`` modulo
    squares are separated by quadratic characters at prime ideals of residue
    degree one (a homomorphism killing squares), so ``a = prod |nu_i|^2 * c^2`` is
    a linear system over ``F_2``; a solution is verified exactly by ``real_sqrt``
    (``c``), the number of characters doubling while a solution fails to verify.
    Sums of two, then three, then more roots of unity are admitted in turn, and
    several processing orders are tried, to keep the product short enough for
    the floating-point step of ``real_sqrt``. No solution with the longest sums
    proves that ``a`` is outside the span of those norms; a solution that is
    too long to verify is reported as not found.
    """
    root = real_sqrt(K, a)
    if root is not None:
        return root
    M = K.conductor
    longest = 5 if M <= 24 else 4 if M <= 60 else 3
    for k_max in range(2, longest + 1):
        rows = 8 * K.degree
        while rows <= 128 * K.degree:
            candidates, _norms, characters, vectors = _candidate_table(M, k_max, rows)
            mask, target = 0, 0
            for i, (p, r) in enumerate(characters):
                sign = _legendre(a, p, r)
                if sign != 0:
                    mask |= 1 << i
                if sign < 0:
                    target |= 1 << i
            masked = [v & mask for v in vectors]
            natural = list(range(len(masked)))
            solutions = []
            for seed in range(4):
                order = list(natural)
                if seed:
                    random.Random(seed).shuffle(order)
                combo = _solve_f2(masked, order, target)
                if combo is None:
                    break  # the same span for every order: not in the span
                solutions.append(combo)
            if not solutions:
                break  # admit longer sums
            for combo in sorted(solutions, key=lambda c: c.bit_count()):
                mu = K.one
                for index, nu in enumerate(candidates):
                    if combo >> index & 1:
                        mu = K.mul(mu, nu)
                quotient = K.divide(a, K.norm_squared(mu))
                if _denominator_lcm(quotient) > 10**30:
                    continue  # too long for the floating-point reconstruction
                c = real_sqrt(K, quotient)
                if c is not None:
                    return K.mul(mu, c)
            rows *= 2  # the characters did not separate the classes: refine
    return None


# ---- Clebsch-Gordan norms ----------------------------------------------------------


@dataclass
class CGNorm:
    alpha: str
    beta: str
    gamma: str
    copy: int  # multiplicity index (1-based)
    component: tuple[int, int]  # the designated product component (j1, j2)
    norm: Value


def _matrix_at(irrep: Irrep, t: int, order: int) -> tuple[int, Matrix]:
    """``(sign, D(g))`` of the element ``t`` of ``G-dagger`` (``t >= order`` barred)."""
    sign = -1 if (t >= order and irrep.projective) else 1
    return sign, irrep.matrices[t % order]


def multiplicity(model: Model, alpha: Irrep, beta: Irrep, gamma: Irrep) -> int:
    F = model.field
    total = F.total(
        F.mul(
            F.mul(alpha.characters[g], beta.characters[g]),
            F.conjugate(gamma.characters[g]),
        )
        for g in range(model.order)
    )
    m = F.scale(total, Fraction(1, model.order))
    if not F.is_rational(m) or m[0].denominator != 1:
        raise ArithmeticError(f"non-integral multiplicity {m}")
    return int(m[0])


def _projector_entries(
    model: Model,
    alpha: Irrep,
    beta: Irrep,
    gamma: Irrep,
    rows: list[tuple[int, int]],
    j: tuple[int, int],
) -> list[Value]:
    """``P[(i1, i2), j]`` of ``P = (d_gamma / |G|) sum_g conj(D_gamma_11(g))
    D_alpha(g) x D_beta(g)`` for the given rows and the column ``j``."""
    F = model.field
    order = model.order
    scale = Fraction(gamma.dimension, order)
    j1, j2 = j
    weights = [F.conjugate(gamma.matrices[g][0][0]) for g in range(order)]
    out = []
    for i1, i2 in rows:
        total = F.total(
            F.mul(
                weights[g],
                F.mul(alpha.matrices[g][i1][j1], beta.matrices[g][i2][j2]),
            )
            for g in range(order)
            if not F.is_zero(weights[g])
        )
        out.append(F.scale(total, scale))
    return out


def cg_norms(model: Model) -> list[CGNorm]:
    """The squared norms of the Clebsch-Gordan vectors of every irrep pair."""
    F = model.field
    irreps = model.irreps
    out: list[CGNorm] = []
    for ia, alpha in enumerate(irreps):
        for beta in irreps[ia:]:
            da, db = alpha.dimension, beta.dimension
            components = [(j1, j2) for j1 in range(da) for j2 in range(db)]
            for gamma in irreps:
                if gamma.projective != (alpha.projective != beta.projective):
                    continue
                m = multiplicity(model, alpha, beta, gamma)
                if m == 0:
                    continue
                found: list[tuple[list[Value], Value]] = []  # (vector, squared norm)
                for j in components:
                    if len(found) == m:
                        break
                    diagonal = _projector_entries(model, alpha, beta, gamma, [j], j)[0]
                    if F.is_zero(diagonal):
                        continue
                    if m == 1:
                        found.append(([], diagonal))
                        out.append(
                            CGNorm(
                                alpha.symbol.text,
                                beta.symbol.text,
                                gamma.symbol.text,
                                1,
                                j,
                                diagonal,
                            )
                        )
                        break
                    w = _projector_entries(model, alpha, beta, gamma, components, j)
                    for u, nu in found:
                        overlap = F.total(
                            F.mul(F.conjugate(x), y) for x, y in zip(u, w, strict=True)
                        )
                        factor = F.divide(overlap, nu)
                        w = [
                            F.sub(y, F.mul(factor, x))
                            for x, y in zip(u, w, strict=True)
                        ]
                    if all(F.is_zero(x) for x in w):
                        continue
                    norm = w[components.index(j)]
                    if norm != F.total(F.norm_squared(x) for x in w):
                        raise ArithmeticError(
                            "the designated component is not the squared norm"
                        )
                    found.append((w, norm))
                    out.append(
                        CGNorm(
                            alpha.symbol.text,
                            beta.symbol.text,
                            gamma.symbol.text,
                            len(found),
                            j,
                            norm,
                        )
                    )
                if len(found) != m:
                    raise ArithmeticError(
                        f"found {len(found)} of {m} copies of {gamma.symbol.text} in "
                        f"{alpha.symbol.text} x {beta.symbol.text}"
                    )
    return out


# ---- chain-adapted bases -------------------------------------------------------------


def closure(table: list[list[int]], generators: list[int]) -> list[int]:
    reached = {0}
    frontier = [0]
    while frontier:
        g = frontier.pop()
        for s in generators:
            h = table[g][s]
            if h not in reached:
                reached.add(h)
                frontier.append(h)
    return sorted(reached)


def embed_by_words(
    H_table: list[list[int]],
    H_gens: list[int],
    G_table: list[list[int]],
    G_gens: list[int],
) -> list[int] | None:
    """The map ``H -> G`` sending the generators as given and extended along
    words, when it is an injective homomorphism; ``None`` otherwise."""
    image = {0: 0}
    frontier = [0]
    while frontier:
        h = frontier.pop()
        for s, t in zip(H_gens, G_gens, strict=True):
            hs = H_table[h][s]
            gt = G_table[image[h]][t]
            if hs in image:
                if image[hs] != gt:
                    return None
            else:
                image[hs] = gt
                frontier.append(hs)
    if len(image) != len(H_table):
        return None
    phi = [image[h] for h in range(len(H_table))]
    if len(set(phi)) != len(phi):
        return None
    for a in range(len(H_table)):
        for b in range(len(H_table)):
            if phi[H_table[a][b]] != G_table[phi[a]][phi[b]]:
                return None
    return phi


def label_index(model: Model, text: str) -> int:
    for k, label in enumerate(model.labels):
        if label.text == text:
            return k
    raise KeyError(text)


def element_order(table: list[list[int]], g: int) -> int:
    n, h = 1, g
    while h != 0:
        h = table[h][g]
        n += 1
    return n


def subgroup_generators(model: Model, rule: str, within: list[int]) -> list[int]:
    """Generators (element indices of ``G``) of the subgroup a chain step names,
    chosen inside ``within`` (the previous subgroup of the chain, or ``G``)."""
    table = model.table
    inverses = model.inverses
    binary = [
        g
        for g in within
        if model.group.elements[g][1] == 1 and element_order(table, g) == 2
    ]
    if rule.startswith("axis:"):
        a = label_index(model, rule[5:])
        if a not in within:
            raise ValueError(f"{rule[5:]} is not in the previous subgroup")
        return [a]
    if rule.startswith("dihedral:"):
        a = label_index(model, rule[9:])
        for c in binary:
            if table[table[c][a]][c] == inverses[a]:
                return [a, c]
        raise ValueError(f"no binary rotation inverts {rule[9:]}")
    if rule == "tetrahedral":
        a = label_index(model, "C3^1(1)")
        for c in binary:
            if len(closure(table, [a, c])) == 12:
                return [a, c]
        raise ValueError("no tetrahedral subgroup found")
    if rule == "d2":
        a = binary[0]
        for c in binary[1:]:
            if table[a][c] == table[c][a] and len(closure(table, [a, c])) == 4:
                return [a, c]
        raise ValueError("no D2 subgroup found")
    if rule == "c2":
        return [binary[0]]
    raise ValueError(rule)


def model_generators(H: Model) -> list[int]:
    """Generators of the library model of the subgroup, matching ``subgroup_generators``."""
    name = H.spec.symbol
    if name == "T":
        return [label_index(H, "C3^1(1)"), label_index(H, "C2(1)")]
    n = H.spec.family_n
    assert n is not None
    if name.startswith("D"):
        return [1, n]  # (1, 0) and (0, 1) of the dihedral rotation group
    return [1]


def embed_double(G: Model, H: Model, G_gens: list[int]) -> list[int]:
    """``H-dagger -> G-dagger`` on the lifted generators, trying the barred lifts."""
    H_gens = model_generators(H)
    for bars in product((0, 1), repeat=len(G_gens)):
        images = [g + bar * G.order for g, bar in zip(G_gens, bars, strict=True)]
        phi = embed_by_words(H.double_table, H_gens, G.double_table, images)
        if phi is not None:
            return phi
    raise ValueError(
        f"{H.spec.symbol} does not embed in {G.spec.symbol} on these generators"
    )


def isotypic_projector(
    F: CyclotomicField,
    G: Model,
    gamma: Irrep,
    lifted: list[Matrix],
    H: Model,
    sigma: Irrep,
    phi: list[int],
) -> Matrix:
    """``(d_sigma / |H-dagger|) sum_h conj(chi_sigma(h)) D_gamma(phi(h))`` over
    ``H-dagger``, in the chain field ``F`` (``lifted`` the matrices of ``gamma``
    lifted into ``F``; the characters of ``sigma`` are lifted from the field of
    ``H``)."""
    d = gamma.dimension
    out = [[F.zero] * d for _ in range(d)]
    if sigma.projective != gamma.projective:
        return out
    for t, g in enumerate(phi):
        chi = F.conjugate(H.field.lift(sigma.character(t, H.order), F))
        if F.is_zero(chi):
            continue
        sign = -1 if (g >= G.order and gamma.projective) else 1
        D = lifted[g % G.order]
        if sign < 0:
            chi = F.neg(chi)
        for i in range(d):
            for j in range(d):
                if not F.is_zero(D[i][j]):
                    out[i][j] = F.add(out[i][j], F.mul(chi, D[i][j]))
    scale = Fraction(sigma.dimension, 2 * H.order)
    return [[F.scale(x, scale) for x in row] for row in out]


def rank(F: CyclotomicField, A: Matrix) -> int:
    rows = [list(r) for r in A]
    r = 0
    cols = len(A[0]) if A else 0
    for c in range(cols):
        pivot = next(
            (i for i in range(r, len(rows)) if not F.is_zero(rows[i][c])), None
        )
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inv = F.inverse(rows[r][c])
        rows[r] = [F.mul(inv, x) for x in rows[r]]
        for i in range(len(rows)):
            if i != r and not F.is_zero(rows[i][c]):
                f = rows[i][c]
                rows[i] = [
                    F.sub(x, F.mul(f, y)) for x, y in zip(rows[i], rows[r], strict=True)
                ]
        r += 1
    return r


@dataclass
class ChainNorm:
    gamma: str
    labels: tuple[str, ...]  # the subgroup irreps along the chain
    component: int
    norm: Value


@dataclass
class ChainResult:
    chain: str
    conductor: int  # the chain field: lcm of the document conductors along the chain
    norms: list[ChainNorm] = field(default_factory=list)
    non_canonical: list[str] = field(default_factory=list)
    # the products N_i N_k whose square roots the chain-adapted matrices need
    # (entries <v_i|D(g)|v_k> / sqrt(N_i N_k) that are non-zero for some g)
    matrix_products: dict[Value, list[str]] = field(default_factory=dict)


def chain_norms(G: Model, steps: list[tuple[str, str]]) -> ChainResult:
    """The squared norms of the chain-adapted vectors ``P e_j`` of every irrep of
    ``G`` for the chain ``G > H_1 > H_2 > ...``, computed in the chain field
    ``Q(zeta_lcm(M_G, M_H1, ...))`` (the isotypic projectors of a subgroup need
    its characters, which may leave the field of ``G``)."""
    subgroups: list[tuple[Model, list[int]]] = []
    conductor = G.conductor
    within = list(range(G.order))
    for name, rule in steps:
        H = build_model(parse_group(name))
        gens = subgroup_generators(G, rule, within)
        within = closure(G.table, gens)
        if len(within) != H.order:
            raise ValueError(
                f"the generators of {name} in {G.spec.symbol} close on the wrong order"
            )
        subgroups.append((H, embed_double(G, H, gens)))
        conductor = math.lcm(conductor, H.conductor)
    F = cyclotomic_field(conductor)
    result = ChainResult(
        " > ".join([G.spec.symbol] + [name for name, _ in steps]), conductor
    )
    for gamma in G.irreps:
        d = gamma.dimension
        lifted = [
            [[G.field.lift(x, F) for x in row] for row in D] for D in gamma.matrices
        ]
        total_rank = 0
        vectors: list[tuple[list[Value], Value, tuple[str, ...]]] = []
        for chosen in product(*(H.irreps for H, _ in subgroups)):
            if any(s.projective != gamma.projective for s in chosen):
                continue
            P: Matrix | None = None
            for (H, phi), sigma in zip(subgroups, chosen, strict=True):
                Q = isotypic_projector(F, G, gamma, lifted, H, sigma, phi)
                if P is None:
                    P = Q
                else:
                    P = [
                        [
                            F.total(F.mul(P[i][k], Q[k][j]) for k in range(d))
                            for j in range(d)
                        ]
                        for i in range(d)
                    ]
            assert P is not None
            r = rank(F, P)
            total_rank += r
            if r == 0:
                continue
            labels = tuple(s.symbol.text for s in chosen)
            if r > 1:
                result.non_canonical.append(f"{gamma.symbol.text} > {labels}: rank {r}")
            j = next(i for i in range(d) if not F.is_zero(P[i][i]))
            result.norms.append(ChainNorm(gamma.symbol.text, labels, j, P[j][j]))
            vectors.append(([P[i][j] for i in range(d)], P[j][j], labels))
        if total_rank != d:
            raise ArithmeticError(
                f"the chain projectors of {gamma.symbol.text} do not sum to 1"
            )
        if len(vectors) == d:
            _matrix_products(F, gamma, lifted, vectors, result)
    return result


def _matrix_products(
    F: CyclotomicField,
    gamma: Irrep,
    lifted: list[Matrix],
    vectors: list[tuple[list[Value], Value, tuple[str, ...]]],
    result: ChainResult,
) -> None:
    """Record ``N_i N_k`` for every pair of chain-adapted vectors with a non-zero
    matrix element ``<v_i|D(g)|v_k>`` for some ``g``: the chain-adapted matrices
    with real positive normalisation are these elements over ``sqrt(N_i N_k)``."""
    d = gamma.dimension
    needed: set[tuple[int, int]] = set()
    for D in lifted:
        for i, (v, _, _) in enumerate(vectors):
            if all((i, k) in needed for k in range(d)):
                continue
            Dv = [F.total(F.mul(D[r][c], v[c]) for c in range(d)) for r in range(d)]
            for k, (w, _, _) in enumerate(vectors):
                if (i, k) in needed:
                    continue
                element = F.total(F.mul(F.conjugate(w[r]), Dv[r]) for r in range(d))
                if not F.is_zero(element):
                    needed.add((i, k))
    for i, k in sorted(needed):
        _, ni, li = vectors[i]
        _, nk, lk = vectors[k]
        result.matrix_products.setdefault(F.mul(ni, nk), []).append(
            f"{gamma.symbol.text}: <{' > '.join(li)}|D|{' > '.join(lk)}>"
        )


# ---- the survey ----------------------------------------------------------------------


def _report(K: CyclotomicField, grouped: dict[Value, list[str]], what: str) -> None:
    """One line per distinct value, then the closure summary: the smallest
    cyclotomic field holding every square root (the field of the real-positive
    phase rule) and whether every value is a Hermitian norm from ``K``."""
    closure: int | None = K.conductor
    factors = 0
    for value in sorted(grouped, key=lambda v: embed(K, v).real):
        analysis = analyse_norm(K.conductor, value)
        examples = grouped[value]
        print(f"  {analysis.describe(K)}  [{len(examples)}x, e.g. {examples[0]}]")
        if analysis.conductor is None or closure is None:
            closure = None
        else:
            closure = math.lcm(closure, analysis.conductor)
        factors += analysis.hermitian_factor is not None
    where = "no cyclotomic field" if closure is None else f"Q(zeta_{closure})"
    print(
        f"  => {what}: real-positive normalisation needs {where}; "
        f"Hermitian-norm factors in Q(zeta_{K.conductor}) found for "
        f"{factors} of {len(grouped)} values"
    )


def survey_group(name: str) -> None:
    started = time.perf_counter()
    model = build_model(parse_group(name))
    K = model.field
    print(f"\n## {name}  (|G| = {model.order}, field Q(zeta_{K.conductor}))")
    norms = cg_norms(model)
    grouped: dict[Value, list[str]] = {}
    for record in norms:
        copy = "" if record.copy == 1 else f" #{record.copy}"
        grouped.setdefault(record.norm, []).append(
            f"{record.alpha} x {record.beta} > {record.gamma}{copy}"
        )
    print(
        f"Clebsch-Gordan: {len(norms)} vectors (first partners) over "
        f"{len({(r.alpha, r.beta) for r in norms})} irrep pairs, "
        f"{len(grouped)} distinct squared norms ({time.perf_counter() - started:.1f} s)"
    )
    _report(K, grouped, "Clebsch-Gordan vectors")
    if model.group.type == "proper" and name in CHAINS:
        for steps in CHAINS[name]:
            started = time.perf_counter()
            result = chain_norms(model, steps)
            chain_grouped: dict[Value, list[str]] = {}
            for cn in result.norms:
                chain_grouped.setdefault(cn.norm, []).append(
                    f"{cn.gamma} > {' > '.join(cn.labels)}"
                )
            print(
                f"Chain {result.chain} (chain field Q(zeta_{result.conductor})): "
                f"{len(result.norms)} chain-adapted vectors, "
                f"{len(chain_grouped)} distinct squared norms "
                f"({time.perf_counter() - started:.1f} s)"
            )
            for line in result.non_canonical:
                print(f"  non-canonical: {line}")
            Fc = cyclotomic_field(result.conductor)
            _report(Fc, chain_grouped, "chain-adapted vectors")
            print(
                f"  products N_i N_k under the non-zero matrix elements "
                f"<v_i|D(g)|v_k>: {len(result.matrix_products)} distinct"
            )
            _report(Fc, result.matrix_products, "chain-adapted matrices")


def main(argv: list[str]) -> int:
    groups = argv[1:] or DEFAULT_GROUPS
    print("# Cyclotomic closure survey (docs/research/scripts/cyclotomic_closure.py)")
    for name in groups:
        survey_group(name)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
