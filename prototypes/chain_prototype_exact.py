"""Prototype: exact arithmetic for chain-adapted bases (wayfinder ticket #50).

Throwaway code for the "Exact computation" section of docs/spec/group-chains.md:

* ``real_sqrt`` -- the positive real square root of a totally positive element
  of ``Q(zeta_M)`` by the modular decision procedure of the chapter (sign search
  over residues modulo a prime ``p = 1 (mod m)``, exact verification), in the
  smallest cyclotomic subfield holding the element and then along the quadratic
  sub-extensions up to ``Q(zeta_M)``;
* ``hermitian_factor`` -- a search for ``mu`` with ``mu * conj(mu) = N``
  (products of short sums of roots of unity, quadratic characters over ``F_2``),
  ported from the research script of ticket #40 with the floating-point root
  replaced by ``real_sqrt``; it generates the shipped norm-factor entries;
* ``norm_factor`` -- the canonical norm factor ``mu_F(N)``: the real root when
  it exists, ``(1 + i)/2`` for ``N = 1/2`` without ``sqrt2``, else a table entry.

Nothing here is production code. Floating point enters in one place only: the
bound ``B`` on the coefficients of a root (with a factor-two margin); the
decision itself is exact.
"""

from __future__ import annotations

import cmath
import math
import random
from fractions import Fraction
from functools import cache
from itertools import combinations, product

from pointgroup_tables.field import CyclotomicField, Value, cyclotomic_field

TWO_PI = 2 * math.pi


# ---- helpers --------------------------------------------------------------------------


def units(M: int) -> list[int]:
    return [k for k in range(1, M + 1) if math.gcd(k, M) == 1]


def galois(F: CyclotomicField, x: Value, k: int) -> Value:
    """The Galois conjugate ``sigma_k: zeta -> zeta^k`` of ``x`` (``gcd(k, M) = 1``)."""
    return F.total(F.scale(F.zeta(j * k), c) for j, c in enumerate(x) if c)


def embed(F: CyclotomicField, x: Value, k: int = 1) -> complex:
    """The complex embedding ``zeta -> exp(2 pi i k / M)`` (bounds and signs only)."""
    M = F.conductor
    return sum(
        (float(c) * cmath.exp(1j * TWO_PI * j * k / M) for j, c in enumerate(x) if c),
        0j,
    )


def denominator_lcm(x: Value) -> int:
    out = 1
    for c in x:
        out = math.lcm(out, c.denominator)
    return out


def smallest_subfield(F: CyclotomicField, x: Value) -> tuple[CyclotomicField, Value]:
    """``x`` in the smallest ``Q(zeta_m)``, ``m | M``, that holds it."""
    M = F.conductor
    for m in sorted(d for d in range(1, M + 1) if M % d == 0):
        if m % 2 == 0 and (m // 2) % 2 == 1:
            continue  # Q(zeta_{2k}) = Q(zeta_k) for odd k: skip the duplicate
        K = cyclotomic_field(m)
        y = F.descend(x, K)
        if y is not None:
            return K, y
    raise AssertionError("unreachable: x lies in F itself")


# ---- primes and modular roots ---------------------------------------------------------


def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_one_mod(m: int, above: int) -> int:
    """The first prime ``p = 1 (mod m)`` with ``p > above``."""
    p = above - above % m + 1
    if p <= above:
        p += m
    while not is_probable_prime(p):
        p += m
    return p


def primitive_root_of_unity(m: int, p: int) -> int:
    """A primitive ``m``-th root of unity modulo the prime ``p = 1 (mod m)``."""
    prime_factors = [q for q in range(2, m + 1) if m % q == 0 and is_probable_prime(q)]
    for g in range(2, p):
        r = pow(g, (p - 1) // m, p)
        if all(pow(r, m // q, p) != 1 for q in prime_factors):
            return r
    raise ArithmeticError("no primitive root of unity")


def sqrt_mod(a: int, p: int) -> int | None:
    """A square root of ``a`` modulo the odd prime ``p`` (Tonelli-Shanks), or ``None``."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def _solve_mod(matrix: list[list[int]], rhs: list[int], p: int) -> list[int]:
    """Gaussian elimination modulo ``p`` (the matrix is invertible)."""
    n = len(matrix)
    rows = [[x % p for x in row] + [b % p] for row, b in zip(matrix, rhs, strict=True)]
    for col in range(n):
        pivot = next(r for r in range(col, n) if rows[r][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        inv = pow(rows[col][col], -1, p)
        rows[col] = [x * inv % p for x in rows[col]]
        for r in range(n):
            if r != col and rows[r][col]:
                f = rows[r][col]
                rows[r] = [
                    (x - f * y) % p for x, y in zip(rows[r], rows[col], strict=True)
                ]
    return [row[n] for row in rows]


@cache
def _inverse_vandermonde_bound(m: int) -> float:
    """``max_j sum_k |Vinv[j][k]|`` for ``V[k][j] = exp(2 pi i k j / m)``, ``k`` a unit:
    the coefficient of a real element is bounded by this times the largest
    embedding modulus (floating point, used for the bound only)."""
    ks = units(m)
    n = len(ks)
    rows = [
        [cmath.exp(1j * TWO_PI * k * j / m) for j in range(n)]
        + [1.0 + 0j if r == c else 0j for c in range(n)]
        for r, k in enumerate(ks)
    ]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(rows[r][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        piv = rows[col][col]
        rows[col] = [v / piv for v in rows[col]]
        for r in range(n):
            if r != col:
                f = rows[r][col]
                if f:
                    rows[r] = [
                        a - f * b for a, b in zip(rows[r], rows[col], strict=True)
                    ]
    return max(sum(abs(v) for v in row[n:]) for row in rows)


# ---- the exact real square root ---------------------------------------------------------


def _sqrt_in_own_field(K: CyclotomicField, a: Value) -> Value | None:
    """``sqrt(a)`` in ``K = Q(zeta_m)`` for a real ``a`` of ``K`` (``m`` not 2 mod 4),
    or ``None`` when ``a`` is not a square in ``K``."""
    m = K.conductor
    if K.is_rational(a):
        r = a[0]
        if r < 0:
            return None
        pn, qn = math.isqrt(r.numerator), math.isqrt(r.denominator)
        if pn * pn == r.numerator and qn * qn == r.denominator:
            return K.rational(Fraction(pn, qn))
        if m <= 2:
            return None
    D = denominator_lcm(a)
    A = K.scale(a, D * D)
    ints = [int(c) for c in A]
    ks = units(m)
    n = len(ks)
    if n == 1:
        return None  # rational, handled above
    # |sigma_k(x)|^2 = sigma_k(A) <= sum |a_j| (|zeta| = 1), so every coefficient of x is
    # bounded by the inverse-Vandermonde row sum times sqrt(sum |a_j|): an integer bound
    # that never overflows (A may have hundreds of digits after a norm-factor search)
    size = sum(abs(c) for c in ints)
    B = 2 * math.ceil(_inverse_vandermonde_bound(m)) * (math.isqrt(size) + 1) + 2
    attempt = 0
    p = prime_one_mod(m, 2 * B)
    while True:
        r = primitive_root_of_unity(m, p)
        sigma = [sum(c * pow(r, k * j, p) for j, c in enumerate(ints)) % p for k in ks]
        if all(sigma):
            break
        attempt += 1
        p = prime_one_mod(m, p)
        if attempt > 20:
            raise ArithmeticError("no prime avoiding the conjugates")
    roots: list[int] = []
    for s in sigma:
        root = sqrt_mod(s, p)
        if root is None:
            return None  # sigma_k(A) is a non-residue: A is not a square in K
        roots.append(root)
    position = {k: i for i, k in enumerate(ks)}
    pairs = [k for k in ks if k < m - k]
    V = [[pow(r, k * j, p) for j in range(n)] for k in ks]
    for signs in product((1, -1), repeat=len(pairs) - 1):
        s = [0] * n
        for sign, k in zip((1,) + signs, pairs, strict=True):
            value = roots[position[k]] if sign > 0 else (p - roots[position[k]]) % p
            s[position[k]] = value
            s[position[m - k]] = value
        c = _solve_mod(V, s, p)
        lifted = [Fraction(x if x <= p // 2 else x - p, D) for x in c]
        x = K.reduce(lifted)
        if K.mul(x, x) == a:
            if embed(K, x).real < 0:  # x has moderate coefficients: |x|^2 = a
                x = K.neg(x)
            return x
    return (
        None  # no sign pattern verifies: a is not a square in K (nor totally positive)
    )


def _quadratic_characters(H: list[int], modulus: int) -> list[dict[int, int]]:
    """The non-trivial characters ``H -> {+1, -1}`` of a finite abelian group of
    residues (coordinates of ``H / H^2`` over ``F_2``, built greedily)."""
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


@cache
def real_sqrt(F: CyclotomicField, a: Value) -> Value | None:
    """The positive real square root of ``a`` in ``F``, or ``None`` when ``a`` is not
    a square in ``F`` (exact; see the module docstring)."""
    if not F.is_real(a):
        return None
    K, b = smallest_subfield(F, a)
    direct = _sqrt_in_own_field(K, b)
    if direct is not None:
        return K.lift(direct, F)
    if K.conductor == F.conductor:
        return None
    # along the quadratic sub-extensions K(y_psi) of F / K: sqrt(a) in F iff
    # a / y_psi^2 is a square in K for some order-two character psi of Gal(F/K)
    M, m = F.conductor, K.conductor
    H = [t for t in units(M) if t % m == 1 % m]
    for psi in _quadratic_characters(H, M):
        if (M - 1) in psi and psi[M - 1] == -1:
            continue  # a real root lies in a real sub-extension: psi(-1) = 1
        y: Value | None = None
        for size in range(1, 8):
            theta = F.total(
                F.scale(F.add(F.zeta(j), F.zeta(-j)), j) for j in range(1, size + 1)
            )
            candidate = F.total(F.scale(galois(F, theta, t), psi[t]) for t in H)
            if not F.is_zero(candidate):
                y = candidate
                break
        if y is None:
            raise ArithmeticError("no generic element for the Gauss sum")
        delta = F.descend(F.mul(y, y), K)
        if delta is None:
            raise ArithmeticError("y_psi^2 is not in K")
        c = _sqrt_in_own_field(K, K.divide(b, delta))
        if c is not None:
            root = F.mul(y, K.lift(c, F))
            if F.mul(root, root) != a:
                raise ArithmeticError("the reconstructed root does not square to a")
            if embed(F, root).real < 0:
                root = F.neg(root)
            return root
    return None


# ---- Hermitian-norm factors (the search behind the shipped entries) ---------------------


def _short_sums(K: CyclotomicField, k_max: int) -> list[Value]:
    M = K.conductor
    out: list[Value] = []
    if M % 2:
        out += [K.sub(K.one, K.zeta(j)) for j in range(1, M)]
    for k in range(2, k_max + 1):
        for exponents in combinations(range(1, M), k - 1):
            x = K.total([K.one] + [K.zeta(j) for j in exponents])
            if not K.is_zero(x):
                out.append(x)
    extra = k_max - 1 if M <= 24 else 2
    for k in range(1, extra + 1):
        for exponents in combinations(range(1, M), k):
            x = K.total([K.rational(2)] + [K.zeta(j) for j in exponents])
            if not K.is_zero(x):
                out.append(x)
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


@cache
def _character_primes(M: int, count: int) -> list[tuple[int, int]]:
    """``count`` pairs ``(p, r)``: a prime ``p = 1 (mod M)`` and a primitive ``M``-th
    root of unity modulo ``p``, one prime ideal of residue degree one each."""
    out: list[tuple[int, int]] = []
    p = 1000 - 1000 % M + 1
    while len(out) < count:
        if is_probable_prime(p):
            r = primitive_root_of_unity(M, p)
            for k in units(M):
                out.append((p, pow(r, k, p)))
        p += M
    return out[:count]


def _legendre(x: Value, p: int, r: int) -> int:
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
def _candidate_table(
    M: int, k_max: int, rows: int
) -> tuple[list[Value], list[tuple[int, int]], list[int]]:
    K = cyclotomic_field(M)
    characters = _character_primes(M, rows)
    candidates: list[Value] = []
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
            vectors.append(bits)
    return candidates, characters, vectors


def _solve_f2(vectors: list[int], order: list[int], target: int) -> int | None:
    pivots: dict[int, tuple[int, int]] = {}
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


def hermitian_factor(K: CyclotomicField, a: Value, seeds: int = 6) -> Value | None:
    """A ``mu`` in ``K`` with ``mu * conj(mu) = a``, or ``None`` when the search found
    none (not a proof). ``mu = c * prod nu_i``, ``nu_i`` short sums of roots of
    unity, ``c`` real, found as an ``F_2`` system over quadratic characters and
    verified by ``real_sqrt``."""
    root = real_sqrt(K, a)
    if root is not None:
        return root
    M = K.conductor
    longest = 5 if M <= 24 else 4 if M <= 60 else 3
    for k_max in range(2, longest + 1):
        rows = 8 * K.degree
        while rows <= 128 * K.degree:
            candidates, characters, vectors = _candidate_table(M, k_max, rows)
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
            for seed in range(seeds):
                order = list(natural)
                if seed:
                    random.Random(seed).shuffle(order)
                combo = _solve_f2(masked, order, target)
                if combo is None:
                    break
                solutions.append(combo)
            if not solutions:
                break
            for combo in sorted(set(solutions), key=lambda c: c.bit_count()):
                mu = K.one
                for index, nu in enumerate(candidates):
                    if combo >> index & 1:
                        mu = K.mul(mu, nu)
                quotient = K.divide(a, K.norm_squared(mu))
                c = real_sqrt(K, quotient)
                if c is not None:
                    return K.mul(mu, c)
            rows *= 2
    return None


# ---- the canonical norm factor ----------------------------------------------------------

NormFactorTable = dict[tuple[int, Value], Value]


def norm_factor(F: CyclotomicField, N: Value, table: NormFactorTable) -> Value | None:
    """``mu_F(N)`` of the chapter: the real positive root in ``F``; ``(1 + i)/2`` for
    ``N = 1/2`` without ``sqrt2``; otherwise the table entry ``(M, N)``, or ``None``."""
    root = real_sqrt(F, N)
    if root is not None:
        return root
    M = F.conductor
    if N == F.rational(Fraction(1, 2)) and M % 4 == 0:
        i = F.zeta(M // 4)
        return F.scale(F.add(F.one, i), Fraction(1, 2))
    return table.get((M, N))
