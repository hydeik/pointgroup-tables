"""The version-1 canonical basis versus the chain-adapted basis through the
principal element.

Computational evidence for ``docs/research/canonical-versus-chain.md``
(wayfinder research ticket #49). Run with ``uv run python3
docs/research/scripts/canonical_versus_chain.py [GROUP>H1>H2 ...]``; without
arguments the whole survey of the note is run (the axial ``sweep`` set and the
polyhedral chains through the principal element). The script uses only the
standard library and the exact arithmetic of ``pointgroup_tables.field``.

What it computes
----------------

For every irrep of a group ``G`` and a chain ``G > H_1 > ... > H_m`` it builds
the chain-adapted basis exactly as ``docs/spec/group-chains.md`` ("Chain-adapted
basis" and "Exact computation") defines it:

- the step subgroups are resolved by the embedding rule of the spec (admissible
  isomorphisms of double groups preserving parity and label invariants,
  filtered by canonicity, the lexicographically smallest image chosen);
- the label spaces are the ranges of products of exact isotypic projectors,
  checked to have rank ``0`` or ``1`` and to sum to the dimension;
- the first vector of every block is ``P_t e_j / mu(N)`` with ``j`` the
  designated canonical component and ``mu(N)`` the canonical norm factor
  (the exact real square root when it lies in the chain field, ``(1 + i)/2``
  for ``N = 1/2`` otherwise);
- the other vectors of a block are transferred from the first one by the
  matrix-element projector built from the sub-chain realisation of the block's
  irrep, computed recursively by the same definition in the library model of
  ``H_1``.

It then compares the chain-adapted basis ``(w_s)`` with the canonical basis
``(e_k)`` of the version-1 realisation: the component permutation, the
diagonal unitary where every ``w_s`` is a phase times a canonical vector, the
exact mixing block where it is not, whether the chain-adapted matrices equal the
canonical ones entry by entry, up to the permutation, or need the phases too,
and the squared norm ``N`` of every designated component.
"""

from __future__ import annotations

import cmath
import math
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cache
from itertools import product

from pointgroup_tables.document import Model, build_model
from pointgroup_tables.field import CyclotomicField, Matrix, Value, cyclotomic_field
from pointgroup_tables.irreps import Irrep
from pointgroup_tables.symbols import parse_group
from pointgroup_tables.values import pretty

TWO_PI = 2 * math.pi

# the survey: every chain through the principal element of the version-2 map
AXIAL_FAMILIES = [
    ("Dn", "D{n}", "C{n}"),
    ("Cnv", "C{n}v", "C{n}"),
    ("Dnh", "D{n}h", "C{n}h"),
    ("Dnd", "D{n}d", "S{m}"),
]
SWEEP = range(2, 13)  # n = 1 resolves to the aliases C2, C1h, C2v, C2h (see the note)
POLYHEDRAL = [
    "T>C3",
    "Th>S6",
    "Td>S4",
    "O>C4",
    "Oh>C4h",
    "O>D4>C4",
    "Oh>D4h>C4h",
    "Td>D2d>S4",
    "I>D5>C5",
    "Ih>D5d>S10",
]


# ---- exact square roots and the norm factor -----------------------------------------


def embed(F: CyclotomicField, x: Value, k: int = 1) -> complex:
    """The complex embedding ``zeta -> exp(2 pi i k / M)`` of ``x``."""
    M = F.conductor
    return sum(
        (float(c) * cmath.exp(1j * TWO_PI * j * k / M) for j, c in enumerate(x) if c),
        0j,
    )


def units(M: int) -> list[int]:
    return [k for k in range(1, M + 1) if math.gcd(k, M) == 1]


def _denominator_lcm(x: Value) -> int:
    out = 1
    for c in x:
        out = math.lcm(out, c.denominator)
    return out


@cache
def _inverse_vandermonde(M: int) -> list[list[complex]]:
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
                rows[r] = [v - f * w for v, w in zip(rows[r], rows[col], strict=True)]
    return [row[n:] for row in rows]


def real_sqrt(F: CyclotomicField, a: Value) -> Value | None:
    """The positive real square root of ``a`` in ``F`` or ``None`` (the
    floating-point reconstruction of the cyclotomic-closure script with exact
    verification; it suffices for the small denominators met here)."""
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
    D = _denominator_lcm(a)
    A = F.scale(a, D * D)
    ks = units(M)
    n = len(ks)
    images = [embed(F, A, k) for k in ks]
    if any(abs(z.imag) > 1e-6 or z.real <= 0 for z in images):
        return None
    roots = [math.sqrt(z.real) for z in images]
    position = {k: i for i, k in enumerate(ks)}
    pairs = [k for k in ks if k < M - k]
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


@dataclass
class NormFactor:
    N: Value
    mu: Value
    how: str  # "square", "(1+i)/2"


def norm_factor(F: CyclotomicField, N: Value) -> NormFactor:
    """The canonical norm factor ``mu_F(N)`` of the spec for the cases a chain
    through the principal element produces; anything else is reported."""
    root = real_sqrt(F, N)
    if root is not None:
        return NormFactor(N, root, "square")
    if F.is_rational(N) and F.as_rational(N) == Fraction(1, 2) and F.conductor % 4 == 0:
        i = F.zeta(F.conductor // 4)
        return NormFactor(N, F.scale(F.add(F.one, i), Fraction(1, 2)), "(1+i)/2")
    raise ArithmeticError(
        f"no canonical norm factor known for N = {pretty(F, N)} in Q(zeta_{F.conductor})"
    )


# ---- group helpers --------------------------------------------------------------------


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


def element_order(table: list[list[int]], g: int) -> int:
    n, h = 1, g
    while h != 0:
        h = table[h][g]
        n += 1
    return n


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


def invariants(model: Model, g: int) -> tuple[str, int, int | None, int]:
    """The invariants of the canonical label of element ``g`` that an admissible
    isomorphism preserves: kind, rotation order, modulus of the exponent, parity."""
    label = model.labels[g]
    power = None if label.power is None else abs(label.power)
    return (label.kind, label.order, power, model.group.elements[g][1])


def generating_set(model: Model) -> list[int]:
    """A small generating set of ``H`` in canonical element order: the first
    element outside the closure so far, repeated."""
    gens: list[int] = []
    reached = [0]
    for g in range(1, model.order):
        if g not in reached:
            gens.append(g)
            reached = closure(model.table, gens)
            if len(reached) == model.order:
                break
    return gens


@dataclass
class Step:
    H: Model
    phi: list[int]  # H-dagger index -> G-dagger index (the embedding into G)


def restriction_multiplicities(
    K: Model, H: Model, psi: list[int]
) -> list[tuple[str, str, int]]:
    """``(tau, sigma, m)`` for every irrep ``tau`` of ``K`` and ``sigma`` of ``H``
    of the same kind, ``psi: H-dagger -> K-dagger`` the embedding."""
    L = cyclotomic_field(math.lcm(K.conductor, H.conductor))
    out = []
    for tau in K.irreps:
        for sigma in H.irreps:
            if tau.projective != sigma.projective:
                continue
            total = L.zero
            for h, k in enumerate(psi):
                a = K.field.lift(tau.character(k, K.order), L)
                b = L.conjugate(H.field.lift(sigma.character(h, H.order), L))
                total = L.add(total, L.mul(a, b))
            m = L.scale(total, Fraction(1, 2 * H.order))
            assert L.is_rational(m) and L.as_rational(m).denominator == 1, (
                tau.symbol.text,
                sigma.symbol.text,
            )
            out.append((tau.symbol.text, sigma.symbol.text, int(L.as_rational(m))))
    return out


def resolve_step(G: Model, K: Step | None, H: Model) -> Step:
    """The embedding rule of ``docs/spec/group-chains.md``: among the admissible
    isomorphisms ``H-dagger -> K-dagger`` (``K = G`` for the first step) keep the
    canonical ones and take the lexicographically smallest image."""
    within = (
        set(range(G.order)) if K is None else {g % G.order for g in K.phi[: K.H.order]}
    )
    gens = generating_set(H)
    candidates: list[list[int]] = []
    for h in gens:
        inv = invariants(H, h)
        half_turn = element_order(H.table, h) == 2
        images = []
        for g in sorted(within):
            if invariants(G, g) != inv:
                continue
            images.append(g)
            if half_turn:
                images.append(g + G.order)
        candidates.append(images)
    best: list[int] | None = None
    for choice in product(*candidates):
        phi = embed_by_words(H.double_table, gens, G.double_table, list(choice))
        if phi is None:
            continue
        ok = True
        for h in range(H.order):
            g = phi[h]
            if g % G.order not in within or invariants(G, g % G.order) != invariants(
                H, h
            ):
                ok = False
                break
            if element_order(H.table, h) != 2 and g >= G.order:
                ok = False  # the lifted preimage must go to the lifted preimage
                break
        if not ok:
            continue
        if best is not None and phi >= best:
            continue
        # canonicity for both kinds: restriction of K_lib to H through psi = K.phi^-1 . phi
        if K is None:
            psi = phi
        else:
            back = {g: h for h, g in enumerate(K.phi)}
            if any(g not in back for g in phi):
                continue
            psi = [back[g] for g in phi]
        K_model = G if K is None else K.H
        if any(m > 1 for _, _, m in restriction_multiplicities(K_model, H, psi)):
            continue
        best = phi
    if best is None:
        raise ValueError(
            f"no canonical admissible embedding of {H.spec.symbol} in {G.spec.symbol}"
        )
    return Step(H, best)


def resolve_chain(G: Model, names: list[str]) -> list[Step]:
    steps: list[Step] = []
    for name in names:
        H = build_model(parse_group(name))
        steps.append(resolve_step(G, steps[-1] if steps else None, H))
    return steps


# ---- the chain-adapted basis --------------------------------------------------------


def matrix_mul(F: CyclotomicField, A: Matrix, B: Matrix) -> Matrix:
    n, m, p = len(A), len(B), len(B[0])
    return [
        [F.total(F.mul(A[i][k], B[k][j]) for k in range(m)) for j in range(p)]
        for i in range(n)
    ]


def apply(F: CyclotomicField, A: Matrix, v: list[Value]) -> list[Value]:
    return [F.total(F.mul(A[i][k], v[k]) for k in range(len(v))) for i in range(len(A))]


def inner(F: CyclotomicField, u: list[Value], v: list[Value]) -> Value:
    """``<u|v> = sum conj(u_i) v_i``."""
    return F.total(F.mul(F.conjugate(a), b) for a, b in zip(u, v, strict=True))


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


def lifted_matrices(F: CyclotomicField, G: Model, gamma: Irrep) -> list[Matrix]:
    """``D_gamma`` on ``G-dagger`` (both preimages), lifted into ``F``."""
    base = [[[G.field.lift(x, F) for x in row] for row in D] for D in gamma.matrices]
    if gamma.projective:
        barred = [[[F.neg(x) for x in row] for row in D] for D in base]
    else:
        barred = base
    return base + barred


def isotypic_projector(
    F: CyclotomicField, D: list[Matrix], H: Model, sigma: Irrep, phi: list[int]
) -> Matrix:
    """``(d_sigma / |H-dagger|) sum_h conj(chi_sigma(h)) D(phi(h))`` over ``H-dagger``."""
    d = len(D[0])
    out = [[F.zero] * d for _ in range(d)]
    for h, g in enumerate(phi):
        chi = F.conjugate(H.field.lift(sigma.character(h, H.order), F))
        if F.is_zero(chi):
            continue
        M = D[g]
        for i in range(d):
            for j in range(d):
                if not F.is_zero(M[i][j]):
                    out[i][j] = F.add(out[i][j], F.mul(chi, M[i][j]))
    scale = Fraction(sigma.dimension, 2 * H.order)
    return [[F.scale(x, scale) for x in row] for row in out]


def transfer_projector(
    F: CyclotomicField,
    D: list[Matrix],
    H: Model,
    sigma_matrices: list[Matrix],
    dim_sigma: int,
    phi: list[int],
    s: int,
    t: int,
) -> Matrix:
    """``(d_sigma / |H-dagger|) sum_h conj(D'_sigma(h)_{s t}) D(phi(h))``."""
    d = len(D[0])
    out = [[F.zero] * d for _ in range(d)]
    for h, g in enumerate(phi):
        c = F.conjugate(sigma_matrices[h][s][t])
        if F.is_zero(c):
            continue
        M = D[g]
        for i in range(d):
            for j in range(d):
                if not F.is_zero(M[i][j]):
                    out[i][j] = F.add(out[i][j], F.mul(c, M[i][j]))
    scale = Fraction(dim_sigma, 2 * H.order)
    return [[F.scale(x, scale) for x in row] for row in out]


@dataclass
class Adapted:
    """The chain-adapted basis of one irrep."""

    gamma: Irrep
    labels: list[tuple[str, ...]]  # chain labels in component order
    vectors: list[list[Value]]  # w_s in canonical coordinates
    matrices: list[Matrix]  # D'_gamma on G-dagger (both preimages), chain field
    designated: list[
        tuple[tuple[str, ...], int, NormFactor]
    ]  # (first label of a block, j, mu)


def chain_basis(F: CyclotomicField, G: Model, steps: list[Step]) -> dict[str, Adapted]:
    """The chain-adapted basis of every irrep of ``G`` for the resolved chain,
    by the recursive definition of the spec, computed in the chain field ``F``."""
    out: dict[str, Adapted] = {}
    if not steps:
        for gamma in G.irreps:
            D = lifted_matrices(F, G, gamma)
            d = gamma.dimension
            e = [[F.one if i == k else F.zero for i in range(d)] for k in range(d)]
            out[gamma.symbol.text] = Adapted(gamma, [() for _ in range(d)], e, D, [])
        return out
    first = steps[0]
    # the pulled-back chain H_1 > phi_1^-1(H_2) > ...
    back = {g: h for h, g in enumerate(first.phi)}
    pulled = [Step(s.H, [back[g] for g in s.phi]) for s in steps[1:]]
    sub = chain_basis(F, first.H, pulled)
    for gamma in G.irreps:
        D = lifted_matrices(F, G, gamma)
        d = gamma.dimension
        # label spaces
        spaces: list[tuple[tuple[str, ...], Matrix]] = []
        total = 0
        kinds = [
            [s for s in step.H.irreps if s.projective == gamma.projective]
            for step in steps
        ]
        projectors = [
            {s.symbol.text: isotypic_projector(F, D, step.H, s, step.phi) for s in kind}
            for step, kind in zip(steps, kinds, strict=True)
        ]
        for chosen in product(*kinds):
            P = None
            for step_projectors, s in zip(projectors, chosen, strict=True):
                Q = step_projectors[s.symbol.text]
                P = Q if P is None else matrix_mul(F, P, Q)
            assert P is not None
            r = rank(F, P)
            if r > 1:
                raise ArithmeticError(
                    f"{G.spec.symbol} {gamma.symbol.text}: label space of rank {r}"
                )
            total += r
            if r == 1:
                spaces.append((tuple(s.symbol.text for s in chosen), P))
        if total != d:
            raise ArithmeticError(
                f"{G.spec.symbol} {gamma.symbol.text}: ranks sum to {total}, not {d}"
            )
        labels = [s for s, _ in spaces]
        P_of = dict(spaces)
        vectors: dict[tuple[str, ...], list[Value]] = {}
        designated = []
        # blocks by the first label
        blocks: dict[str, list[tuple[str, ...]]] = {}
        for s in labels:
            blocks.setdefault(s[0], []).append(s)
        for sigma_text, members in blocks.items():
            t = members[0]
            P_t = P_of[t]
            j = next(i for i in range(d) if not F.is_zero(P_t[i][i]))
            mu = norm_factor(F, P_t[j][j])
            inv_mu = F.inverse(mu.mu)
            w_t = [F.mul(P_t[i][j], inv_mu) for i in range(d)]
            assert inner(F, w_t, w_t) == F.one
            vectors[t] = w_t
            designated.append((t, j, mu))
            realisation = sub[sigma_text]
            sub_index = {s: k for k, s in enumerate(realisation.labels)}
            assert [s[1:] for s in members] == realisation.labels, (
                members,
                realisation.labels,
            )
            for s in members[1:]:
                T = transfer_projector(
                    F,
                    D,
                    first.H,
                    realisation.matrices,
                    realisation.gamma.dimension,
                    first.phi,
                    sub_index[s[1:]],
                    sub_index[t[1:]],
                )
                w_s = apply(F, T, w_t)
                assert inner(F, w_s, w_s) == F.one, (gamma.symbol.text, s)
                assert apply(F, P_of[s], w_s) == w_s, (gamma.symbol.text, s)
                vectors[s] = w_s
        W = [vectors[s] for s in labels]
        adapted = [
            [[inner(F, W[a], apply(F, M, W[b])) for b in range(d)] for a in range(d)]
            for M in D
        ]
        out[gamma.symbol.text] = Adapted(gamma, labels, W, adapted, designated)
    # within-block check: the step-1 restriction equals the sub-chain realisation
    for a in out.values():
        for sigma_text, realisation in sub.items():
            rows = [k for k, s in enumerate(a.labels) if s[0] == sigma_text]
            if not rows:
                continue
            for h, g in enumerate(first.phi):
                block = [[a.matrices[g][r][c] for c in rows] for r in rows]
                assert block == realisation.matrices[h], (
                    G.spec.symbol,
                    a.gamma.symbol.text,
                    sigma_text,
                )
    return out


# ---- comparison with the canonical basis ---------------------------------------------


@dataclass
class Comparison:
    gamma: str
    labels: list[tuple[str, ...]]
    permutation: list[int | None]  # pi(s): the canonical index, None if mixed
    phases: list[str]  # lambda_s pretty, or "mixed"
    mixing: list[tuple[tuple[str, ...], list[str]]]  # non-monomial vectors, pretty
    monomial: bool
    identity_order: bool
    identity_phases: bool
    matrices: str  # "equal", "permuted", "phased", "mixed"
    designated: list[str]


def compare(F: CyclotomicField, G: Model, a: Adapted) -> Comparison:
    d = a.gamma.dimension
    D = lifted_matrices(F, G, a.gamma)
    permutation: list[int | None] = []
    phases: list[str] = []
    mixing = []
    for s, w in zip(a.labels, a.vectors, strict=True):
        support = [k for k in range(d) if not F.is_zero(w[k])]
        if len(support) == 1:
            permutation.append(support[0])
            phases.append(pretty(F, w[support[0]]))
        else:
            permutation.append(None)
            phases.append("mixed")
            mixing.append((s, [pretty(F, x) for x in w]))
    monomial = all(p is not None for p in permutation)
    identity_order = permutation == list(range(d))
    identity_phases = monomial and all(
        w[k] == F.one
        for w, k in zip(a.vectors, permutation, strict=True)
        if k is not None
    )
    if not monomial:
        relation = "mixed"
    elif a.matrices == D:
        relation = "equal"
    elif all(
        a.matrices[g][r][c] == D[g][permutation[r]][permutation[c]]  # type: ignore[index]
        for g in range(len(D))
        for r in range(d)
        for c in range(d)
    ):
        relation = "permuted"
    else:
        relation = "phased"
    designated = [
        f"{'>'.join(t)}: j={j} N={pretty(F, mu.N)} mu={pretty(F, mu.mu)} ({mu.how})"
        for t, j, mu in a.designated
    ]
    return Comparison(
        a.gamma.symbol.text,
        a.labels,
        permutation,
        phases,
        mixing,
        monomial,
        identity_order,
        identity_phases,
        relation,
        designated,
    )


@dataclass
class ChainReport:
    chain: str
    conductor: int  # the chain field L
    comparisons: list[Comparison] = field(default_factory=list)
    embeddings: list[str] = field(default_factory=list)

    @property
    def summary(self) -> str:
        rel = {c.gamma: c.matrices for c in self.comparisons}
        counts = {
            k: sum(1 for v in rel.values() if v == k)
            for k in ("equal", "permuted", "phased", "mixed")
        }
        norms = sorted(
            {
                x.split("N=")[1].split(" mu=")[0]
                for c in self.comparisons
                for x in c.designated
            }
        )
        return (
            f"{self.chain}: L={self.conductor}; matrices equal {counts['equal']}, permuted {counts['permuted']}, "
            f"phased {counts['phased']}, mixed {counts['mixed']}; N in {{{', '.join(norms)}}}"
        )


def survey_chain(spelling: str) -> ChainReport:
    names = spelling.split(">")
    G = build_model(parse_group(names[0]))
    steps = resolve_chain(G, names[1:])
    L = G.conductor
    for step in steps:
        L = math.lcm(L, step.H.conductor)
    F = cyclotomic_field(L)
    report = ChainReport(
        " > ".join([G.spec.symbol] + [s.H.spec.symbol for s in steps]), L
    )
    for step in steps:
        gens = generating_set(step.H)
        images = ", ".join(
            f"{step.H.labels[h].text} -> {G.double_labels[step.phi[h]].text}"
            for h in gens
        )
        report.embeddings.append(f"{step.H.spec.symbol}: {images}")
    for adapted in chain_basis(F, G, steps).values():
        report.comparisons.append(compare(F, G, adapted))
    return report


def print_report(report: ChainReport, verbose: bool) -> None:
    print(f"== {report.summary}")
    for e in report.embeddings:
        print(f"   embedding {e}")
    for c in report.comparisons:
        if not verbose and c.matrices == "equal" and c.identity_phases:
            continue
        pi = ", ".join("-" if p is None else str(p) for p in c.permutation)
        print(
            f"   {c.gamma}: labels {['>'.join(s) for s in c.labels]}; pi = ({pi})"
            f"{' = id' if c.identity_order else ''}; phases {c.phases}; matrices {c.matrices}"
        )
        for t in c.designated:
            print(f"      designated {t}")
        for s, w in c.mixing:
            print(f"      w_{'>'.join(s)} = ({', '.join(w)})")


def main(argv: list[str]) -> int:
    verbose = "-v" in argv
    args = [a for a in argv if a != "-v"]
    if args:
        chains = args
    else:
        chains = []
        for _, G, H in AXIAL_FAMILIES:
            for n in SWEEP:
                chains.append(f"{G.format(n=n)}>{H.format(n=n, m=2 * n)}")
        chains += POLYHEDRAL
    summaries = []
    for spelling in chains:
        report = survey_chain(spelling)
        print_report(report, verbose)
        summaries.append(report.summary)
    print()
    print("Summary")
    for s in summaries:
        print(s)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
