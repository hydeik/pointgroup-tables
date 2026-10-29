"""Prototype: chain-adapted bases by exact projectors (wayfinder ticket #50).

Throwaway code for docs/spec/group-chains.md, "Chain-adapted basis" and "Exact
computation", on top of the library:

* the subgroup catalogue and the chain classes of a group (ported from the
  research script of ticket #41, ``docs/research/scripts/subgroup_chains.py``
  on ``research/subgroup-chains``), the admissible isomorphisms of the
  embedding rule and the "lexicographically smallest image" choice;
* the label spaces as products of exact isotypic projectors, the designated
  vector ``P_t e_j / mu_F(N)``, the within-block transfer by the matrix-element
  projector of the sub-chain realisation (recursive), the chain-adapted
  matrices, and the checks of the "Properties" section;
* the chain conductor by the rule of the chapter and the norm-factor entries a
  chain needs (``hermitian_factor`` search for the missing ones);
* the ``subgroup_chain`` table and the ``components`` of ``matrices``.

Nothing here is production code.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cache
from itertools import combinations
from typing import Any

from chain_prototype_exact import (
    NormFactorTable,
    hermitian_factor,
    norm_factor,
    real_sqrt,
)

from pointgroup_tables.document import Model, build_model
from pointgroup_tables.field import (
    CyclotomicField,
    Matrix,
    Value,
    cyclotomic_field,
    field_metadata,
)
from pointgroup_tables.irreps import Irrep, named_elements
from pointgroup_tables.labels import REFLECTION, ROTATION, reduced_label
from pointgroup_tables.symbols import parse_group

TERMINAL_FAMILIES = ("Cn", "Cnh", "S2n")
Subgroup = frozenset[int]


class ChainError(Exception):
    """A chain the prototype cannot realise (the internal error of the chapter)."""


# ---- subgroups (ported from the #41 research script) -----------------------------------


def closure(table: list[list[int]], generators: list[int]) -> Subgroup:
    members = {0}
    frontier = [0]
    while frontier:
        following = []
        for x in frontier:
            for g in generators:
                y = table[x][g]
                if y not in members:
                    members.add(y)
                    following.append(y)
        frontier = following
    return frozenset(members)


def all_subgroups(table: list[list[int]]) -> dict[Subgroup, list[int]]:
    cyclic: dict[Subgroup, int] = {}
    for g in range(len(table)):
        cyclic.setdefault(closure(table, [g]), g)
    subgroups: dict[Subgroup, list[int]] = {S: [g] for S, g in cyclic.items()}
    frontier = list(subgroups)
    while frontier:
        following = []
        for S in frontier:
            for C, c in cyclic.items():
                if C <= S:
                    continue
                joined = closure(table, subgroups[S] + [c])
                if joined not in subgroups:
                    subgroups[joined] = subgroups[S] + [c]
                    following.append(joined)
        frontier = following
    return subgroups


def conjugate(
    table: list[list[int]], inverses: list[int], S: Subgroup, g: int
) -> Subgroup:
    return frozenset(table[table[g][x]][inverses[g]] for x in S)


def sort_key(S: Subgroup) -> tuple[int, ...]:
    return tuple(sorted(S))


def rotation_kind(model: Model, rotations: set[Any]) -> tuple[str, int]:
    R = model.group.rotation_group
    m = len(rotations)
    orders = [R.order_of(r) for r in rotations]
    if max(orders) == m:
        return ("C", m)
    if m % 2 == 0 and max(orders) == m // 2:
        return ("D", m // 2)
    return {12: ("T", 12), 24: ("O", 24), 60: ("I", 60)}[m]


def identify(model: Model, S: Subgroup) -> str:
    """The normalised Schoenflies symbol of a subgroup."""
    G = model.group
    elements = [G.elements[i] for i in S]
    proper = {r for r, p in elements if p == 1}
    improper = {r for r, p in elements if p == -1}
    kind, n = rotation_kind(model, proper)
    if not improper:
        symbol = {"C": f"C{n}", "D": f"D{n}", "T": "T", "O": "O", "I": "I"}[kind]
    elif (G.rotation_group.identity, -1) in elements:
        if kind == "C":
            symbol = f"C{n}h" if n % 2 == 0 else f"S{2 * n}"
        elif kind == "D":
            symbol = f"D{n}h" if n % 2 == 0 else f"D{n}d"
        else:
            symbol = kind + "h"
    else:
        full_kind, full_n = rotation_kind(model, proper | improper)
        if full_kind == "C":
            symbol = f"C{n}h" if n % 2 else f"S{2 * n}"
        elif full_kind == "D" and kind == "C":
            symbol = f"C{full_n}v"
        elif full_kind == "D":
            symbol = f"D{n}h" if n % 2 else f"D{n}d"
        else:
            symbol = "Td"
    return parse_group(symbol).symbol


def invariant(model: Model, index: int) -> tuple[Any, ...]:
    """Parity, label kind, rotation order, |exponent| and, unless a half turn,
    whether the element of ``G-dagger`` is barred."""
    N = model.order
    base, barred = index % N, index >= N
    label = model.labels[base]
    parity = model.group.elements[base][1]
    half_turn = label.kind in (ROTATION, REFLECTION) and label.order == 2
    power = abs(label.power) if label.power is not None else 0
    return (parity, label.kind, label.order, power, None if half_turn else barred)


def element_orders(table: list[list[int]]) -> list[int]:
    orders = []
    for g in range(len(table)):
        k, x = 1, g
        while x != 0:
            x = table[x][g]
            k += 1
        orders.append(k)
    return orders


def generating_set(table: list[list[int]]) -> list[int]:
    orders = element_orders(table)
    generators: list[int] = []
    generated = closure(table, [0])
    for g in sorted(range(1, len(table)), key=lambda g: (-orders[g], g)):
        if g in generated:
            continue
        generators.append(g)
        generated = closure(table, generators)
        if len(generated) == len(table):
            break
    return generators


def extend(
    lib_table: list[list[int]],
    table: list[list[int]],
    generators: list[int],
    images: list[int],
) -> dict[int, int] | None:
    phi = {0: 0}
    frontier = [0]
    while frontier:
        following = []
        for x in frontier:
            for g, y in zip(generators, images, strict=True):
                xg, image = lib_table[x][g], table[phi[x]][y]
                if xg in phi:
                    if phi[xg] != image:
                        return None
                else:
                    phi[xg] = image
                    following.append(xg)
        frontier = following
    if len(phi) != len(lib_table) or len(set(phi.values())) != len(phi):
        return None
    return phi


# The designated canonical component of a block-first vector (ticket #58):
# "invariant": among the components with (P_t)_jj != 0, the first that is an eigenvector of
# the principal element of G with eigenvalue 1, else the first (the chapter's rule);
# "first": the rule of ticket #43 (the first component with (P_t)_jj != 0), kept for comparison.
DESIGNATED_RULE = "invariant"

# The embedding order rule (ticket #58):
# "principal": the images of the elements of H_lib-dagger compared in descending element order,
# ties by canonical index (the chapter's rule); "sequence": the whole canonical element list in
# its own order (the rule of ticket #42), kept for comparison.
ORDER_RULE = "principal"


def image_key(lib: Model, phi: list[int]) -> list[int]:
    """The sequence the order rule compares."""
    if ORDER_RULE == "sequence":
        return phi
    orders = element_orders(lib.double_table)
    return [phi[x] for x in sorted(range(len(phi)), key=lambda x: (-orders[x], x))]


def admissible_isomorphisms(model: Model, lib: Model, S: Subgroup) -> list[list[int]]:
    """Every admissible isomorphism ``H_lib-dagger -> preimage of S in G-dagger``
    (the embedding rule), each as the image list over the canonical element list
    of ``H_lib-dagger``."""
    N = model.order
    double = sorted(S) + [i + N for i in sorted(S)]
    targets: dict[tuple[Any, ...], list[int]] = {}
    for i in double:
        targets.setdefault(invariant(model, i), []).append(i)
    lib_table, table = lib.double_table, model.double_table
    generators = generating_set(lib_table)
    candidates = [targets.get(invariant(lib, g), []) for g in generators]
    lib_invariants = [invariant(lib, x) for x in range(len(lib_table))]
    found: list[list[int]] = []

    def search(k: int, images: list[int]) -> None:
        if k == len(generators):
            phi = extend(lib_table, table, generators, images)
            if phi is None:
                return
            if all(lib_invariants[x] == invariant(model, i) for x, i in phi.items()):
                found.append([phi[x] for x in range(len(lib_table))])
            return
        for y in candidates[k]:
            if y not in images:
                search(k + 1, images + [y])

    search(0, [])
    return found


def branching(model: Model, lib: Model, phi: list[int]) -> list[dict[int, int]]:
    """Per irrep of ``G`` the multiplicities of the irreps of ``H_lib`` in its restriction."""
    F = cyclotomic_field(math.lcm(model.conductor, lib.conductor))
    FG, FH = model.field, lib.field
    order = 2 * lib.order
    lifted_H = [
        [F.conjugate(FH.lift(gamma.character(x, lib.order), F)) for x in range(order)]
        for gamma in lib.irreps
    ]
    out = []
    for irrep in model.irreps:
        chi = [FG.lift(irrep.character(phi[x], model.order), F) for x in range(order)]
        row: dict[int, int] = {}
        total = 0
        for k, gamma in enumerate(lib.irreps):
            value = F.scale(
                F.total(F.mul(a, b) for a, b in zip(chi, lifted_H[k], strict=True)),
                Fraction(1, order),
            )
            if not F.is_rational(value) or value[0].denominator != 1:
                raise RuntimeError(f"non-integral multiplicity {value}")
            m = int(value[0])
            if m:
                row[k] = m
                total += m * gamma.dimension
        if total != irrep.dimension:
            raise RuntimeError("the restriction does not exhaust the irrep")
        out.append(row)
    return out


@dataclass
class SubgroupClass:
    symbol: str
    members: list[Subgroup]
    phi: list[int]
    multiplicities: list[dict[int, int]]

    @property
    def representative(self) -> Subgroup:
        return self.members[0]

    @property
    def order(self) -> int:
        return len(self.members[0])

    def max_multiplicity(self, model: Model, projective: bool) -> int:
        return max(
            (
                max(row.values())
                for irrep, row in zip(model.irreps, self.multiplicities, strict=True)
                if irrep.projective == projective
            ),
            default=1,
        )


@cache
def catalogue(symbol: str) -> list[SubgroupClass]:
    """The proper subgroups of the group up to conjugacy, each identified and
    embedded (the representative is the conjugate with the smallest image)."""
    model = build_model(parse_group(symbol))
    table, inverses = model.table, model.inverses
    N = model.order
    seen: set[Subgroup] = set()
    classes: list[SubgroupClass] = []
    for S in sorted(all_subgroups(table), key=lambda S: (-len(S), sort_key(S))):
        if S in seen or len(S) == N:
            continue
        orbit = {conjugate(table, inverses, S, g) for g in range(N)}
        seen |= orbit
        members = sorted(orbit, key=sort_key)
        lib = build_model(parse_group(identify(model, members[0])))
        best: tuple[list[int], Subgroup] | None = None
        for member in members:
            for phi in admissible_isomorphisms(model, lib, member):
                if best is None or image_key(lib, phi) < image_key(lib, best[0]):
                    best = (phi, member)
        assert best is not None
        phi, representative = best
        members.remove(representative)
        members.insert(0, representative)
        classes.append(
            SubgroupClass(lib.spec.symbol, members, phi, branching(model, lib, phi))
        )
    return classes


def is_terminal(symbol: str) -> bool:
    return parse_group(symbol).family in TERMINAL_FAMILIES


@dataclass
class ChainClass:
    symbols: list[str]  # G first
    subgroups: list[Subgroup]  # as subsets of G
    steps: list[tuple[int, int]]  # per step: max multiplicity (ordinary, projective)

    @property
    def both_canonical(self) -> bool:
        return all(step == (1, 1) for step in self.steps)


@cache
def chain_classes(symbol: str) -> list[ChainClass]:
    """Every chain class of the group ending at the cyclic shape through steps
    multiplicity-free for at least one kind, up to simultaneous conjugation."""
    model = build_model(parse_group(symbol))
    N = model.order
    table, inverses = model.table, model.inverses
    found: list[ChainClass] = []

    def walk(top: str, into_G: dict[int, int], chain: ChainClass) -> None:
        if is_terminal(top) and chain.symbols[1:]:
            found.append(chain)
            return
        sub_model = build_model(parse_group(top))
        for cls in catalogue(top):
            if cls.order == 1:
                continue
            step = (
                cls.max_multiplicity(sub_model, False),
                cls.max_multiplicity(sub_model, True),
            )
            if 1 in step:
                mapping = {
                    x: into_G[cls.phi[x] % sub_model.order] for x in range(cls.order)
                }
                S = frozenset(mapping.values())
                walk(
                    cls.symbol,
                    mapping,
                    ChainClass(
                        chain.symbols + [cls.symbol],
                        chain.subgroups + [S],
                        chain.steps + [step],
                    ),
                )

    walk(symbol, {i: i for i in range(N)}, ChainClass([symbol], [], []))
    distinct: dict[tuple[tuple[int, ...], ...], ChainClass] = {}
    for chain in found:
        key = min(
            tuple(sort_key(conjugate(table, inverses, S, g)) for S in chain.subgroups)
            for g in range(N)
        )
        distinct.setdefault(key, chain)
    return list(distinct.values())


def accepted_chains(symbol: str) -> list[ChainClass]:
    """The chain classes ``--chain`` accepts: canonical for both kinds (the
    terminal is of the cyclic shape by construction), sorted by number of steps
    and then by spelling."""
    out = [c for c in chain_classes(symbol) if c.both_canonical]
    out.sort(key=lambda c: (len(c.symbols), ">".join(c.symbols[1:])))
    return out


def class_reps_of(model: Model, S: Subgroup) -> list[str]:
    """The reduced class labels of ``G`` whose classes meet ``S``, in class order."""
    out = []
    for cls in model.classes:
        if any(i in S for i in cls.members):
            label = model.labels[cls.representative]
            out.append(reduced_label(label, cls.inverse != cls.index))
    return out


def signature(model: Model, S: Subgroup) -> tuple[str, ...]:
    """The class-label signature of a subgroup: the reduced class labels of ``G``
    whose classes it meets (``E`` excluded), in class order."""
    return tuple(r for r in class_reps_of(model, S) if r != "E")


def spelling(model: Model, chain: ChainClass, siblings: list[ChainClass]) -> str:
    """The canonical spelling of the chapter ("Fully qualified spelling", ticket #59):
    at step ``i`` an anchor exactly when another accepted class with the same symbols
    and the same signatures before step ``i`` has a different signature at step ``i``;
    the anchor is the shortest list of reps (canonical class order, then
    lexicographic in that order) that excludes every such class; when no list does
    (the signature is contained in every competitor's) the step stays bare and the
    order rule reaches it."""
    sigs = [signature(model, S) for S in chain.subgroups]
    parts = []
    for i, symbol in enumerate(chain.symbols[1:]):
        others = {
            signature(model, s.subgroups[i])
            for s in siblings
            if s.symbols[: i + 2] == chain.symbols[: i + 2]
            and [signature(model, S) for S in s.subgroups[:i]] == sigs[:i]
            and signature(model, s.subgroups[i]) != sigs[i]
        }
        anchor = ""
        if others:
            found: tuple[str, ...] | None = None
            for size in range(1, len(sigs[i]) + 1):
                for combo in combinations(sigs[i], size):
                    if all(any(r not in o for r in combo) for o in others):
                        found = combo
                        break
                if found:
                    break
            if found:
                anchor = "[" + ",".join(found) + "]"
        parts.append(symbol + anchor)
    return ">".join(parts)


def parse_spelling(text: str) -> list[tuple[str, list[str]]]:
    """``"D4>C2[C2']"`` -> ``[("D4", []), ("C2", ["C2'"])]``."""
    steps = []
    for part in text.split(">"):
        m = re.fullmatch(r"([^\[]+)(?:\[(.*)\])?", part)
        if m is None:
            raise ChainError(f"malformed chain {text!r}")
        steps.append((m.group(1), m.group(2).split(",") if m.group(2) else []))
    return steps


def same_chain(model: Model, c: ChainClass, chosen: list[Subgroup]) -> bool:
    """Whether ``chosen`` (subgroups of ``G``, step by step) lies in the chain class ``c``."""
    return any(
        all(
            conjugate(model.table, model.inverses, c.subgroups[j], g) == chosen[j]
            for j in range(len(c.subgroups))
        )
        for g in range(model.order)
    )


@dataclass
class TwinSet:
    """Accepted chain classes sharing one canonical spelling (twin chain classes,
    ticket #59): the embedding rule reaches exactly one of them."""

    spelling: str
    classes: list[ChainClass]
    reached: int  # index into ``classes``


@cache
def twin_sets(symbol: str) -> list[TwinSet]:
    model = build_model(parse_group(symbol))
    chains = accepted_chains(symbol)
    by: dict[str, list[ChainClass]] = {}
    for c in chains:
        by.setdefault(spelling(model, c, chains), []).append(c)
    out = []
    for text, cs in by.items():
        if len(cs) > 1:
            resolved = resolve_text(symbol, text)
            chosen = [frozenset(st.members) for st in resolved.steps]
            reached = next(k for k, c in enumerate(cs) if same_chain(model, c, chosen))
            out.append(TwinSet(text, cs, reached))
    return out


@cache
def reachable_chains(symbol: str) -> list[tuple[str, ChainClass]]:
    """What ``list chains`` holds: one entry per canonical spelling, with the chain
    class the embedding rule resolves the spelling to; sorted by number of steps and
    then by spelling in ASCII order."""
    model = build_model(parse_group(symbol))
    chains = accepted_chains(symbol)
    unreachable = {
        id(c)
        for t in twin_sets(symbol)
        for k, c in enumerate(t.classes)
        if k != t.reached
    }
    out = [(spelling(model, c, chains), c) for c in chains if id(c) not in unreachable]
    out.sort(key=lambda tc: (len(tc[1].symbols), tc[0]))
    return out


def differing_step(model: Model, a: ChainClass, b: ChainClass) -> int:
    """The first step (1-based) at which ``b`` differs from ``a`` up to simultaneous
    conjugation in ``G`` (``len + 1`` when they are the same class)."""
    best = 0
    for g in range(model.order):
        conj = [conjugate(model.table, model.inverses, S, g) for S in b.subgroups]
        d = 0
        while d < len(conj) and conj[d] == a.subgroups[d]:
            d += 1
        best = max(best, d)
    return best + 1


# ---- resolved chains ------------------------------------------------------------------


@dataclass
class Step:
    lib: Model  # H_lib
    phi: list[int]  # H_lib-dagger index -> G-dagger index
    members: Subgroup  # G indices of H


@dataclass
class ResolvedChain:
    G: Model
    steps: list[Step]
    symbols: list[str]  # G first
    text: str  # the spelling used

    @property
    def chain_field(self) -> int:
        return (
            math.lcm(self.G.conductor, *(s.lib.conductor for s in self.steps))
            if self.steps
            else self.G.conductor
        )


def resolve(symbol: str, chain: ChainClass, text: str) -> ResolvedChain:
    """The resolved chain of a chain class: step by step the admissible
    isomorphism with the lexicographically smallest image, over the subgroups
    of the class compatible with the steps already chosen."""
    model = build_model(parse_group(symbol))
    N = model.order
    table, inverses = model.table, model.inverses
    chosen: list[Subgroup] = []
    steps: list[Step] = []
    for i, name in enumerate(chain.symbols[1:]):
        lib = build_model(parse_group(name))
        candidates: set[Subgroup] = set()
        for g in range(N):
            if all(
                conjugate(table, inverses, chain.subgroups[j], g) == chosen[j]
                for j in range(i)
            ):
                candidates.add(conjugate(table, inverses, chain.subgroups[i], g))
        best: tuple[list[int], Subgroup] | None = None
        for S in candidates:
            for phi in admissible_isomorphisms(model, lib, S):
                if best is None or image_key(lib, phi) < image_key(lib, best[0]):
                    best = (phi, S)
        assert best is not None
        phi, S = best
        chosen.append(S)
        steps.append(Step(lib, phi, S))
    return ResolvedChain(model, steps, list(chain.symbols), text)


def resolve_text(symbol: str, text: str) -> ResolvedChain:
    """Resolve a chain given by a spelling (``"D4>C4"``, ``"D4>D2[C2']>C2"``,
    ``"default"``) as the chapter does: over every accepted class whose symbols
    match and whose signatures meet the anchors step by step, the lexicographically
    smallest image sequence (twin chain classes share a spelling; one is reached)."""
    model = build_model(parse_group(symbol))
    chains = accepted_chains(symbol)
    if text == "default":
        text = default_chain(symbol)
    steps = parse_spelling(text)
    candidates = [
        resolve(symbol, c, text)
        for c in chains
        if c.symbols[1:] == [s for s, _ in steps]
        and all(
            all(r in signature(model, c.subgroups[i]) for r in reps)
            for i, (_, reps) in enumerate(steps)
        )
    ]
    if candidates:
        return min(
            candidates, key=lambda rc: [image_key(s.lib, s.phi) for s in rc.steps]
        )
    raise ChainError(f"no accepted chain '{text}' for {symbol}")


DEFAULT_POLYHEDRAL = {
    "T": "C3",
    "Th": "S6",
    "O": "C4",
    "Oh": "C4h",
    "Td": "S4",
    "I": "D5>C5",
    "Ih": "D5d>S10",
}


def default_chain(symbol: str) -> str:
    spec = parse_group(symbol)
    if spec.polyhedral:
        return DEFAULT_POLYHEDRAL[spec.symbol]
    n = spec.family_n
    return {
        "Cn": "",
        "Cnh": "",
        "S2n": "",
        "Dn": f"C{n}",
        "Cnv": f"C{n}",
        "Dnh": f"C{n}h",
        "Dnd": f"S{2 * n}",
    }[spec.family]


def pulled_back(chain: ResolvedChain) -> ResolvedChain:
    """The chain ``H_1,lib > phi_1^-1(H_2) > ...`` inside the first step's library group."""
    first = chain.steps[0]
    inverse = {g: h for h, g in enumerate(first.phi)}
    steps = []
    for step in chain.steps[1:]:
        phi = [inverse[g] for g in step.phi]
        members = frozenset(h % first.lib.order for h in phi)
        steps.append(Step(step.lib, phi, members))
    return ResolvedChain(
        first.lib, steps, chain.symbols[1:], ">".join(chain.symbols[2:])
    )


# ---- projectors --------------------------------------------------------------------


def matrix_of(irrep: Irrep, t: int, order: int) -> tuple[int, Matrix]:
    """``(sign, D(g))`` at the ``G-dagger`` index ``t``."""
    sign = -1 if (t >= order and irrep.projective) else 1
    return sign, irrep.matrices[t % order]


def lift_matrix(F0: CyclotomicField, F: CyclotomicField, D: Matrix) -> Matrix:
    return [[F0.lift(x, F) for x in row] for row in D]


def mat_mul(F: CyclotomicField, A: Matrix, B: Matrix) -> Matrix:
    n, m, k = len(A), len(B[0]), len(B)
    return [
        [
            F.total(F.mul(A[i][t], B[t][j]) for t in range(k) if not F.is_zero(A[i][t]))
            for j in range(m)
        ]
        for i in range(n)
    ]


def adjoint(F: CyclotomicField, A: Matrix) -> Matrix:
    return [[F.conjugate(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def isotypic_projector(
    F: CyclotomicField,
    G: Model,
    gamma: Irrep,
    lifted: list[Matrix],
    step: Step,
    sigma: Irrep,
) -> Matrix:
    """``(d_sigma / |H-dagger|) sum_h conj(chi_sigma(h)) D_gamma(phi(h))``."""
    d = gamma.dimension
    out = [[F.zero] * d for _ in range(d)]
    if sigma.projective != gamma.projective:
        return out
    H = step.lib
    for t, g in enumerate(step.phi):
        chi = F.conjugate(H.field.lift(sigma.character(t, H.order), F))
        if F.is_zero(chi):
            continue
        sign, D = matrix_of(gamma, g, G.order)
        if sign < 0:
            chi = F.neg(chi)
        D = lifted[g % G.order]
        for i in range(d):
            for j in range(d):
                if not F.is_zero(D[i][j]):
                    out[i][j] = F.add(out[i][j], F.mul(chi, D[i][j]))
    scale = Fraction(sigma.dimension, 2 * H.order)
    return [[F.scale(x, scale) for x in row] for row in out]


def rank(F: CyclotomicField, A: Matrix) -> int:
    rows = [list(r) for r in A]
    r = 0
    for c in range(len(A[0]) if A else 0):
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


# ---- the chain-adapted basis ------------------------------------------------------------


@dataclass
class ChainIrrep:
    irrep: Irrep
    labels: list[
        tuple[int, ...]
    ]  # chain labels in component order (irrep indices per step)
    norms: list[
        Value | None
    ]  # the designated squared norm per component (None: transferred)
    designated: list[
        int | None
    ]  # the designated canonical index per block-first vector
    W: Matrix  # columns: basis vectors in canonical components, in F_doc
    matrices: list[Matrix]  # D'(g) on the elements of G (canonical section), in F_doc

    @property
    def dimension(self) -> int:
        return len(self.labels)


@dataclass
class ChainBasis:
    chain: ResolvedChain
    conductor: int  # M_chain
    irreps: list[ChainIrrep]
    field: CyclotomicField = field(init=False)

    def __post_init__(self) -> None:
        self.field = cyclotomic_field(self.conductor)


class MissingNormFactor(ChainError):
    def __init__(self, N: Value, M: int, where: str) -> None:
        super().__init__(f"no canonical norm factor for N in Q(zeta_{M}) at {where}")
        self.N, self.M, self.where = N, M, where


@cache
def step_order(symbol: str) -> list[int]:
    """The irrep indices of a step group in its step order (#57): exponent order at
    the group's own principal element for the cyclic shape, canonical order else."""
    lib = build_model(parse_group(symbol))
    if lib.spec.family not in TERMINAL_FAMILIES:
        return list(range(len(lib.irreps)))
    F = lib.field
    c = named_elements(lib.group).principal
    N = element_orders(lib.table)[c]
    keys = []
    for k, sigma in enumerate(lib.irreps):
        if sigma.projective:
            candidates = [(e, F.root_of_unity(2 * N, e)) for e in range(-N + 1, N + 1)]
        else:
            candidates = [
                (e, F.root_of_unity(N, e))
                for e in range(-(N // 2) + (1 if N % 2 == 0 else 0), N // 2 + 1)
            ]
        chi = sigma.character(c, lib.order)
        e = next(e for e, value in candidates if value == chi)
        keys.append((-e, k))
    return [k for _, k in sorted(keys)]


def _label_spaces(
    F: CyclotomicField, chain: ResolvedChain, gamma: Irrep, lifted: list[Matrix]
) -> list[tuple[tuple[int, ...], Matrix]]:
    """``(label, P_label)`` for every chain label of ``gamma`` (rank one), in
    lexicographic component order; checks ranks ``<= 1`` summing to ``d``."""
    G = chain.G
    d = gamma.dimension
    out: list[tuple[tuple[int, ...], Matrix]] = []
    total = 0

    def descend(i: int, P: Matrix | None, label: tuple[int, ...]) -> None:
        nonlocal total
        if i == len(chain.steps):
            assert P is not None
            r = rank(F, P)
            total += r
            if r > 1:
                raise ChainError(
                    f"label space of rank {r}: {gamma.symbol.text} {label}"
                )
            if r == 1:
                out.append((label, P))
            return
        step = chain.steps[i]
        for k in step_order(step.lib.spec.symbol):
            sigma = step.lib.irreps[k]
            if sigma.projective != gamma.projective:
                continue
            Q = isotypic_projector(F, G, gamma, lifted, step, sigma)
            if all(F.is_zero(x) for row in Q for x in row):
                continue
            PQ = Q if P is None else mat_mul(F, P, Q)
            if all(F.is_zero(x) for row in PQ for x in row):
                continue
            descend(i + 1, PQ, label + (k,))

    if not chain.steps:
        for k in range(d):
            P = [[F.one if i == j == k else F.zero for j in range(d)] for i in range(d)]
            out.append(((), P))
        return out
    descend(0, None, ())
    if total != d:
        raise ChainError(f"ranks sum to {total}, not {d}: {gamma.symbol.text}")
    return out


def designated_index(
    F: CyclotomicField, G: Model, lifted: list[Matrix], P: Matrix
) -> int:
    """The designated canonical index of a block-first vector with projector ``P``."""
    d = len(P)
    candidates = [i for i in range(d) if not F.is_zero(P[i][i])]
    if DESIGNATED_RULE == "invariant":
        Dc = lifted[named_elements(G.group).principal]
        fixed = [i for i in candidates if F.is_zero(F.add(Dc[i][i], F.neg(F.one)))]
        if fixed:
            return fixed[0]
    return candidates[0]


def chain_norms(
    chain: ResolvedChain,
) -> tuple[int, list[tuple[str, tuple[int, ...], Value]]]:
    """``(L, [(irrep, label, N)])``: the designated squared norms of every block-first
    vector of every irrep, in the chain field ``Q(zeta_L)``."""
    G = chain.G
    L = chain.chain_field
    F = cyclotomic_field(L)
    out = []
    for gamma in G.irreps:
        lifted = [lift_matrix(G.field, F, D) for D in gamma.matrices]
        spaces = _label_spaces(F, chain, gamma, lifted)
        seen_blocks: set[int] = set()
        for label, P in spaces:
            block = label[0] if label else -1
            if block in seen_blocks:
                continue
            seen_blocks.add(block)
            j = designated_index(F, G, lifted, P)
            out.append((gamma.symbol.text, label, P[j][j]))
    return L, out


def sub_chain_conductors(chain: ResolvedChain, table: NormFactorTable) -> list[int]:
    if not chain.steps or len(chain.steps) == 1:
        return []
    sub = pulled_back(chain)
    return [chain_conductor(sub, table)[0]] + sub_chain_conductors(sub, table)


_CONDUCTOR_CACHE: dict[
    tuple[Any, ...], tuple[int, list[tuple[str, tuple[int, ...], Value]]]
] = {}


def chain_key(chain: ResolvedChain) -> tuple[Any, ...]:
    return (chain.G.spec.symbol,) + tuple(
        (s.lib.spec.symbol, tuple(s.phi)) for s in chain.steps
    )


def chain_conductor(
    chain: ResolvedChain, table: NormFactorTable
) -> tuple[int, list[tuple[str, tuple[int, ...], Value]]]:
    """The chain conductor by the rule of the chapter: the smallest multiple of
    ``L`` dividing ``lcm(L, 24)``, a multiple of every sub-chain conductor, in
    which every designated squared norm has a canonical norm factor (a square,
    ``1/2``, a table entry, or an entry the ``hermitian_factor`` search adds to
    the table). Returns the conductor and the norms with their labels."""
    key = chain_key(chain)
    if key in _CONDUCTOR_CACHE:
        return _CONDUCTOR_CACHE[key]
    L, norms = chain_norms(chain)
    FL = cyclotomic_field(L)
    base = math.lcm(L, *sub_chain_conductors(chain, table))
    top = math.lcm(L, 24)
    candidates = sorted(m for m in range(base, top + 1, base) if top % m == 0)
    distinct = {N for _, _, N in norms}
    for M in candidates:
        F = cyclotomic_field(M)
        ok = True
        for N in sorted(distinct):
            NM = FL.lift(N, F)
            if norm_factor(F, NM, table) is not None:
                continue
            mu = hermitian_factor(F, NM)
            if mu is None:
                ok = False
                break
            table[(M, NM)] = mu
        if ok:
            _CONDUCTOR_CACHE[key] = (M, norms)
            return M, norms
    raise ChainError(f"no conductor up to {top} for {chain.G.spec.symbol} {chain.text}")


_BASIS_CACHE: dict[tuple[Any, ...], ChainBasis] = {}


def label_space_agreement(a: ResolvedChain, b: ResolvedChain) -> dict[str, list[str]]:
    """Per kind, the irreps of ``G`` whose rank-one label spaces coincide as sets
    of lines under the two chains (``same_ordinary``, ``same_projective``) and
    those whose do not (``different_*``); twin chain classes agree on every
    ordinary irrep and differ on every projective one (ticket #59)."""
    G = a.G
    out: dict[str, list[str]] = {
        "same_ordinary": [],
        "different_ordinary": [],
        "same_projective": [],
        "different_projective": [],
    }
    for gamma in G.irreps:
        lines = []
        for chain in (a, b):
            F = cyclotomic_field(chain.chain_field)
            lifted = [lift_matrix(G.field, F, D) for D in gamma.matrices]
            spaces = _label_spaces(F, chain, gamma, lifted)
            lines.append(frozenset(tuple(tuple(row) for row in P) for _, P in spaces))
        kind = "projective" if gamma.projective else "ordinary"
        key = ("same_" if lines[0] == lines[1] else "different_") + kind
        out[key].append(gamma.symbol.text)
    return out


def chain_basis(
    chain: ResolvedChain, table: NormFactorTable, M_doc: int | None = None
) -> ChainBasis:
    """The chain-adapted basis of every irrep of ``chain.G`` by the definition of
    the chapter, computed in ``Q(zeta_{M_doc})`` (``M_chain`` by default)."""
    if M_doc is None:
        M_doc = chain_conductor(chain, table)[0]
    key = chain_key(chain) + (M_doc,)
    if key in _BASIS_CACHE:
        return _BASIS_CACHE[key]
    G = chain.G
    F = cyclotomic_field(M_doc)
    sub: ChainBasis | None = None
    if chain.steps:
        sub = chain_basis(pulled_back(chain), table)  # in its own chain conductor
        if M_doc % sub.conductor:
            raise ChainError(
                f"the sub-chain conductor {sub.conductor} does not divide {M_doc}"
            )
    irreps: list[ChainIrrep] = []
    for gamma in G.irreps:
        d = gamma.dimension
        lifted = [lift_matrix(G.field, F, D) for D in gamma.matrices]
        spaces = _label_spaces(F, chain, gamma, lifted)
        labels = [label for label, _ in spaces]
        columns: list[list[Value] | None] = [None] * d
        norms: list[Value | None] = [None] * d
        designated: list[int | None] = [None] * d
        if not chain.steps:
            W = [[F.one if i == j else F.zero for j in range(d)] for i in range(d)]
            irreps.append(ChainIrrep(gamma, labels, norms, designated, W, lifted))
            continue
        first = chain.steps[0]
        assert sub is not None
        # blocks by the first label
        blocks: dict[int, list[int]] = {}
        for c, label in enumerate(labels):
            blocks.setdefault(label[0], []).append(c)
        for sigma_index, components in blocks.items():
            sigma_basis = sub.irreps[sigma_index]
            # the sub-chain realisation's components are its labels (sigma_2, ...)
            t = components[0]
            P_t = spaces[t][1]
            j = designated_index(F, G, lifted, P_t)
            N = P_t[j][j]
            mu = norm_factor(F, N, table)
            if mu is None:
                raise MissingNormFactor(
                    N,
                    M_doc,
                    f"{G.spec.symbol} {chain.text} {gamma.symbol.text} {labels[t]}",
                )
            inv_mu = F.inverse(mu)
            w_t = [F.mul(P_t[i][j], inv_mu) for i in range(d)]
            columns[t] = w_t
            norms[t] = N
            designated[t] = j
            for s in components[1:]:
                # P^sigma_{s t} = (dim sigma / |H-dagger|) sum_h conj(D'_sigma(h)_{s t}) D_Gamma(phi_1(h))
                s_sub = sigma_basis.labels.index(labels[s][1:])
                t_sub = sigma_basis.labels.index(labels[t][1:])
                w_s = [F.zero] * d
                for h, g in enumerate(first.phi):
                    sign_h, Dp = matrix_of(sigma_basis.irrep, h, first.lib.order)
                    Dp = sigma_basis.matrices[h % first.lib.order]
                    coeff = F.conjugate(sub.field.lift(Dp[s_sub][t_sub], F))
                    if F.is_zero(coeff):
                        continue
                    if sign_h < 0:
                        coeff = F.neg(coeff)
                    sign_g, _ = matrix_of(gamma, g, G.order)
                    D = lifted[g % G.order]
                    if sign_g < 0:
                        coeff = F.neg(coeff)
                    Dw = [
                        F.total(
                            F.mul(D[r][c], w_t[c])
                            for c in range(d)
                            if not F.is_zero(D[r][c])
                        )
                        for r in range(d)
                    ]
                    w_s = [
                        F.add(a, F.mul(coeff, b)) for a, b in zip(w_s, Dw, strict=True)
                    ]
                scale = Fraction(sigma_basis.dimension, 2 * first.lib.order)
                w_s = [F.scale(x, scale) for x in w_s]
                columns[s] = w_s
                nn = F.total(F.norm_squared(x) for x in w_s)
                if nn != F.one:
                    raise ChainError(
                        f"transferred vector not unit: {gamma.symbol.text} {labels[s]}"
                    )
        filled = [col for col in columns if col is not None]
        assert len(filled) == d
        W = [[filled[c][i] for c in range(d)] for i in range(d)]
        Wd = adjoint(F, W)
        matrices = [mat_mul(F, mat_mul(F, Wd, D), W) for D in lifted]
        irreps.append(ChainIrrep(gamma, labels, norms, designated, W, matrices))
    basis = ChainBasis(chain, M_doc, irreps)
    _BASIS_CACHE[key] = basis
    return basis


# ---- checks (the "Properties" of the chapter) -------------------------------------------


def check_basis(basis: ChainBasis) -> list[str]:
    """Unitarity, the multiplication table with the factor system, block structure
    per step, step-1 blocks equal to the sub-chain realisation, terminal diagonal
    with the characters. Returns the failures."""
    F = basis.field
    chain = basis.chain
    G = chain.G
    failures: list[str] = []
    sub = chain_basis(pulled_back(chain), {}) if chain.steps else None
    for ci in basis.irreps:
        d = ci.dimension
        name = f"{G.spec.symbol} {chain.text} {ci.irrep.symbol.text}"
        # unitarity of W
        WdW = mat_mul(F, adjoint(F, ci.W), ci.W)
        if any(
            WdW[i][j] != (F.one if i == j else F.zero)
            for i in range(d)
            for j in range(d)
        ):
            failures.append(f"{name}: W not unitary")
        # multiplication table with the factor system
        table = G.double_table
        for a in range(G.order):
            for b in range(G.order):
                c = table[a][b]
                sign = -1 if (c >= G.order and ci.irrep.projective) else 1
                lhs = mat_mul(F, ci.matrices[a], ci.matrices[b])
                rhs = ci.matrices[c % G.order]
                if sign < 0:
                    rhs = [[F.neg(x) for x in row] for row in rhs]
                if lhs != rhs:
                    failures.append(f"{name}: D(a)D(b) != D(ab) at ({a},{b})")
                    break
            else:
                continue
            break
        # block structure per step and the terminal characters
        for i, step in enumerate(chain.steps):
            for h, g in enumerate(step.phi):
                sign_g, _ = matrix_of(ci.irrep, g, G.order)
                D = ci.matrices[g % G.order]
                for r in range(d):
                    for c in range(d):
                        if ci.labels[r][: i + 1] != ci.labels[c][
                            : i + 1
                        ] and not F.is_zero(D[r][c]):
                            failures.append(
                                f"{name}: not block-diagonal at step {i + 1}"
                            )
                            break
                    else:
                        continue
                    break
                if i == len(chain.steps) - 1:
                    for r in range(d):
                        sigma = step.lib.irreps[ci.labels[r][-1]]
                        chi = step.lib.field.lift(sigma.character(h, step.lib.order), F)
                        value = D[r][r] if sign_g > 0 else F.neg(D[r][r])
                        if value != chi:
                            failures.append(
                                f"{name}: terminal diagonal is not the character at {h}"
                            )
                            break
        # step-1 blocks equal the sub-chain realisation
        if sub is not None:
            first = chain.steps[0]
            for h, g in enumerate(first.phi):
                if h >= first.lib.order:
                    break
                D = ci.matrices[g % G.order]
                sign_g, _ = matrix_of(ci.irrep, g, G.order)
                for r in range(d):
                    for c in range(d):
                        if ci.labels[r][0] != ci.labels[c][0]:
                            continue
                        sb = sub.irreps[ci.labels[r][0]]
                        rs, cs = (
                            sb.labels.index(ci.labels[r][1:]),
                            sb.labels.index(ci.labels[c][1:]),
                        )
                        value = D[r][c] if sign_g > 0 else F.neg(D[r][c])
                        if value != sub.field.lift(sb.matrices[h][rs][cs], F):
                            failures.append(
                                f"{name}: step-1 block differs from the sub-chain realisation at {h}"
                            )
                            break
                    else:
                        continue
                    break
    return failures


# ---- the subgroup_chain table and the matrices components --------------------------------


def group_record(model: Model) -> dict[str, Any]:
    meta = dict(model.group_metadata())
    meta.pop("input_symbol", None)
    return meta


def subgroup_chain_table(chain: ResolvedChain) -> dict[str, Any]:
    G = chain.G
    N = G.order
    table = G.table
    items = []
    previous_members: list[int] = list(range(N))
    previous_model = G
    previous_phi: list[int] | None = None
    for step in chain.steps:
        H = sorted(step.members)
        Hset = set(H)
        # left cosets g H in K
        cosets: list[list[int]] = []
        seen: set[int] = set()
        for g in previous_members:
            if g in seen:
                continue
            coset = sorted(table[g][h] for h in H)
            seen |= set(coset)
            cosets.append(coset)
        cosets.sort(key=lambda c: c[0])
        # double cosets H g H
        double_cosets = []
        covered: set[int] = set()
        for index, coset in enumerate(cosets):
            if coset[0] in covered:
                continue
            members = {table[h][x] for h in H for x in coset}
            covered |= members
            parts = sorted(i for i, c in enumerate(cosets) if c[0] in members)
            double_cosets.append(
                {"representative": min(members), "size": len(members), "cosets": parts}
            )
        double_cosets.sort(key=lambda r: r["representative"])
        # embedding: G-dagger indices of the lifted preimages of H_lib's canonical elements
        embedding = [step.phi[k] for k in range(step.lib.order)]
        # branching from the previous group's irreps through phi (pulled back when i > 1)
        if previous_phi is None:
            rows = branching(G, step.lib, step.phi)
        else:
            inverse = {g: h for h, g in enumerate(previous_phi)}
            rows = branching(previous_model, step.lib, [inverse[g] for g in step.phi])
        items.append(
            {
                "group": group_record(step.lib),
                "embedding": embedding,
                "cosets": cosets,
                "double_cosets": double_cosets,
                "irreps": [irrep.record() for irrep in step.lib.irreps],
                "branching": [sorted(row) for row in rows],
            }
        )
        assert Hset <= set(previous_members)
        previous_members = H
        previous_model = step.lib
        previous_phi = step.phi
    return {"irreps": [irrep.record() for irrep in G.irreps], "items": items}


def components_of(basis: ChainBasis, ci: ChainIrrep) -> list[dict[str, Any]]:
    out = []
    for k, label in enumerate(ci.labels):
        symbols = [
            basis.chain.steps[j].lib.irreps[idx].symbol.text
            for j, idx in enumerate(label)
        ]
        out.append({"index": k, "labels": list(label), "symbols": symbols})
    return out


def chain_metadata(basis: ChainBasis, input_text: str) -> dict[str, Any]:
    return {
        "input": input_text,
        "canonical": basis.chain.text,
        "symbols": list(basis.chain.symbols),
    }


def field_record(basis: ChainBasis) -> dict[str, Any]:
    return field_metadata(basis.conductor)


__all__ = [
    "ChainBasis",
    "ChainClass",
    "ChainError",
    "ChainIrrep",
    "MissingNormFactor",
    "ResolvedChain",
    "Step",
    "TwinSet",
    "accepted_chains",
    "chain_basis",
    "chain_conductor",
    "chain_metadata",
    "chain_norms",
    "check_basis",
    "components_of",
    "default_chain",
    "differing_step",
    "field_record",
    "label_space_agreement",
    "parse_spelling",
    "pulled_back",
    "reachable_chains",
    "real_sqrt",
    "resolve",
    "resolve_text",
    "same_chain",
    "signature",
    "spelling",
    "subgroup_chain_table",
    "twin_sets",
]
