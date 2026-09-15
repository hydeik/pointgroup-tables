"""Irreps of the axial families in the canonical basis, their characters and
their canonical Mulliken symbols (docs/spec/basis-conventions.md "Canonical
basis: axial families", docs/spec/irrep-labels.md).

An irrep is built from the generator images on the rotation group ``R``
(``a = (1, 0)``, ``b = (0, 1)``), evaluated at the canonical section for the
projective irreps, with the parity sign of a type II group attached. Its
characters are the traces; its symbol is derived from the characters on the
named elements by the rules of the irrep-labels chapter, so that the same
function recomputes it from a character table (``IL-2``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any

from .abstract import CENTRO, PointGroup
from .double import DoubleGroup
from .field import (
    CyclotomicField,
    Matrix,
    Value,
    matrix_identity,
    matrix_mul,
    matrix_neg,
    matrix_trace,
)

SYMBOL = re.compile(r"^([12])?([ABETGHI])(\d+|\d*[13579]/2)?(g|u|'|'')?$")
LETTER_RANK = {letter: rank for rank, letter in enumerate("ABETGHI")}
SUFFIX_RANK = {None: 0, "g": 0, "u": 1, "'": 0, "''": 1}


# ---- symbols -------------------------------------------------------------------------


@dataclass(frozen=True)
class Symbol:
    """A canonical Mulliken symbol decomposed into its fields."""

    letter: str
    subscript: int | str | None = None
    suffix: str | None = None
    component: int | None = None

    @property
    def text(self) -> str:
        return (
            (str(self.component) if self.component else "")
            + self.letter
            + ("" if self.subscript is None else str(self.subscript))
            + (self.suffix or "")
        )

    @property
    def partner(self) -> Symbol | None:
        if self.component is None:
            return None
        return Symbol(self.letter, self.subscript, self.suffix, 3 - self.component)

    def sort_key(self, projective: bool) -> tuple[Any, ...]:
        """The canonical irrep order (docs/spec/irrep-labels.md)."""
        if self.subscript is None:
            value = Fraction(0)
        elif isinstance(self.subscript, str):
            value = Fraction(int(self.subscript.split("/")[0]), 2)
        else:
            value = Fraction(self.subscript)
        return (
            projective,
            SUFFIX_RANK[self.suffix],
            LETTER_RANK[self.letter],
            value,
            self.component or 0,
        )

    @staticmethod
    def parse(text: str) -> Symbol:
        match = SYMBOL.match(text)
        if match is None:
            raise ValueError(f"'{text}' is not a canonical Mulliken symbol")
        component, letter, subscript, suffix = match.groups()
        parsed: int | str | None
        if subscript is None:
            parsed = None
        elif "/" in subscript:
            parsed = subscript
        else:
            parsed = int(subscript)
        return Symbol(letter, parsed, suffix, int(component) if component else None)


# ---- the named elements --------------------------------------------------------------


@dataclass(frozen=True)
class NamedElements:
    """The elements the symbol rules refer to, as indices into the canonical
    element order of ``G`` (docs/spec/irrep-labels.md, "Named elements")."""

    principal: int  # c = (1, 0) with its parity; order N
    proper_principal: int  # c0: c, or (2, 0) for Cnh / Dnh with odd n
    order: int  # N, the order of c
    proper_order: int  # the order of c0
    secondary: int | None  # s = (0, 1) with its parity
    double_prime: int | None  # (1, 1) of D2: C2''(1), rule 6 only
    inversion: int | None
    horizontal: int | None  # sigma_h of Cnh / Dnh with odd n
    d2: bool  # D2 and D2h: rule 6


def named_elements(G: PointGroup) -> NamedElements:
    R = G.rotation_group
    N = R.axis_order
    if N is None:
        raise ValueError("the named elements of the polyhedral groups are not axial")
    dihedral = R.name.startswith("D")
    family, n = G.spec.family, G.spec.family_n
    odd_h = family in ("Cnh", "Dnh") and n % 2 == 1

    def find(r: Any, parity: int | None = None) -> int | None:
        for index, (element, p) in enumerate(G.elements):
            if element == r and (parity is None or p == parity):
                return index
        return None

    def rotation(k: int, s: int = 0) -> Any:
        return (k % N, s) if dihedral else k % N

    def element(r: Any) -> int:
        proper = find(r, +1)
        index = proper if proper is not None else find(r)
        assert index is not None, r
        return index

    principal = element(rotation(1))
    proper_principal = element(rotation(2)) if odd_h else principal
    secondary = element(rotation(0, 1)) if dihedral else None
    inversion = find(R.identity, -1) if G.type == CENTRO else None
    horizontal = find(rotation(N // 2), -1) if odd_h else None
    d2 = family in ("Dn", "Dnh") and n == 2  # rule 6: D2 and D2h only
    return NamedElements(
        principal=principal,
        proper_principal=proper_principal,
        order=N,
        proper_order=n if odd_h else N,
        secondary=secondary,
        double_prime=element(rotation(1, 1)) if d2 else None,
        inversion=inversion,
        horizontal=horizontal,
        d2=d2,
    )


# ---- symbol rules ----------------------------------------------------------------------


def _exponent(
    F: CyclotomicField, x: Value, m: int, candidates: range
) -> tuple[int, int] | None:
    """``(k, sign)`` with ``x = zeta_m^(sign k)`` for ``k`` in ``candidates``, or ``None``."""
    for k in candidates:
        if x == F.root_of_unity(m, k):
            return k, +1
        if x == F.root_of_unity(m, -k):
            return k, -1
    return None


def _two_cosine(F: CyclotomicField, x: Value, m: int, candidates: range) -> int | None:
    """``k`` with ``x = zeta_m^k + zeta_m^-k`` for ``k`` in ``candidates``, or ``None``."""
    for k in candidates:
        if x == F.add(F.root_of_unity(m, k), F.root_of_unity(m, -k)):
            return k
    return None


def derive_symbol(
    F: CyclotomicField,
    named: NamedElements,
    chi: list[Value],
    projective: bool,
) -> Symbol:
    """The canonical Mulliken symbol of an irrep from its characters on the
    elements of ``G`` (``chi[g]``; for a projective irrep the characters of the
    lifted elements), by the axial rules 1-11 of the irrep-labels chapter."""
    one = F.one
    minus_one = F.neg(one)
    d = int(F.as_rational(chi[0]))

    def is_plus(index: int | None) -> bool:
        assert index is not None
        return chi[index] == F.rational(d)

    if not projective:
        N0 = named.proper_order
        x = chi[named.proper_principal]
        component: int | None = None
        subscript: int | str | None
        if d == 1:
            if x == one:
                letter = "A"
            elif x == minus_one:
                letter = "B"
            else:
                found = _exponent(F, x, N0, range(1, (N0 + 1) // 2))
                if found is None:
                    raise ValueError(
                        "a one-dimensional character off the roots of unity"
                    )
                letter = "E"
                k, sign = found
                component = 1 if sign > 0 else 2
        elif d == 2:
            k = _two_cosine(F, x, N0, range(1, (N0 + 1) // 2))
            if k is None:
                raise ValueError("a two-dimensional character off 2 cos(2 pi k / N)")
            letter = "E"
        else:
            raise ValueError(f"an axial irrep of dimension {d}")
        if letter == "E":
            subscript = None if N0 <= 4 else k
        elif named.d2:
            # rule 6: the three non-trivial irreps of D2 are B1, B2, B3 by the
            # binary rotation they are symmetric under: C2, C2''(1), C2'(1)
            if chi[named.principal] == one and chi[named.secondary or 0] == one:
                subscript = None
            elif chi[named.principal] == one:
                letter, subscript = "B", 1
            elif chi[named.double_prime or 0] == one:
                letter, subscript = "B", 2
            else:
                letter, subscript = "B", 3
        elif named.secondary is not None:
            subscript = 1 if chi[named.secondary] == one else 2
        else:
            subscript = None
        suffix: str | None = None
        if named.inversion is not None:
            suffix = "g" if is_plus(named.inversion) else "u"
        elif named.horizontal is not None:
            suffix = "'" if is_plus(named.horizontal) else "''"
        return Symbol(letter, subscript, suffix, component)

    N = named.order
    x = chi[named.principal]
    component = None
    if d == 2:
        k = _two_cosine(F, x, 2 * N, range(1, N, 2))
        if k is None:
            raise ValueError(
                "a two-dimensional projective character off 2 cos(pi k / N)"
            )
        letter, subscript = "E", f"{k}/2"
    elif d == 1:
        if x == one:
            letter, subscript = "A", "1/2"
        elif x == minus_one:
            subscript = f"{N}/2"
            if named.secondary is not None:
                letter = "E"
                y = chi[named.secondary]
                if y == F.constant("i"):
                    component = 1
                elif y == F.neg(F.constant("i")):
                    component = 2
                else:
                    raise ValueError(
                        "the secondary character of a k = N pair is not +-i"
                    )
            else:
                letter = "B"
        else:
            found = _exponent(F, x, 2 * N, range(1, N, 2))
            if found is None:
                raise ValueError(
                    "a one-dimensional projective character off the roots of unity"
                )
            k, sign = found
            letter, subscript = "E", f"{k}/2"
            component = 1 if sign > 0 else 2
    else:
        raise ValueError(f"an axial projective irrep of dimension {d}")
    suffix = None
    if named.inversion is not None:
        suffix = "g" if is_plus(named.inversion) else "u"
    return Symbol(letter, subscript, suffix, component)


# ---- irreps ---------------------------------------------------------------------------


@dataclass
class Irrep:
    """One irrep of ``G``: its matrices on the elements of ``G`` in the canonical
    basis (on the lifted elements for a projective irrep, with the parity sign
    of a type II group attached), its characters and its symbol."""

    symbol: Symbol
    projective: bool
    matrices: list[Matrix]
    characters: list[Value]
    index: int = -1
    partner: int | None = None
    frobenius_schur: int = 0
    bethe: str | None = None
    reference_symbol: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def dimension(self) -> int:
        return len(self.matrices[0])

    def character(self, index: int, order: int) -> Value:
        """The character at an element of ``G`` (``index < |G|``) or of ``G-dagger``
        (a barred index ``>= |G|`` negates a projective character)."""
        value = self.characters[index % order]
        if index >= order and self.projective:
            return tuple(-a for a in value)
        return value

    def record(self) -> dict[str, Any]:
        record: dict[str, Any] = {"index": self.index, "symbol": self.symbol.text}
        if self.reference_symbol is not None:
            record["reference_symbol"] = self.reference_symbol
        if self.bethe is not None:
            record["bethe"] = self.bethe
        record["dimension"] = self.dimension
        record["projective"] = self.projective
        record["letter"] = self.symbol.letter
        if self.symbol.component is not None:
            record["component"] = self.symbol.component
            record["partner"] = self.partner
        if self.symbol.subscript is not None:
            record["subscript"] = self.symbol.subscript
        if self.symbol.suffix is not None:
            record["suffix"] = self.symbol.suffix
        record["frobenius_schur"] = self.frobenius_schur
        return record


def _generator_images(
    F: CyclotomicField, N: int, dihedral: bool
) -> list[tuple[bool, Matrix, Matrix | None]]:
    """``(projective, D(a), D(b))`` of every irrep of ``R`` or ``R-dagger``
    (``D(b)`` is ``None`` for a cyclic ``R``), from the closed forms."""
    one = [[F.one]]
    minus = [[F.neg(F.one)]]
    i = F.constant("i") if dihedral else None
    images: list[tuple[bool, Matrix, Matrix | None]] = []
    if dihedral:
        images.append((False, one, one))
        images.append((False, one, minus))
        if N % 2 == 0:
            images.append((False, minus, one))
            images.append((False, minus, minus))
        for k in range(1, N):
            if 2 * k >= N:
                break
            Da = [[F.root_of_unity(N, k), F.zero], [F.zero, F.root_of_unity(N, -k)]]
            Db = [[F.zero, F.one], [F.one, F.zero]]
            images.append((False, Da, Db))
        assert i is not None
        for k in range(1, N, 2):
            Da = [
                [F.root_of_unity(2 * N, k), F.zero],
                [F.zero, F.root_of_unity(2 * N, -k)],
            ]
            Db = [[F.zero, i], [i, F.zero]]
            images.append((True, Da, Db))
        if N % 2 == 1:
            images.append((True, minus, [[i]]))
            images.append((True, minus, [[F.neg(i)]]))
    else:
        for k in range(N):
            images.append((False, [[F.root_of_unity(N, k)]], None))
        for k in range(1, 2 * N, 2):
            images.append((True, [[F.root_of_unity(2 * N, k)]], None))
    return images


def build_irreps(
    G: PointGroup,
    D: DoubleGroup,
    F: CyclotomicField,
    table: list[list[int]],
    factor: list[list[int]],
) -> list[Irrep]:
    """Every irrep of an axial group in canonical irrep order, with symbols,
    partners and Frobenius-Schur indicators."""
    R = G.rotation_group
    N = R.axis_order
    if N is None:
        raise ValueError("build_irreps handles the axial families only")
    dihedral = R.name.startswith("D")
    named = named_elements(G)
    # the lifted rotation part of every element, as a word a^k b^s in R-dagger
    words: list[tuple[int, int]] = []
    for r, _ in G.elements:
        q = D.cover.section(r)
        words.append((q, 0) if isinstance(q, int) else (q[0], q[1]))
    plain_words: list[tuple[int, int]] = [
        (r, 0) if isinstance(r, int) else (r[0], r[1]) for r, _ in G.elements
    ]
    parities = [p for _, p in G.elements]
    suffix_signs = [1] if G.type != CENTRO else [1, -1]

    irreps: list[Irrep] = []
    for projective, Da, Db in _generator_images(F, N, dihedral):
        d = len(Da)
        identity = matrix_identity(F, d)
        period = 2 * N if projective else N
        powers_a: list[Matrix] = [identity]
        for _ in range(1, period):
            powers_a.append(matrix_mul(F, powers_a[-1], Da))
        base: list[Matrix] = []
        for k, s in words if projective else plain_words:
            matrix = powers_a[k % period]
            if s and Db is not None:
                matrix = matrix_mul(F, matrix, Db)
            base.append(matrix)
        for sign in suffix_signs:
            matrices = [
                matrix if (p == 1 or sign == 1) else matrix_neg(F, matrix)
                for matrix, p in zip(base, parities, strict=True)
            ]
            characters = [matrix_trace(F, matrix) for matrix in matrices]
            symbol = derive_symbol(F, named, characters, projective)
            irreps.append(Irrep(symbol, projective, matrices, characters))

    irreps.sort(key=lambda irrep: irrep.symbol.sort_key(irrep.projective))
    by_symbol = {irrep.symbol.text: k for k, irrep in enumerate(irreps)}
    if len(by_symbol) != len(irreps):
        raise ValueError("two irreps received the same canonical symbol")
    for k, irrep in enumerate(irreps):
        irrep.index = k
        partner = irrep.symbol.partner
        if partner is not None:
            irrep.partner = by_symbol[partner.text]
        irrep.frobenius_schur = frobenius_schur(F, irrep, table, factor)
    return irreps


def frobenius_schur(
    F: CyclotomicField, irrep: Irrep, table: list[list[int]], factor: list[list[int]]
) -> int:
    """``(1/|G|) sum_g omega(g, g) chi(g^2)``, with ``omega = 1`` for an ordinary irrep."""
    order = len(table)
    total = F.total(
        F.neg(irrep.characters[table[g][g]])
        if irrep.projective and factor[g][g] == -1
        else irrep.characters[table[g][g]]
        for g in range(order)
    )
    value = F.scale(total, Fraction(1, order))
    if not F.is_rational(value) or value[0].denominator != 1:
        raise ValueError(f"the Frobenius-Schur sum of {irrep.symbol.text} is {value}")
    return int(value[0])


def select(irreps: list[Irrep], which: str) -> list[Irrep]:
    """The rows of ``--irreps ordinary|projective|all``; indices stay global."""
    if which == "all":
        return list(irreps)
    projective = which == "projective"
    return [irrep for irrep in irreps if irrep.projective == projective]
