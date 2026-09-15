"""Canonical element labels (docs/spec/label-conventions.md, canonical rules 1-7).

A label is derived by a fixed formula from ``(r, parity)``, the family, the
group type and ``H``. It is kept as a structure (kind, order, power, axis) and
rendered to the ASCII grammar on demand, so element records never parse text.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .abstract import PointGroup

IDENTITY, ROTATION, INVERSION, REFLECTION, ROTOREFLECTION = (
    "identity",
    "rotation",
    "inversion",
    "reflection",
    "rotoreflection",
)
PRINCIPAL, SECONDARY, SECONDARY_DOUBLE_PRIME, POLYHEDRAL = (
    "principal",
    "secondary",
    "secondary_double_prime",
    "polyhedral",
)


@dataclass(frozen=True)
class Label:
    """The decomposed canonical label of an element."""

    kind: str
    order: int
    power: int | None = None
    axis_class: str | None = None
    axis_index: int | None = None
    plane: str | None = None  # "h", "v", "d" or "" for reflections
    primed: bool = False  # C2'(k) of a polyhedral group (rule 6)
    barred: bool = False

    @property
    def primes(self) -> str:
        if self.primed:
            return "'"
        return {SECONDARY: "'", SECONDARY_DOUBLE_PRIME: "''"}.get(
            self.axis_class or "", ""
        )

    @property
    def text(self) -> str:
        return ("~" if self.barred else "") + self._base()

    def _base(self) -> str:
        suffix = f"({self.axis_index})" if self.axis_index is not None else ""
        if self.kind == IDENTITY:
            return "E"
        if self.kind == INVERSION:
            return "i"
        if self.kind == REFLECTION:
            letter = "sigma" if self.plane == "" else f"sigma_{self.plane}"
            return letter + suffix
        letter = "C" if self.kind == ROTATION else "S"
        if self.order == 2 and self.kind == ROTATION:
            return f"C2{self.primes}{suffix}"
        return f"{letter}{self.order}^{self.power}{suffix}"


def reduced_turn(turn: Fraction) -> tuple[int, int]:
    """``turn`` in ``(0, 1)`` -> ``(m, j)`` in lowest terms with ``|j| <= m/2``."""
    m, j = turn.denominator, turn.numerator
    if 2 * j > m:
        j -= m
    return m, j


def _rotation(turn: Fraction, axis_class: str, axis_index: int | None = None) -> Label:
    m, j = reduced_turn(turn)
    if m == 2:
        return Label(ROTATION, 2, axis_class=axis_class, axis_index=axis_index)
    return Label(ROTATION, m, power=j, axis_class=axis_class, axis_index=axis_index)


def _improper(turn: Fraction, axis_class: str, axis_index: int | None = None) -> Label:
    """The improper element ``i * C(turn)``: rule 3 and, for polyhedral groups, rule 7."""
    g = (turn + Fraction(1, 2)) % 1
    if g == 0:
        return Label(
            REFLECTION, 2, axis_class=axis_class, axis_index=axis_index, plane="h"
        )
    if g == Fraction(1, 2):
        return Label(INVERSION, 2)
    m, j = reduced_turn(g)
    return Label(
        ROTOREFLECTION, m, power=j, axis_class=axis_class, axis_index=axis_index
    )


def _secondary(k_prime: int, n_prime: int) -> Label:
    """Rule 2: the secondary axis with index ``k'`` within ``G+`` of axis order ``n'``."""
    if n_prime % 2:
        return Label(ROTATION, 2, axis_class=SECONDARY, axis_index=k_prime + 1)
    if k_prime % 2 == 0:
        return Label(ROTATION, 2, axis_class=SECONDARY, axis_index=k_prime // 2 + 1)
    return Label(
        ROTATION,
        2,
        axis_class=SECONDARY_DOUBLE_PRIME,
        axis_index=(k_prime - 1) // 2 + 1,
    )


def _plane(letter: str, index: int) -> Label:
    axis_class = SECONDARY if letter == "v" else SECONDARY_DOUBLE_PRIME
    return Label(REFLECTION, 2, axis_class=axis_class, axis_index=index, plane=letter)


def _axial_labels(G: PointGroup) -> list[Label]:
    family, n = G.spec.family, G.spec.n
    assert n is not None
    N = G.rotation_group.axis_order
    assert N is not None
    if family == "S2n":
        n_prime = n // 2
        step = 1 if n_prime % 2 else 2
    else:
        n_prime = n
        step = N // n
    labels = []
    for r, parity in G.elements:
        k, s = (r, 0) if isinstance(r, int) else r
        turn = Fraction(k, N)
        if parity == +1 and s == 0:
            labels.append(Label(IDENTITY, 1) if k == 0 else _rotation(turn, PRINCIPAL))
        elif parity == +1:
            labels.append(_secondary(k // step, n_prime))
        elif s == 0:
            labels.append(_improper(turn, PRINCIPAL))
        else:
            labels.append(_vertical_plane(family, n, N, step, n_prime, k))
    return labels


def _vertical_plane(
    family: str, n: int, N: int, step: int, n_prime: int, k: int
) -> Label:
    """Rule 4: the plane ``(k, 1, -1)`` perpendicular to secondary direction ``k``."""
    if family == "Cnv":
        if n % 2:
            return _plane("v", k + 1)
        return _plane("v", k // 2 + 1) if k % 2 == 0 else _plane("d", (k - 1) // 2 + 1)
    if family == "Dnh":
        # the plane contains direction d = k + N/2; it is named by the axis there
        d = (k + N // 2) % N
        axis = _secondary(d // step, n_prime)
        assert axis.axis_index is not None
        return _plane("v" if axis.axis_class == SECONDARY else "d", axis.axis_index)
    if family == "Dnd":
        return _plane("d", k + 1) if n % 2 else _plane("d", (k - 1) // 2 + 1)
    raise ValueError(f"family {family} has no vertical planes")


def _polyhedral_labels(G: PointGroup) -> list[Label]:
    family = G.spec.family
    R = G.rotation_group
    elements = R.elements
    identity = R.identity
    order = {r: R.order_of(r) for r in elements}

    def cyclic_subgroup(p: Any) -> list[Any]:
        out, q = [identity], p
        while q != identity:
            out.append(q)
            q = R.mul(q, p)
        return out

    # rule 5: the axis of r is the maximal cyclic subgroup containing r
    axis_of: dict[Any, frozenset[Any]] = {}
    for r in elements[1:]:
        best = None
        for g in elements[1:]:
            if r in cyclic_subgroup(g) and (best is None or order[g] > order[best]):
                best = g
        assert best is not None
        axis_of[r] = frozenset(cyclic_subgroup(best))
    axis_order = {axis: max(order[x] for x in axis) for axis in set(axis_of.values())}
    # axes numbered per maximal order, in order of first appearance
    axis_number: dict[frozenset[Any], int] = {}
    for r in elements[1:]:
        axis = axis_of[r]
        if axis not in axis_number:
            same = sum(1 for a in axis_number if axis_order[a] == axis_order[axis])
            axis_number[axis] = same + 1
    # the positive generator: the maximal-order element with the smallest index
    position = {r: i for i, r in enumerate(elements)}
    generator = {
        axis: min(
            (x for x in axis if order[x] == axis_order[axis]), key=lambda x: position[x]
        )
        for axis in axis_number
    }
    has_c2_on_higher_axes = any(
        order[r] == 2 and axis_order[axis_of[r]] > 2 for r in elements[1:]
    )

    def turn_and_axis(r: Any) -> tuple[Fraction, int, bool]:
        axis = axis_of[r]
        powers = cyclic_subgroup(generator[axis])
        turn = Fraction(powers.index(r), axis_order[axis])
        on_binary_axis = axis_order[axis] == 2
        return turn, axis_number[axis], on_binary_axis

    labels = []
    for r, parity in G.elements:
        if r == identity:
            labels.append(Label(IDENTITY, 1) if parity == +1 else Label(INVERSION, 2))
            continue
        turn, index, on_binary_axis = turn_and_axis(r)
        if parity == +1:
            # rule 6: C2'(k) on an order-2 axis when higher axes also carry C2;
            # every axis of a polyhedral group is a polyhedral axis
            if order[r] == 2 and on_binary_axis and has_c2_on_higher_axes:
                labels.append(
                    Label(
                        ROTATION,
                        2,
                        axis_class=POLYHEDRAL,
                        axis_index=index,
                        primed=True,
                    )
                )
            else:
                labels.append(_rotation(turn, POLYHEDRAL, index))
        elif order[r] == 2:
            labels.append(
                _polyhedral_plane(
                    family, index, on_binary_axis and has_c2_on_higher_axes
                )
            )
        else:
            labels.append(_improper(turn, POLYHEDRAL, index))
    return labels


def _polyhedral_plane(family: str, index: int, on_secondary_axis: bool) -> Label:
    """Rule 7: reflections named by family."""
    if family == "Th":
        plane = "h"
    elif family == "Oh":
        plane = "d" if on_secondary_axis else "h"
    elif family == "Td":
        plane = "d"
    elif family == "Ih":
        plane = ""
    else:
        raise ValueError(f"family {family} has no reflections")
    return Label(REFLECTION, 2, axis_class=POLYHEDRAL, axis_index=index, plane=plane)


def canonical_labels(G: PointGroup) -> list[Label]:
    """The canonical label of every element of ``G``, in canonical order."""
    return _polyhedral_labels(G) if G.rotation_group.polyhedral else _axial_labels(G)


def reduced_label(label: Label, has_inverse: bool) -> str:
    """The reduced form of a label for a class label (docs/spec/conjugacy-classes.md)."""
    bar = "~" if label.barred else ""
    if label.kind in (IDENTITY, INVERSION):
        return bar + label._base()
    if label.kind == REFLECTION:
        return bar + ("sigma" if label.plane == "" else f"sigma_{label.plane}")
    if label.order == 2:
        return f"{bar}C2{label.primes}"
    letter = "C" if label.kind == ROTATION else "S"
    assert label.power is not None
    j = abs(label.power) if has_inverse else label.power
    return (
        f"{bar}{letter}{label.order}" if j == 1 else f"{bar}{letter}{label.order}^{j}"
    )
