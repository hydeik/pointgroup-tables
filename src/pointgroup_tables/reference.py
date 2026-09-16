"""Reference labels: a label convention applied to the canonical labels of a
model (docs/spec/label-conventions.md, docs/spec/conjugacy-classes.md
"Reference class labels", docs/spec/irrep-labels.md "Reference irrep labels",
docs/spec/double-groups.md "Section overrides").

``ReferenceLabelling`` holds, for one model under one convention, the
reference label of every element of ``G`` and ``G-dagger``, of every class
and of every irrep. ``apply`` writes them into an assembled document next to
the canonical ones: labels are output only, so no table changes shape and no
index moves.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, Any

from .conventions import LabelConvention
from .labels import IDENTITY, INVERSION, REFLECTION, ROTATION, Label

if TYPE_CHECKING:
    from .document import Model
    from .irreps import Irrep

PLUSMINUS = "plusminus"
_POWER = re.compile(r"^([CS])(\d+)\^(-?)(\d+)(\(\d+\))?$")


# ---- spelling rules ---------------------------------------------------------------------


def spell_exponent(label: str, style: str) -> str:
    """``exponent_style``: ``C3^1 -> C3+``, ``C3^-1 -> C3-``, ``C8^3 -> C8^3+``,
    ``S10^-3(2) -> S10^3-(2)``; every other label, and ``"canonical"``, unchanged."""
    match = _POWER.match(label)
    if style != PLUSMINUS or match is None:
        return label
    letter, order, sign, power, axis = match.groups()
    sign = "-" if sign else "+"
    body = (
        f"{letter}{order}{sign}" if power == "1" else f"{letter}{order}^{power}{sign}"
    )
    return body + (axis or "")


def spell_barred(label: str, barred: bool, style: str) -> str:
    """``barred_format``: the reference's spelling of a barred element."""
    return style.replace("{label}", label) if barred else label


def spell_reduced(label: Label, has_inverse: bool, style: str) -> str:
    """The reduced label of a class representative (docs/spec/conjugacy-classes.md)
    spelled by ``exponent_style``; the bar is not part of it."""
    if label.kind in (IDENTITY, INVERSION):
        return label.text.lstrip("~")
    if label.kind == REFLECTION:
        return "sigma" if label.plane == "" else f"sigma_{label.plane}"
    if label.order == 2:
        return f"C2{label.primes}"
    letter = "C" if label.kind == ROTATION else "S"
    assert label.power is not None
    j = abs(label.power)
    body = f"{letter}{label.order}" if j == 1 else f"{letter}{label.order}^{j}"
    if has_inverse or style != PLUSMINUS:
        return body if has_inverse else f"{letter}{label.order}^{label.power}"
    return body + ("-" if label.power < 0 else "+")


# ---- axial secondary directions ------------------------------------------------------------


@dataclass(frozen=True)
class Secondary:
    """The master group ``D_Mh`` of an axial group with ``R = D_N`` and the
    reference's numbering of its secondary directions."""

    N: int
    M: int
    twist: bool  # Cnv with n = 2 (mod 4): direction 0 is the master's C2''1 direction
    primed: tuple[int, ...]
    double_primed: tuple[int, ...]

    def lookup(self, k: int) -> tuple[str, int, int]:
        """``(primes, index, pole)`` of the master's secondary axis in direction ``k``."""
        if self.twist:
            k = (k + self.N // 2) % self.N
        u = k * (self.M // self.N)
        for primes, poles in (("'", self.primed), ("''", self.double_primed)):
            for i, pole in enumerate(poles):
                if pole % self.M == u:
                    return primes, i + 1, pole
        raise ValueError(f"direction {k} has no entry in the D{self.M}h table")

    def swapped(self, k: int) -> bool:
        """Whether the reference's pole of direction ``k`` is the antipode of the
        pole of our lift (docs/spec/label-conventions.md, "Reference section")."""
        ours = (k if 2 * k <= self.N else k - self.N) * (self.M // self.N)
        if self.twist:
            ours += self.M // 2
        _, _, pole = self.lookup(k)
        return (pole - ours) % (2 * self.M) == self.M


def _secondary(model: Model, rules: LabelConvention) -> Secondary | None:
    spec, R = model.spec, model.group.rotation_group
    if R.polyhedral or not R.name.startswith("D"):
        return None
    N = R.axis_order
    assert N is not None and spec.n is not None
    M = N if N % 2 == 0 else 2 * N
    if M not in rules.secondary:
        return None
    primed, double_primed = rules.secondary[M]
    twist = spec.family == "Cnv" and spec.n % 4 == 2
    return Secondary(N, M, twist, primed, double_primed)


# ---- the labelling ------------------------------------------------------------------------


@dataclass(frozen=True)
class ReferenceLabelling:
    """The reference labels of one model under one label convention."""

    model: Model
    rules: LabelConvention

    @cached_property
    def secondary(self) -> Secondary | None:
        return _secondary(self.model, self.rules)

    # ---- elements ---------------------------------------------------------------------------

    @cached_property
    def element_labels(self) -> list[str]:
        """The reference label of every element of ``G``, in canonical order."""
        override = self.rules.overrides.get(self.model.spec.symbol, {})
        out = []
        for element, label in zip(
            self.model.group.elements, self.model.labels, strict=True
        ):
            canonical = label.text
            if canonical in override:
                out.append(override[canonical])
                continue
            direction = self._direction(element)
            if direction is None or self.secondary is None:
                out.append(spell_exponent(canonical, self.rules.exponent_style))
                continue
            primes, index, _ = self.secondary.lookup(direction)
            if label.kind == ROTATION:
                out.append(f"C2{primes}{index}")
            else:
                out.append(f"sigma_{label.plane}{index}")
        return out

    def _direction(self, element: tuple[Any, int]) -> int | None:
        """The secondary direction ``k`` of an element ``((k, 1), parity)`` of an
        axial group with ``R = D_N``; ``None`` for every other element."""
        r, _ = element
        if self.model.group.rotation_group.polyhedral or isinstance(r, int):
            return None
        k, s = r
        return k if s == 1 else None

    @cached_property
    def swapped(self) -> frozenset[int]:
        """The elements of ``G`` whose reference lift is our barred element."""
        listed = self.rules.swapped.get(self.model.spec.symbol)
        if listed is not None:
            return frozenset(
                i for i, label in enumerate(self.model.labels) if label.text in listed
            )
        if self.secondary is None:
            return frozenset()
        return frozenset(
            i
            for i, element in enumerate(self.model.group.elements)
            if (k := self._direction(element)) is not None and self.secondary.swapped(k)
        )

    def is_barred(self, index: int) -> bool:
        """Whether the element ``index`` of ``G-dagger`` is barred in the reference section."""
        N = self.model.order
        return (index >= N) != (index % N in self.swapped)

    @cached_property
    def double_element_labels(self) -> list[str]:
        """The reference label of every element of ``G-dagger``."""
        N = self.model.order
        return [
            spell_barred(
                self.element_labels[i % N], self.is_barred(i), self.rules.barred_format
            )
            for i in range(2 * N)
        ]

    # ---- classes ----------------------------------------------------------------------------

    def _class_label(
        self, size: int, representative: int, has_inverse: bool, double: bool
    ) -> str:
        """A class of one element carries its reference label; a larger class its
        size and the reduced canonical label of its representative spelled by
        the convention. In ``G-dagger`` the bar is the reference section's."""
        N = self.model.order
        if size == 1:
            labels = self.double_element_labels if double else self.element_labels
            return labels[representative]
        reduced = spell_reduced(
            self.model.labels[representative % N],
            has_inverse,
            self.rules.exponent_style,
        )
        barred = double and self.is_barred(representative)
        return str(size) + spell_barred(reduced, barred, self.rules.barred_format)

    @cached_property
    def class_labels(self) -> list[str]:
        override = self.rules.class_overrides.get(self.model.spec.symbol, {})
        inverses = self.model.inverses
        out = []
        for cls in self.model.classes:
            has_inverse = inverses[cls.representative] in cls.members
            label = self._class_label(
                cls.size, cls.representative, has_inverse, double=False
            )
            out.append(override.get(cls.label, label))
        return out

    @cached_property
    def double_class_labels(self) -> list[str]:
        override = self.rules.class_overrides.get(self.model.spec.symbol, {})
        N = self.model.order
        inverses = self.model.double_inverses
        out = []
        for cls in self.model.double_classes:
            rep = cls.representative
            has_inverse = inverses[rep] in cls.members
            if cls.split:
                label = self._class_label(cls.size, rep, has_inverse, double=True)
            else:
                # the two halves joined by "+", the reference's unbarred half first
                unbarred = rep if not self.is_barred(rep) else (rep + N) % (2 * N)
                barred = (unbarred + N) % (2 * N)
                half = cls.size // 2
                label = (
                    self._class_label(half, unbarred, has_inverse, double=True)
                    + "+"
                    + self._class_label(half, barred, has_inverse, double=True)
                )
            out.append(override.get(cls.label, label))
        return out

    # ---- irreps -----------------------------------------------------------------------------

    @cached_property
    def irrep_symbols(self) -> list[str]:
        """The reference symbol of every irrep, in canonical irrep order."""
        rules = self.rules.irreps
        override = rules.get("overrides", {}).get(self.model.spec.symbol, {})
        return [
            override.get(irrep.symbol.text, reference_symbol(irrep, rules))
            for irrep in self.model.irreps
        ]

    # ---- the document -----------------------------------------------------------------------

    def apply(self, document: dict[str, Any]) -> None:
        """Insert the reference labels into an assembled document, in place."""
        for name in ("elements", "cayley", "factor", "characters", "matrices"):
            table = document.get(name)
            if table is None:
                continue
            if "items" in table:
                double = table["group"] == "G_double"
                labels = self.double_element_labels if double else self.element_labels
                table["items"] = [
                    _with(record, "reference_label", labels[record["index"]])
                    for record in table["items"]
                ]
            if "classes" in table:
                double = table.get("class_group", table.get("group")) == "G_double"
                labels = self.double_class_labels if double else self.class_labels
                table["classes"] = [
                    _with(record, "reference_label", labels[record["index"]])
                    for record in table["classes"]
                ]
            if "irreps" in table:
                table["irreps"] = [
                    _with(
                        record, "reference_symbol", self.irrep_symbols[record["index"]]
                    )
                    for record in table["irreps"]
                ]


def reference_symbol(irrep: Irrep, rules: dict[str, Any]) -> str:
    """The spelling rules of the ``[irreps]`` table applied to a canonical symbol."""
    symbol = irrep.symbol
    letters = rules.get("projective_letters" if irrep.projective else "letters", {})
    letter = letters.get(symbol.letter, symbol.letter)
    component = ""
    if symbol.component is not None:
        style = rules.get("component_style", "prefix")
        component = (
            f"^{symbol.component}" if style == "superscript" else str(symbol.component)
        )
    subscript = ""
    if symbol.subscript is not None:
        text = str(symbol.subscript)
        if "/" in text and rules.get("fraction_style", "plain") == "subscript":
            text = f"_{{{text}}}"
        subscript = text
    suffix = symbol.suffix or ""
    if irrep.projective and suffix in ("g", "u"):
        suffix = rules.get("projective_parity_separator", "") + suffix
    return component + letter + subscript + suffix


def _with(record: dict[str, Any], key: str, value: Any) -> dict[str, Any]:
    """The record with ``key`` inserted after the canonical ``label`` / ``symbol`` key."""
    anchor = "label" if key == "reference_label" else "symbol"
    out: dict[str, Any] = {}
    for k, v in record.items():
        out[k] = v
        if k == anchor:
            out[key] = value
    return out
