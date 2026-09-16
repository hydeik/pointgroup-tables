"""The shipped label and basis conventions (docs/spec/label-conventions.md,
docs/spec/basis-conventions.md).

Version 1 ships the canonical basis and the Altmann & Herzig label and basis
conventions; a name outside these tables is exit code 4. The rules of a label
convention, the axial rule and coverage of a basis convention and its
polyhedral generator data (the images of a lifted generating set of ``2T``,
``2O``, ``2I`` in every irrep of dimension two or more) are read from the
convention's TOML file on first use. How a basis convention is applied to a
group (its conductor, the axial generator images, the lifted generators) is
``basis.py``.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from functools import cache
from importlib import resources
from importlib.resources.abc import Traversable
from typing import Any

from .errors import CONVENTIONS_HINT, UnknownConventionError
from .field import Matrix, cyclotomic_field
from .values import parse

CANONICAL = "canonical"
DATA = resources.files(__package__) / "data" / "conventions"


def _shipped(kind: str) -> dict[str, Traversable]:
    """``{ name: file }`` of the convention files under ``data/conventions/<kind>``,
    the canonical convention first, the others by name."""
    files: dict[str, Traversable] = {}
    for entry in DATA.joinpath(kind).iterdir():
        if entry.name.endswith(".toml"):
            head = tomllib.loads(entry.read_text(encoding="utf-8"))
            files[head["name"]] = entry
    return {
        name: files[name]
        for name in sorted(files, key=lambda name: (name != CANONICAL, name))
    }


def _headers(files: dict[str, Traversable]) -> dict[str, dict[str, Any]]:
    """``{ name: { name, source } }`` (no ``source`` for the canonical basis)."""
    table: dict[str, dict[str, Any]] = {}
    for name, entry in files.items():
        head = tomllib.loads(entry.read_text(encoding="utf-8"))
        table[name] = {"name": head["name"]}
        if "source" in head:
            table[name]["source"] = head["source"]
    return table


# name -> the convention's file
LABEL_FILES: dict[str, Traversable] = _shipped("labels")
BASIS_FILES: dict[str, Traversable] = _shipped("basis")
# name -> metadata block ({ name, source }); reference conventions add coverage per group
LABEL_CONVENTIONS: dict[str, dict[str, Any]] = _headers(LABEL_FILES)
BASIS_CONVENTIONS: dict[str, dict[str, Any]] = _headers(BASIS_FILES)


def label_convention(name: str) -> dict[str, Any]:
    if name not in LABEL_CONVENTIONS:
        raise UnknownConventionError(
            f"unknown label convention '{name}'", CONVENTIONS_HINT
        )
    return dict(LABEL_CONVENTIONS[name])


def basis_convention(name: str) -> dict[str, Any]:
    if name not in BASIS_CONVENTIONS:
        raise UnknownConventionError(
            f"unknown basis convention '{name}'", CONVENTIONS_HINT
        )
    return dict(BASIS_CONVENTIONS[name])


# ---- label convention rules ------------------------------------------------------------


@dataclass(frozen=True)
class IrrepRules:
    """The ``[irreps]`` table of a label convention (docs/spec/irrep-labels.md,
    "Reference irrep labels in a label convention")."""

    component_style: str = "prefix"  # "prefix" (1E) or "superscript" (^1E)
    fraction_style: str = "plain"  # "plain" (E1/2) or "subscript" (E_{1/2})
    projective_parity_separator: str = ""  # "," gives E1/2,g
    letters: dict[str, str] = field(default_factory=dict)  # ordinary irreps
    projective_letters: dict[str, str] = field(default_factory=dict)
    overrides: dict[str, dict[str, str]] = field(default_factory=dict)  # per group

    @staticmethod
    def from_table(table: dict[str, Any]) -> IrrepRules:
        return IrrepRules(
            component_style=table.get("component_style", "prefix"),
            fraction_style=table.get("fraction_style", "plain"),
            projective_parity_separator=table.get("projective_parity_separator", ""),
            letters=dict(table.get("letters", {})),
            projective_letters=dict(table.get("projective_letters", {})),
            overrides={
                group: dict(entries)
                for group, entries in table.get("overrides", {}).items()
            },
        )


@dataclass(frozen=True)
class LabelConvention:
    """The rule vocabulary of a label convention file (docs/spec/label-conventions.md,
    "Label conventions (reference labels)"), as read from the shipped TOML."""

    name: str
    source: str
    exponent_style: str  # "plusminus" or "canonical"
    barred_format: str  # "~{label}" or "{label}~"
    axial_max_n: int  # [coverage]: every axial family up to this n is verified
    polyhedral: bool  # [coverage]: whether the polyhedral groups are verified
    # [secondary."M"]: the poles of the master group's C2'1, C2'2, ... and
    # C2''1, C2''2, ... in units of 180/M degrees
    secondary: dict[int, tuple[tuple[int, ...], tuple[int, ...]]]
    overrides: dict[str, dict[str, str]]  # [overrides.<group>]: canonical -> reference
    swapped: dict[str, frozenset[str]]  # [barred.<group>].swapped
    class_overrides: dict[str, dict[str, str]]  # [classes.overrides.<group>]
    irreps: IrrepRules

    def covers(self, family_n: int | None) -> bool:
        """Whether a group (axial by its family ``n``, polyhedral for ``None``) is
        inside the range the reference prints."""
        if family_n is None:
            return self.polyhedral
        return family_n <= self.axial_max_n


@cache
def label_convention_rules(name: str) -> LabelConvention:
    """The rules of a shipped label convention, read once from its file."""
    label_convention(name)
    data = tomllib.loads(LABEL_FILES[name].read_text(encoding="utf-8"))
    coverage = data.get("coverage", {})
    return LabelConvention(
        name=data["name"],
        source=data["source"],
        exponent_style=data.get("exponent_style", "canonical"),
        barred_format=data.get("barred_format", "~{label}"),
        axial_max_n=int(coverage.get("axial_max_n", 0)),
        polyhedral=bool(coverage.get("polyhedral", False)),
        secondary={
            int(M): (tuple(entry["primed"]), tuple(entry["double_primed"]))
            for M, entry in data.get("secondary", {}).items()
        },
        overrides={
            group: dict(table) for group, table in data.get("overrides", {}).items()
        },
        swapped={
            group: frozenset(table.get("swapped", ()))
            for group, table in data.get("barred", {}).items()
        },
        class_overrides={
            group: dict(table)
            for group, table in data.get("classes", {}).get("overrides", {}).items()
        },
        irreps=IrrepRules.from_table(data.get("irreps", {})),
    )


# ---- basis conventions -----------------------------------------------------------------


@dataclass(frozen=True)
class PolyhedralData:
    """The ``[polyhedral.<cover>]`` block of a basis convention: the conductor of
    the matrices, the canonical labels of a lifted generating set (elements of
    the proper group ``T``, ``O`` or ``I``, lifted by the canonical section; a
    ``~`` prefix selects the barred lift) and, per canonical irrep symbol, the
    images ``D(q1), D(q2), ...`` as matrices over ``Q(zeta_conductor)``."""

    cover: str
    conductor: int
    generators: tuple[str, ...]
    irreps: dict[str, tuple[Matrix, ...]]


@dataclass(frozen=True)
class BasisConvention:
    """The rules of a basis convention file (docs/spec/basis-conventions.md,
    "Reference basis conventions"), as read from the shipped TOML: the axial
    rule (``None`` for a file without an ``[axial]`` table, which does not cover
    the axial families), the covers with polyhedral data, and the reference's
    range for the ``coverage`` metadata."""

    name: str
    source: str | None
    axial_rule: str | None  # "canonical", "transform" or "spin-harmonic"
    anchors: (
        str | None
    )  # spin-harmonic: the label convention whose anchoring fixes alpha
    covers: frozenset[str]  # the R-dagger with a [polyhedral.<cover>] block
    axial_max_n: int  # [coverage]: every axial family up to this n is verified
    polyhedral: bool  # [coverage]: whether the polyhedral groups are verified

    @property
    def reference(self) -> bool:
        """Whether this is a reference convention (the canonical basis is not)."""
        return self.name != CANONICAL

    def covers_group(self, family_n: int | None, cover: str | None) -> bool:
        """Whether the convention realises a group at all (axial by its family
        ``n`` with ``cover = None``, polyhedral by its ``R-dagger``)."""
        if cover is None:
            return self.axial_rule is not None
        return cover in self.covers

    def verified(self, family_n: int | None) -> bool:
        """Whether a group (axial by its family ``n``, polyhedral for ``None``) is
        inside the range the reference prints (``coverage = "reference"``)."""
        if family_n is None:
            return self.polyhedral
        return family_n <= self.axial_max_n


def _basis_file(name: str) -> dict[str, Any]:
    basis_convention(name)
    return tomllib.loads(BASIS_FILES[name].read_text(encoding="utf-8"))


@cache
def basis_convention_rules(name: str) -> BasisConvention:
    """The rules of a shipped basis convention, read once from its file (the
    polyhedral matrices are parsed separately, per cover, by ``polyhedral_data``)."""
    data = _basis_file(name)
    axial = data.get("axial")
    coverage = data.get("coverage", {})
    rule = None if axial is None else axial["rule"]
    if rule not in (None, "canonical", "spin-harmonic"):
        raise ValueError(f"basis convention '{name}': unsupported axial rule '{rule}'")
    anchors = None if axial is None else axial.get("anchors")
    if rule == "spin-harmonic" and anchors not in LABEL_CONVENTIONS:
        raise ValueError(
            f"basis convention '{name}': the spin-harmonic rule must anchor to a shipped label convention"
        )
    return BasisConvention(
        name=data["name"],
        source=data.get("source"),
        axial_rule=rule,
        anchors=anchors,
        covers=frozenset(data.get("polyhedral", {})),
        axial_max_n=int(coverage.get("axial_max_n", 0)),
        polyhedral=bool(coverage.get("polyhedral", False)),
    )


@cache
def polyhedral_data(name: str, cover: str) -> PolyhedralData:
    """The polyhedral generator data of a basis convention for ``2T``, ``2O`` or ``2I``."""
    block = _basis_file(name)["polyhedral"][cover]
    F = cyclotomic_field(block["conductor"])
    generators = tuple(block["generators"])
    irreps: dict[str, tuple[Matrix, ...]] = {}
    for symbol, entry in block["irreps"].items():
        matrices = entry["matrices"]
        if len(matrices) != len(generators):
            raise ValueError(
                f"basis convention '{name}': {cover} irrep {symbol} has "
                f"{len(matrices)} matrices for {len(generators)} generators"
            )
        irreps[symbol] = tuple(
            [[parse(F, text) for text in row] for row in matrix] for matrix in matrices
        )
    return PolyhedralData(cover, F.conductor, generators, irreps)
