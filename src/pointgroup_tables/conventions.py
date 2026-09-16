"""The shipped label and basis conventions (docs/spec/label-conventions.md,
docs/spec/basis-conventions.md).

Version 1 ships the canonical basis and the Altmann & Herzig label
convention (the basis convention follows in its own ticket); a name outside
these tables is exit code 4. The rules of a label convention and the
polyhedral generator data of a basis convention (the images of the lifted
generator pair of ``2T``, ``2O``, ``2I`` in every irrep of dimension two or
more) are read from the convention's TOML file on first use.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from functools import cache
from importlib import resources
from typing import Any

from .errors import CONVENTIONS_HINT, UnknownConventionError
from .field import Matrix, cyclotomic_field
from .values import parse

CANONICAL = "canonical"
DATA = resources.files(__package__) / "data" / "conventions"


def _shipped(kind: str) -> dict[str, dict[str, Any]]:
    """``{ name: { name, source } }`` of the convention files under ``data/conventions/<kind>``."""
    table: dict[str, dict[str, Any]] = {}
    for entry in sorted(DATA.joinpath(kind).iterdir(), key=lambda e: e.name):
        if entry.name.endswith(".toml"):
            head = tomllib.loads(entry.read_text(encoding="utf-8"))
            table[head["name"]] = {"name": head["name"], "source": head["source"]}
    return table


# name -> metadata block ({ name, source }); reference conventions add coverage per group
LABEL_CONVENTIONS: dict[str, dict[str, Any]] = _shipped("labels")
BASIS_CONVENTIONS: dict[str, dict[str, Any]] = {CANONICAL: {"name": CANONICAL}}


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
    irreps: dict[str, Any]  # the [irreps] table (spelling rules, letters, overrides)

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
    data = tomllib.loads(
        DATA.joinpath("labels", f"{name}.toml").read_text(encoding="utf-8")
    )
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
        irreps=dict(data.get("irreps", {})),
    )


# ---- polyhedral generator data ------------------------------------------------------


@dataclass(frozen=True)
class PolyhedralData:
    """The ``[polyhedral.<cover>]`` block of a basis convention: the conductor of
    the matrices, the canonical labels of the lifted generator pair (elements of
    the proper group ``T``, ``O`` or ``I``, lifted by the canonical section) and,
    per canonical irrep symbol, the images ``D(q1), D(q2)`` as matrices over
    ``Q(zeta_conductor)``."""

    cover: str
    conductor: int
    generators: tuple[str, str]
    irreps: dict[str, tuple[Matrix, Matrix]]


def _basis_file(name: str) -> dict[str, Any]:
    return tomllib.loads(
        DATA.joinpath("basis", f"{name}.toml").read_text(encoding="utf-8")
    )


@cache
def polyhedral_data(name: str, cover: str) -> PolyhedralData:
    """The polyhedral generator data of a basis convention for ``2T``, ``2O`` or ``2I``."""
    basis_convention(name)
    block = _basis_file(name)["polyhedral"][cover]
    F = cyclotomic_field(block["conductor"])
    irreps: dict[str, tuple[Matrix, Matrix]] = {}
    for symbol, entry in block["irreps"].items():
        first, second = entry["matrices"]
        irreps[symbol] = (
            [[parse(F, text) for text in row] for row in first],
            [[parse(F, text) for text in row] for row in second],
        )
    generators = block["generators"]
    return PolyhedralData(cover, F.conductor, (generators[0], generators[1]), irreps)
