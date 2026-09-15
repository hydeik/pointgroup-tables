"""The document model: everything one invocation computes, and its projection
to the document of one subcommand (docs/spec/output-schema.md).

``build_document(request)`` is the first seam of the tool: it names the
subcommand, the group as typed and every option, and returns the dict exactly
as serialised. ``model_document`` returns the same dict with every table the
model can produce, the input of the checkers.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import cached_property
from importlib import metadata
from typing import Any

from . import conventions
from .abstract import PointGroup, point_group
from .classes import (
    ConjugacyClass,
    DoubleClass,
    conjugacy_classes,
    double_conjugacy_classes,
    inverses,
)
from .double import DoubleGroup, double_group
from .errors import UsageError
from .field import canonical_conductor, field_metadata
from .labels import Label, canonical_labels
from .symbols import GroupSpec, parse_group

SCHEMA_VERSION = 1
GENERATOR_NAME = "pointgroup-tables"
TABLE_ORDER = ("elements", "cayley", "factor", "characters", "matrices")
VALUE_FORMS = ("exact", "numeric", "pretty")
DEFAULT_VALUES = ("exact", "numeric")

# the tables each table subcommand puts in its document, in table order
TABLES_OF_SUBCOMMAND = {
    "elements": ("elements",),
    "cayley": ("elements", "cayley"),
    "factor": ("elements", "factor"),
    "characters": ("elements", "characters"),
    "matrices": ("elements", "matrices"),
    "all": TABLE_ORDER,
}
# the subcommands that accept --double-group (docs/spec/cli.md, "Options")
DOUBLE_GROUP_SUBCOMMANDS = ("elements", "cayley", "all")
# the tables emitted in the index space of G-dagger under --double-group
DOUBLE_GROUP_TABLES = ("elements", "cayley")


def tool_version() -> str:
    try:
        return metadata.version(GENERATOR_NAME)
    except (
        metadata.PackageNotFoundError
    ):  # pragma: no cover - a source tree without metadata
        return "0+unknown"


def parse_values(text: str) -> tuple[str, ...]:
    """``--values``: a comma-separated subset of exact, numeric, both, pretty."""
    requested: set[str] = set()
    for item in text.split(","):
        item = item.strip()
        if item == "both":
            requested.update(("exact", "numeric"))
        elif item in VALUE_FORMS:
            requested.add(item)
        else:
            raise UsageError(
                f"invalid --values item '{item}': expected exact, numeric, both or pretty"
            )
    return tuple(form for form in VALUE_FORMS if form in requested)


@dataclass(frozen=True)
class Request:
    """One invocation of a table subcommand."""

    subcommand: str
    group: str
    values: tuple[str, ...] = DEFAULT_VALUES
    label_convention: str | None = None
    basis_convention: str = conventions.CANONICAL
    double: bool = False
    irreps: str = "all"


@dataclass(frozen=True)
class Model:
    """The computed tables of one group: the source of every document.

    The tables of ``G`` are built eagerly; those of the double group
    ``G-dagger`` (its Cayley table, labels and classes) and the factor table
    are built on first use.
    """

    spec: GroupSpec
    group: PointGroup
    labels: list[Label]
    table: list[list[int]]
    classes: list[ConjugacyClass]
    conductor: int
    inverses: list[int] = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "inverses", inverses(self.table))

    @property
    def order(self) -> int:
        return self.group.order

    # ---- the double group ----------------------------------------------------

    @cached_property
    def double(self) -> DoubleGroup:
        return double_group(self.group)

    @cached_property
    def double_table(self) -> list[list[int]]:
        return self.double.cayley()

    @cached_property
    def double_inverses(self) -> list[int]:
        return inverses(self.double_table)

    @cached_property
    def double_labels(self) -> list[Label]:
        """The labels of ``G-dagger``: those of ``G``, then the same with the bar."""
        return self.labels + [replace(label, barred=True) for label in self.labels]

    @cached_property
    def double_classes(self) -> list[DoubleClass]:
        return double_conjugacy_classes(
            self.double_table, self.double_labels, self.classes
        )

    @cached_property
    def factor(self) -> list[list[int]]:
        return self.double.factor_system()

    # ---- metadata and records --------------------------------------------------

    def group_metadata(self) -> dict[str, Any]:
        meta: dict[str, Any] = {
            "symbol": self.spec.symbol,
            "input_symbol": self.spec.input_symbol,
        }
        meta["family"] = self.spec.family
        if self.spec.n is not None:
            meta["n"] = self.spec.n
        meta["order"] = self.order
        meta["type"] = self.group.type
        meta["rotation_group"] = self.group.rotation_group.name
        return meta

    def element_records(self, double: bool = False) -> list[dict[str, Any]]:
        if double:
            class_of = {
                g: cls.index for cls in self.double_classes for g in cls.members
            }
            labels, inverse = self.double_labels, self.double_inverses
        else:
            class_of = {g: cls.index for cls in self.classes for g in cls.members}
            labels, inverse = self.labels, self.inverses
        N = self.order
        records = []
        for index, label in enumerate(labels):
            parity = self.group.elements[index % N][1]
            record: dict[str, Any] = {
                "index": index,
                "label": label.text,
                "parity": parity,
            }
            record["kind"] = label.kind
            record["order"] = label.order
            if label.power is not None:
                record["power"] = label.power
            if label.axis_class is not None:
                record["axis_class"] = label.axis_class
            if label.axis_index is not None:
                record["axis_index"] = label.axis_index
            record["inverse"] = inverse[index]
            record["class"] = class_of[index]
            if double:
                record["barred"] = index >= N
                record["base_index"] = index % N
            records.append(record)
        return records

    def class_records(self, double: bool = False) -> list[dict[str, Any]]:
        if double:
            return [self._double_class_record(cls) for cls in self.double_classes]
        return [
            {
                "index": cls.index,
                "label": cls.label,
                "size": cls.size,
                "representative": cls.representative,
                "members": list(cls.members),
                "inverse": cls.inverse,
                "regular": cls.regular,
            }
            for cls in self.classes
        ]

    @staticmethod
    def _double_class_record(double_class: DoubleClass) -> dict[str, Any]:
        record: dict[str, Any] = {
            "index": double_class.index,
            "label": double_class.label,
            "size": double_class.size,
            "representative": double_class.representative,
            "members": list(double_class.members),
            "inverse": double_class.inverse,
            "base_class": double_class.base_class,
            "split": double_class.split,
        }
        if double_class.split:
            record["barred"] = double_class.barred
        return record

    # ---- tables --------------------------------------------------------------

    def elements_table(self, double: bool = False) -> dict[str, Any]:
        return {
            "group": "G_double" if double else "G",
            "items": self.element_records(double),
            "classes": self.class_records(double),
        }

    def cayley_table(self, double: bool = False) -> dict[str, Any]:
        table = self.double_table if double else self.table
        return {
            "group": "G_double" if double else "G",
            "table": [list(row) for row in table],
        }

    def factor_table(self) -> dict[str, Any]:
        """``omega(g, h)`` over ``G``, whatever the index space of the other tables."""
        return {"group": "G", "table": [list(row) for row in self.factor]}

    @staticmethod
    def available_tables() -> tuple[str, ...]:
        """The tables a model can produce, in table order."""
        return ("elements", "cayley", "factor")

    def tables(self, double: bool = False) -> dict[str, dict[str, Any]]:
        """Every table this model can produce, keyed by table name."""
        return {
            "elements": self.elements_table(double),
            "cayley": self.cayley_table(double),
            "factor": self.factor_table(),
        }


def build_model(spec: GroupSpec) -> Model:
    G = point_group(spec)
    labels = canonical_labels(G)
    table = G.cayley()
    classes = conjugacy_classes(G, table, labels)
    return Model(spec, G, labels, table, classes, canonical_conductor(G.rotation_group))


def _metadata(
    model: Model, request: Request, tables: tuple[str, ...]
) -> dict[str, Any]:
    label_convention = None
    if request.label_convention is not None:
        label_convention = conventions.label_convention(request.label_convention)
    basis_convention = conventions.basis_convention(request.basis_convention)
    doc: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generator": {"name": GENERATOR_NAME, "version": tool_version()},
        "group": model.group_metadata(),
        "tables": list(tables),
        "double": request.double,
    }
    if "characters" in tables or "matrices" in tables:
        doc["irreps"] = request.irreps
    doc["values"] = list(request.values)
    doc["field"] = field_metadata(model.conductor)
    if label_convention is not None:
        doc["label_convention"] = label_convention
    doc["basis_convention"] = basis_convention
    return doc


def _assemble(request: Request, tables: tuple[str, ...]) -> dict[str, Any]:
    """The metadata block followed by the named tables, in serialisation order."""
    model = build_model(parse_group(request.group))
    doc = _metadata(model, request, tables)
    available = model.tables(request.double)
    for name in tables:
        if name not in available:
            raise NotImplementedError(f"the {name} table is not available yet")
        doc[name] = available[name]
    return doc


def build_document(request: Request) -> dict[str, Any]:
    """The document of one invocation, as a dict in serialisation order."""
    if request.subcommand not in TABLES_OF_SUBCOMMAND:
        raise UsageError(f"unknown subcommand '{request.subcommand}'")
    if request.double and request.subcommand not in DOUBLE_GROUP_SUBCOMMANDS:
        raise UsageError(
            f"--double-group is not accepted by '{request.subcommand}'",
            f"it applies to {', '.join(DOUBLE_GROUP_SUBCOMMANDS)} only",
        )
    return _assemble(request, TABLES_OF_SUBCOMMAND[request.subcommand])


def model_document(group: str, double: bool = False) -> dict[str, Any]:
    """The document holding every table the model can produce: the input of the checkers.

    With ``double`` the element list and the Cayley table are those of
    ``G-dagger``; the factor table is always that of ``G``.
    """
    available = tuple(name for name in TABLE_ORDER if name in Model.available_tables())
    return _assemble(Request("elements", group, double=double), available)
