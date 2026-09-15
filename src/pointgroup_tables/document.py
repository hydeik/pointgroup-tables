"""The document model: everything one invocation computes, and its projection
to the document of one subcommand (docs/spec/output-schema.md).

``build_document(request)`` is the first seam of the tool: it names the
subcommand, the group as typed and every option, and returns the dict exactly
as serialised. ``model_document`` returns the same dict with every table the
model can produce, the input of the checkers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib import metadata
from typing import Any

from . import conventions
from .abstract import PointGroup, point_group
from .classes import ConjugacyClass, conjugacy_classes, inverses
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
    """The computed tables of one group: the source of every document."""

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

    def element_records(self) -> list[dict[str, Any]]:
        class_of = {g: cls.index for cls in self.classes for g in cls.members}
        records = []
        for index, (label, (_, parity)) in enumerate(
            zip(self.labels, self.group.elements, strict=True)
        ):
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
            record["inverse"] = self.inverses[index]
            record["class"] = class_of[index]
            records.append(record)
        return records

    def class_records(self) -> list[dict[str, Any]]:
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

    def elements_table(self) -> dict[str, Any]:
        return {
            "group": "G",
            "items": self.element_records(),
            "classes": self.class_records(),
        }

    def cayley_table(self) -> dict[str, Any]:
        return {"group": "G", "table": [list(row) for row in self.table]}

    @staticmethod
    def available_tables() -> tuple[str, ...]:
        """The tables a model can produce, in table order."""
        return ("elements", "cayley")

    def tables(self) -> dict[str, dict[str, Any]]:
        """Every table this model can produce, keyed by table name."""
        return {"elements": self.elements_table(), "cayley": self.cayley_table()}


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
    available = model.tables()
    for name in tables:
        if name not in available:
            raise NotImplementedError(f"the {name} table is not available yet")
        doc[name] = available[name]
    return doc


def build_document(request: Request) -> dict[str, Any]:
    """The document of one invocation, as a dict in serialisation order."""
    if request.subcommand not in TABLES_OF_SUBCOMMAND:
        raise UsageError(f"unknown subcommand '{request.subcommand}'")
    return _assemble(request, TABLES_OF_SUBCOMMAND[request.subcommand])


def model_document(group: str) -> dict[str, Any]:
    """The document holding every table the model can produce: the input of the checkers."""
    available = tuple(name for name in TABLE_ORDER if name in Model.available_tables())
    return _assemble(Request("elements", group), available)
