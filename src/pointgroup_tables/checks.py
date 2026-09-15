"""Self-consistency checks as pure functions of the document model
(docs/spec/verification.md, check ledger).

Every checker takes the document dict (the tables as they will be
serialised) and returns the failures it found, each with the ledger
identifier and one message. The same functions serve the test suite, the
mutation tests and ``--verify``.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Any

from .seed import SEED

ELEMENT_LABEL = re.compile(
    r"^~?(E|i|C2|C2\(\d+\)|C2''?\(\d+\)|C\d+\^-?\d+(\(\d+\))?|S\d+\^-?\d+(\(\d+\))?"
    r"|sigma_h(\(\d+\))?|sigma_[vd]\(\d+\)|sigma\(\d+\))$"
)
CLASS_LABEL = re.compile(
    r"^(\d*~?(E|i|C2''?|[CS]\d+(\^-?\d+)?|sigma(_[hvd])?))"
    r"(\+\d*~(E|i|C2''?|[CS]\d+(\^-?\d+)?|sigma(_[hvd])?))?$"
)

# the class structures of the polyhedral groups (docs/spec/label-conventions.md,
# docs/spec/conjugacy-classes.md)
POLYHEDRAL_CLASSES = {
    "T": ["E", "3C2", "4C3", "4C3^-1"],
    "Td": ["E", "3C2", "8C3", "6sigma_d", "6S4"],
    "Th": ["E", "3C2", "4C3", "4C3^-1", "i", "3sigma_h", "4S6^-1", "4S6"],
    "O": ["E", "3C2", "6C2'", "8C3", "6C4"],
    "Oh": ["E", "3C2", "6C2'", "8C3", "6C4", "i", "3sigma_h", "6sigma_d", "8S6", "6S4"],
    "I": ["E", "15C2", "20C3", "12C5", "12C5^2"],
    "Ih": [
        "E",
        "15C2",
        "20C3",
        "12C5",
        "12C5^2",
        "i",
        "15sigma",
        "20S6",
        "12S10^3",
        "12S10",
    ],
}


@dataclass(frozen=True)
class Failure:
    check: str
    message: str

    def __str__(self) -> str:
        return f"{self.check}: {self.message}"


Checker = Callable[[dict[str, Any]], list[Failure]]
Document = dict[str, Any]


def _cayley(doc: Document) -> list[list[int]]:
    return doc["cayley"]["table"]


def _items(doc: Document) -> list[dict[str, Any]]:
    return doc["elements"]["items"]


def _classes(doc: Document) -> list[dict[str, Any]]:
    return doc["elements"]["classes"]


def _element_orders(table: list[list[int]]) -> list[int]:
    """The order of every element; ``0`` when its powers never return to the identity
    (a corrupted table)."""
    size = len(table)
    orders = []
    for g in range(size):
        n, h = 1, g
        while h != 0 and n <= size:
            h = table[h][g]
            n += 1
        orders.append(n if h == 0 else 0)
    return orders


# ---- AG: groups ------------------------------------------------------------


def check_ag1(doc: Document) -> list[Failure]:
    """Cayley table of ``G``: closure, identity, inverses, associativity (all triples)."""
    table = _cayley(doc)
    n = len(table)
    failures = []
    rows = range(n)
    if any(len(row) != n for row in table):
        return [Failure("AG-1", "the Cayley table is not square")]
    if any(not isinstance(x, int) or x < 0 or x >= n for row in table for x in row):
        return [Failure("AG-1", "an entry of the Cayley table is not an element index")]
    if any(table[0][g] != g or table[g][0] != g for g in rows):
        failures.append(Failure("AG-1", "element 0 is not the identity"))
    for g in rows:
        if 0 not in table[g]:
            failures.append(Failure("AG-1", f"element {g} has no inverse"))
        else:
            inverse = table[g].index(0)
            if table[inverse][g] != 0:
                failures.append(
                    Failure("AG-1", f"the right inverse of {g} is not a left inverse")
                )
    for a in rows:
        row_a = table[a]
        for b in rows:
            row_ab = table[row_a[b]]
            row_b = table[b]
            for c in rows:
                if row_ab[c] != row_a[row_b[c]]:
                    failures.append(Failure("AG-1", f"({a} {b}) {c} != {a} ({b} {c})"))
                    return failures
    return failures


def check_ag2(doc: Document) -> list[Failure]:
    """``(order, class count, histogram of element orders)`` equals the seed table (seed groups only)."""
    symbol = doc["group"]["symbol"]
    row = SEED.get(symbol)
    if row is None:
        return []
    failures = []
    order = len(_items(doc))
    if order != row.order:
        failures.append(
            Failure(
                "AG-2",
                f"order {order} of {symbol} differs from the seed table ({row.order})",
            )
        )
    classes = len(_classes(doc))
    if classes != row.class_count:
        failures.append(
            Failure(
                "AG-2",
                f"{classes} classes of {symbol} differ from the seed table ({row.class_count})",
            )
        )
    histogram = dict(Counter(_element_orders(_cayley(doc))))
    if histogram != row.element_orders:
        failures.append(
            Failure(
                "AG-2",
                f"element orders {histogram} of {symbol} differ from the seed table",
            )
        )
    return failures


def check_ag3(doc: Document) -> list[Failure]:
    """Canonical element order: proper then improper, classes ordered by first member,
    classes of a polyhedral group contiguous."""
    items = _items(doc)
    classes = _classes(doc)
    failures = []
    parities = [record["parity"] for record in items]
    if parities != sorted(parities, reverse=True):
        failures.append(Failure("AG-3", "an improper element precedes a proper one"))
    firsts = [min(cls["members"]) for cls in classes]
    if firsts != sorted(firsts):
        failures.append(
            Failure("AG-3", "classes are not ordered by their first member")
        )
    for cls in classes:
        if cls["representative"] != min(cls["members"]):
            failures.append(
                Failure(
                    "AG-3",
                    f"class {cls['index']} is not represented by its first member",
                )
            )
    if doc["group"]["rotation_group"] in ("T", "O", "I"):
        for cls in classes:
            members = sorted(cls["members"])
            if members != list(range(members[0], members[-1] + 1)):
                failures.append(
                    Failure("AG-3", f"class {cls['index']} is not contiguous")
                )
    return failures


# ---- CL: conjugacy classes -------------------------------------------------


def check_cl1(doc: Document) -> list[Failure]:
    """Classes partition the elements; sizes divide ``|G|``; ``size * |centraliser| = |G|``."""
    items = _items(doc)
    classes = _classes(doc)
    table = _cayley(doc)
    order = len(items)
    failures = []
    seen: Counter[int] = Counter()
    for cls in classes:
        seen.update(cls["members"])
        if cls["size"] != len(cls["members"]):
            failures.append(
                Failure("CL-1", f"class {cls['index']} declares size {cls['size']}")
            )
        if order % len(cls["members"]):
            failures.append(
                Failure(
                    "CL-1", f"the size of class {cls['index']} does not divide {order}"
                )
            )
        if cls["representative"] not in cls["members"]:
            failures.append(
                Failure(
                    "CL-1",
                    f"the representative of class {cls['index']} is not a member",
                )
            )
        rep = cls["representative"]
        centraliser = sum(1 for h in range(order) if table[h][rep] == table[rep][h])
        if len(cls["members"]) * centraliser != order:
            failures.append(
                Failure(
                    "CL-1",
                    f"class {cls['index']}: size {len(cls['members'])} x centraliser {centraliser} != {order}",
                )
            )
        for g in cls["members"]:
            if 0 <= g < order and items[g]["class"] != cls["index"]:
                failures.append(
                    Failure(
                        "CL-1",
                        f"element {g} names class {items[g]['class']}, not {cls['index']}",
                    )
                )
    if sorted(seen) != list(range(order)) or any(count != 1 for count in seen.values()):
        failures.append(Failure("CL-1", "the classes do not partition the elements"))
    return failures


def check_cl2(doc: Document) -> list[Failure]:
    """Class labels match the grammar and are pairwise distinct."""
    failures = []
    labels = [cls["label"] for cls in _classes(doc)]
    for label in labels:
        if not CLASS_LABEL.match(label):
            failures.append(
                Failure("CL-2", f"class label '{label}' does not match the grammar")
            )
    for label, count in Counter(labels).items():
        if count > 1:
            failures.append(
                Failure("CL-2", f"class label '{label}' appears {count} times")
            )
    return failures


def _binary_rotation_part(record: dict[str, Any]) -> bool:
    """Whether the rotation part of an element is a binary rotation."""
    return record["kind"] in ("rotation", "reflection") and record["order"] == 2


def check_cl3(doc: Document) -> list[Failure]:
    """``regular`` marks exactly the classes that split in the double group.

    Recomputed from the records: the class of ``g`` merges exactly when the
    rotation part of ``g`` is a binary rotation and some element with a
    different binary rotation part commutes with ``g``.
    """
    items = _items(doc)
    table = _cayley(doc)
    order = len(items)
    inversion = next((r["index"] for r in items if r["kind"] == "inversion"), None)
    failures = []
    for cls in _classes(doc):
        rep = cls["representative"]
        record = items[rep]
        merged = False
        if _binary_rotation_part(record):
            same_rotation = {rep} if inversion is None else {rep, table[rep][inversion]}
            merged = any(
                _binary_rotation_part(items[h])
                and h not in same_rotation
                and table[h][rep] == table[rep][h]
                for h in range(order)
            )
        if cls["regular"] == merged:
            failures.append(
                Failure(
                    "CL-3",
                    f"class {cls['index']} ('{cls['label']}') has regular = {cls['regular']}",
                )
            )
    return failures


def check_cl4(doc: Document) -> list[Failure]:
    """``inverse`` is an involution on class indices and names the class of the inverses."""
    items = _items(doc)
    classes = _classes(doc)
    failures = []
    for cls in classes:
        partner = cls["inverse"]
        if not 0 <= partner < len(classes):
            failures.append(
                Failure("CL-4", f"class {cls['index']} has an out-of-range inverse")
            )
            continue
        if classes[partner]["inverse"] != cls["index"]:
            failures.append(
                Failure("CL-4", f"inverse is not an involution at class {cls['index']}")
            )
        inverses = {items[g]["inverse"] for g in cls["members"]}
        if inverses != set(classes[partner]["members"]):
            failures.append(
                Failure(
                    "CL-4",
                    f"class {partner} does not hold the inverses of class {cls['index']}",
                )
            )
    return failures


def check_cl5(doc: Document) -> list[Failure]:
    """Polyhedral class labels reproduce the listed class structures."""
    expected = POLYHEDRAL_CLASSES.get(doc["group"]["symbol"])
    if expected is None:
        return []
    labels = [cls["label"] for cls in _classes(doc)]
    if labels != expected:
        return [Failure("CL-5", f"class labels {labels} differ from {expected}")]
    return []


# ---- registry ----------------------------------------------------------------

CHECKS: dict[str, tuple[Checker, tuple[str, ...]]] = {
    "AG-1": (check_ag1, ("cayley",)),
    "AG-2": (check_ag2, ("elements", "cayley")),
    "AG-3": (check_ag3, ("elements",)),
    "CL-1": (check_cl1, ("elements", "cayley")),
    "CL-2": (check_cl2, ("elements",)),
    "CL-3": (check_cl3, ("elements", "cayley")),
    "CL-4": (check_cl4, ("elements",)),
    "CL-5": (check_cl5, ("elements",)),
}


def applicable_checks(doc: Document) -> list[str]:
    """The checks whose tables the document holds, in ledger order."""
    return [
        name
        for name, (_, needs) in CHECKS.items()
        if all(table in doc for table in needs)
    ]


def check(doc: Document, which: Iterable[str] | None = None) -> list[Failure]:
    """Run the named checks (default: every applicable one) and return every failure."""
    names = list(which) if which is not None else applicable_checks(doc)
    failures: list[Failure] = []
    for name in names:
        checker, needs = CHECKS[name]
        missing = [table for table in needs if table not in doc]
        if missing:
            raise ValueError(f"check {name} needs the {', '.join(missing)} table")
        failures.extend(checker(doc))
    return failures
