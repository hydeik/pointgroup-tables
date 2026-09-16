"""Self-consistency checks as pure functions of the document model
(docs/spec/verification.md, check ledger).

Every checker takes the document dict (the tables as they will be
serialised) and returns the failures it found, each with the ledger
identifier and one message. The same functions serve the test suite, the
mutation tests and ``--verify``.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import pairwise
from typing import Any

from .field import (
    CyclotomicField,
    Matrix,
    Value,
    cyclotomic_field,
    matrix_identity,
    matrix_mul,
    matrix_trace,
)
from .irreps import (
    POLYHEDRAL_PRINCIPAL_ORDER,
    SYMBOL,
    NamedElements,
    Symbol,
    derive_symbol,
)
from .seed import SEED
from .values import ValueSyntaxError, numeric, parse

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


def _factor(doc: Document) -> list[list[int]]:
    return doc["factor"]["table"]


def _is_double(doc: Document) -> bool:
    """Whether the element list and the Cayley table are those of ``G-dagger``."""
    return doc["elements"]["group"] == "G_double"


def _base_order(doc: Document) -> int:
    """``|G|``, whatever the index space of the element list."""
    items = _items(doc)
    return len(items) // 2 if _is_double(doc) else len(items)


def _base_cayley(doc: Document) -> list[list[int]]:
    """The Cayley table of ``G``; from a ``G-dagger`` table, the products modulo ``|G|``."""
    table = _cayley(doc)
    if not _is_double(doc):
        return table
    N = len(table) // 2
    return [[x % N for x in row[:N]] for row in table[:N]]


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


def _group_axioms(table: list[list[int]], check: str) -> list[Failure]:
    """Closure, identity, inverses and associativity (all triples) of a Cayley table."""
    n = len(table)
    failures = []
    rows = range(n)
    if any(len(row) != n for row in table):
        return [Failure(check, "the Cayley table is not square")]
    if any(not isinstance(x, int) or x < 0 or x >= n for row in table for x in row):
        return [Failure(check, "an entry of the Cayley table is not an element index")]
    if any(table[0][g] != g or table[g][0] != g for g in rows):
        failures.append(Failure(check, "element 0 is not the identity"))
    for g in rows:
        if 0 not in table[g]:
            failures.append(Failure(check, f"element {g} has no inverse"))
        else:
            inverse = table[g].index(0)
            if table[inverse][g] != 0:
                failures.append(
                    Failure(check, f"the right inverse of {g} is not a left inverse")
                )
    for a in rows:
        row_a = table[a]
        for b in rows:
            # (ab)c == a(bc) for every c: the row of ab against row_a permuted by row_b
            if table[row_a[b]] != [row_a[x] for x in table[b]]:
                failures.append(
                    Failure(check, f"({a} {b}) c != {a} ({b} c) for some c")
                )
                return failures
    return failures


def check_ag1(doc: Document) -> list[Failure]:
    """Cayley table of ``G`` (or of ``G-dagger``): closure, identity, inverses, associativity."""
    return _group_axioms(_cayley(doc), "AG-1")


def check_ag2(doc: Document) -> list[Failure]:
    """``(order, class count, histogram of element orders)`` equals the seed table (seed groups only)."""
    symbol = doc["group"]["symbol"]
    row = SEED.get(symbol)
    if row is None or _is_double(doc):
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
    classes of a polyhedral group contiguous. A property of ``G``; the order of
    ``G-dagger`` is checked by ``DG-4``."""
    if _is_double(doc):
        return []
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

    On ``G`` it is recomputed from the records: the class of ``g`` merges
    exactly when the rotation part of ``g`` is a binary rotation and some
    element with a different binary rotation part commutes with ``g``. On
    ``G-dagger`` the ``split`` and ``barred`` flags are recomputed from the
    class memberships: a class splits when it does not contain the barred
    element of its representative.
    """
    if _is_double(doc):
        return _check_cl3_double(doc)
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


def _check_cl3_double(doc: Document) -> list[Failure]:
    N = _base_order(doc)
    failures = []
    for cls in _classes(doc):
        rep = cls["representative"]
        members = set(cls["members"])
        split = (rep + N) % (2 * N) not in members
        if cls["split"] is not split:
            failures.append(
                Failure(
                    "CL-3",
                    f"class {cls['index']} ('{cls['label']}') has split = {cls['split']}",
                )
            )
        elif split and cls.get("barred") is not (rep >= N):
            failures.append(
                Failure(
                    "CL-3",
                    f"class {cls['index']} ('{cls['label']}') has barred = {cls.get('barred')}",
                )
            )
    return failures


def check_cl4(doc: Document) -> list[Failure]:
    """``inverse`` is an involution on class indices and names the class of the inverses;
    on ``G-dagger``, ``base_class``, ``split`` and ``barred`` agree with the element records."""
    items = _items(doc)
    classes = _classes(doc)
    failures = []
    if _is_double(doc):
        failures.extend(_check_cl4_double(doc))
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


def _check_cl4_double(doc: Document) -> list[Failure]:
    """The ``G-dagger`` class records against the ``barred`` / ``base_index`` element fields."""
    items = _items(doc)
    classes = _classes(doc)
    failures = []
    previous_base = -1
    covered: dict[int, set[int]] = {}
    for cls in classes:
        index = cls["index"]
        base = cls["base_class"]
        if base < previous_base:
            failures.append(
                Failure("CL-4", f"class {index} is not ordered by base class")
            )
        previous_base = base
        bases = [items[g]["base_index"] for g in cls["members"]]
        if cls["split"]:
            if "barred" not in cls:
                failures.append(
                    Failure("CL-4", f"split class {index} has no barred flag")
                )
            elif items[cls["representative"]]["barred"] is not cls["barred"]:
                failures.append(
                    Failure(
                        "CL-4",
                        f"the representative of class {index} contradicts barred",
                    )
                )
            if len(set(bases)) != len(bases):
                failures.append(
                    Failure("CL-4", f"split class {index} covers a base element twice")
                )
        else:
            if "barred" in cls:
                failures.append(
                    Failure("CL-4", f"merged class {index} has a barred flag")
                )
            if sorted(Counter(bases).values()) != [2] * (len(bases) // 2):
                failures.append(
                    Failure(
                        "CL-4",
                        f"merged class {index} does not cover every base element twice",
                    )
                )
        covered.setdefault(base, set()).update(bases)
        if len(covered[base]) != len(set(bases)):
            failures.append(
                Failure(
                    "CL-4",
                    f"class {index} does not cover the same base elements as its pair",
                )
            )
    for base in covered:
        pair = [cls for cls in classes if cls["base_class"] == base]
        if len(pair) > 2 or (len(pair) == 2 and not all(cls["split"] for cls in pair)):
            failures.append(
                Failure("CL-4", f"base class {base} is covered by {len(pair)} classes")
            )
        if len(pair) == 2 and [cls.get("barred") for cls in pair] != [False, True]:
            failures.append(
                Failure(
                    "CL-4",
                    f"the classes over base class {base} are not lifted then barred",
                )
            )
    return failures


def check_cl5(doc: Document) -> list[Failure]:
    """Polyhedral class labels reproduce the listed class structures (``G`` only)."""
    expected = POLYHEDRAL_CLASSES.get(doc["group"]["symbol"])
    if expected is None or _is_double(doc):
        return []
    labels = [cls["label"] for cls in _classes(doc)]
    if labels != expected:
        return [Failure("CL-5", f"class labels {labels} differ from {expected}")]
    return []


# ---- DG: double group and factor system ------------------------------------


def check_dg1(doc: Document) -> list[Failure]:
    """``omega(E, g) = omega(g, E) = 1`` and every value is ``1`` or ``-1``."""
    omega = _factor(doc)
    n = len(omega)
    if n != _base_order(doc) or any(len(row) != n for row in omega):
        return [Failure("DG-1", "the factor table is not |G| x |G|")]
    failures = []
    if any(x not in (1, -1) or isinstance(x, bool) for row in omega for x in row):
        failures.append(Failure("DG-1", "a factor is not 1 or -1"))
    if any(omega[0][g] != 1 or omega[g][0] != 1 for g in range(n)):
        failures.append(Failure("DG-1", "the factor system is not normalised"))
    return failures


def check_dg2(doc: Document) -> list[Failure]:
    """The cocycle condition ``omega(g, h) omega(gh, k) = omega(g, hk) omega(h, k)``, all triples."""
    omega = _factor(doc)
    table = _base_cayley(doc)
    n = len(table)
    if len(omega) != n or any(len(row) != n for row in omega):
        return [Failure("DG-2", "the factor table does not match the Cayley table")]
    for g in range(n):
        omega_g = omega[g]
        for h in range(n):
            row_gh = omega[table[g][h]]
            left = row_gh if omega_g[h] == 1 else [-x for x in row_gh]
            right = [omega_g[hk] * w for hk, w in zip(table[h], omega[h], strict=True)]
            if left != right:
                k = next(k for k in range(n) if left[k] != right[k])
                return [
                    Failure("DG-2", f"the cocycle condition fails at ({g}, {h}, {k})")
                ]
    return []


def check_dg3(doc: Document) -> list[Failure]:
    """``omega(g, g^-1) = omega(g^-1, g)``."""
    omega = _factor(doc)
    table = _base_cayley(doc)
    failures = []
    for g, row in enumerate(table):
        if 0 not in row:
            continue  # AG-1 reports the missing inverse
        inverse = row.index(0)
        if omega[g][inverse] != omega[inverse][g]:
            failures.append(
                Failure("DG-3", f"omega({g}, {inverse}) != omega({inverse}, {g})")
            )
    return failures


def check_dg4(doc: Document) -> list[Failure]:
    """The Cayley table of ``G-dagger`` satisfies the group axioms; ``|G-dagger| = 2|G|``;
    ``~E`` is central of order 2 and bars every element; the records carry the
    ``barred`` / ``base_index`` structure of the double group. ``G-dagger`` documents only."""
    if not _is_double(doc):
        return []
    items = _items(doc)
    table = _cayley(doc)
    if len(items) % 2 or len(table) != len(items):
        return [
            Failure("DG-4", "the double group does not have an even number of elements")
        ]
    N = len(items) // 2
    failures = _group_axioms(table, "DG-4")
    if failures:
        return failures
    if table[N][N] != 0:
        failures.append(Failure("DG-4", "~E does not have order 2"))
    for g in range(2 * N):
        if table[N][g] != table[g][N]:
            failures.append(Failure("DG-4", f"~E does not commute with element {g}"))
        if table[N][g] != (g + N) % (2 * N):
            failures.append(
                Failure("DG-4", f"~E times element {g} is not its barred partner")
            )
    for record in items:
        i = record["index"]
        base = items[i % N]
        if record.get("base_index") != i % N or record.get("barred") is not (i >= N):
            failures.append(
                Failure("DG-4", f"element {i} does not cover element {i % N}")
            )
        if i >= N and record["label"] != "~" + base["label"]:
            failures.append(
                Failure("DG-4", f"element {i} is not labelled ~{base['label']}")
            )
        if record["parity"] != base["parity"]:
            failures.append(
                Failure("DG-4", f"element {i} changes the parity of {i % N}")
            )
    for g in range(N):
        for h in range(N):
            product = table[g][h] % N
            if table[g + N][h] % N != product or table[g][h + N] % N != product:
                failures.append(
                    Failure(
                        "DG-4", f"the projection is not a homomorphism at ({g}, {h})"
                    )
                )
                return failures
    return failures


def check_dg5(doc: Document) -> list[Failure]:
    """``omega`` depends on rotation parts only (for a centro group it is the pull-back of
    ``omega`` of ``G+``); ``omega(C2, C2) = -1`` for every binary rotation part."""
    items = _items(doc)
    omega = _factor(doc)
    table = _base_cayley(doc)
    N = len(table)
    failures = []
    inversion = next((r["index"] for r in items[:N] if r["kind"] == "inversion"), None)
    if inversion is not None:
        for g in range(N):
            gi = table[g][inversion]
            if omega[g] != omega[gi] or [row[g] for row in omega] != [
                row[gi] for row in omega
            ]:
                failures.append(
                    Failure(
                        "DG-5", f"omega distinguishes element {g} from {gi} = {g} i"
                    )
                )
    for record in items[:N]:
        g = record["index"]
        if _binary_rotation_part(record) and omega[g][g] != -1:
            failures.append(
                Failure("DG-5", f"omega({record['label']}, {record['label']}) != -1")
            )
        if record["kind"] in ("identity", "inversion") and omega[g][g] != 1:
            failures.append(
                Failure("DG-5", f"omega({record['label']}, {record['label']}) != 1")
            )
    return failures


def check_dg6(doc: Document) -> list[Failure]:
    """Class counts of ``G-dagger`` and the number of projective irreps equal the seed table
    (seed groups only): on ``G`` the regular classes count the projective irreps, on
    ``G-dagger`` the classes and the split pairs are counted."""
    symbol = doc["group"]["symbol"]
    row = SEED.get(symbol)
    if row is None:
        return []
    classes = _classes(doc)
    failures = []
    if _is_double(doc):
        if len(classes) != row.double_class_count:
            failures.append(
                Failure(
                    "DG-6",
                    f"{len(classes)} classes of {symbol}-dagger differ from the seed table ({row.double_class_count})",
                )
            )
        split = sum(1 for cls in classes if cls["split"])
        if split != 2 * row.projective_irreps:
            failures.append(
                Failure(
                    "DG-6",
                    f"{split} split classes of {symbol}-dagger, the seed table has {row.projective_irreps} projective irreps",
                )
            )
    else:
        regular = sum(1 for cls in classes if cls["regular"])
        if regular != row.projective_irreps:
            failures.append(
                Failure(
                    "DG-6",
                    f"{regular} regular classes of {symbol}, the seed table has {row.projective_irreps} projective irreps",
                )
            )
    return failures


# ---- the character table: shared readers -------------------------------------


def _field(doc: Document) -> CyclotomicField:
    return cyclotomic_field(doc["field"]["conductor"])


def _character_rows(doc: Document) -> list[list[Value]] | None:
    """The exact character table as field values, or ``None`` without an exact form."""
    table = doc["characters"]
    if "exact" not in table:
        return None
    F = _field(doc)
    try:
        return [[F.from_exact(pairs) for pairs in row] for row in table["exact"]]
    except (ValueError, TypeError):
        return None  # VL-1 reports the malformed value


def _no_exact(check: str) -> list[Failure]:
    return [
        Failure(
            check, "the characters table carries no well-formed exact form to check"
        )
    ]


def _columns(doc: Document) -> list[dict[str, Any]]:
    return doc["characters"]["classes"]


def _irrep_records(doc: Document) -> list[dict[str, Any]]:
    return doc["characters"]["irreps"]


def _column_group_order(doc: Document) -> int:
    """The order of the group the columns belong to (``G`` or ``G-dagger``)."""
    order = doc["group"]["order"]
    return 2 * order if doc["characters"]["class_group"] == "G_double" else order


def _base_size(cls: dict[str, Any]) -> int:
    """The size of the class of ``G`` a column covers."""
    if "split" in cls and not cls["split"]:
        return cls["size"] // 2
    return cls["size"]


def _blocks_present(doc: Document) -> tuple[bool, bool]:
    """Whether the complete ordinary block and the complete projective block are present."""
    which = doc.get("irreps", "all")
    return which in ("all", "ordinary"), which in ("all", "projective")


def _double_cayley(doc: Document) -> list[list[int]]:
    """The Cayley table of ``G-dagger`` in its index space: from the document's own
    table when it is that of ``G-dagger``, else from the table of ``G`` and the
    factor system through ``s(g) s(h) = omega(g, h) s(gh)``."""
    if _is_double(doc):
        return _cayley(doc)
    table, omega = _cayley(doc), _factor(doc)
    N = len(table)
    out = []
    for x in range(2 * N):
        g, bar_x = x % N, x >= N
        row = []
        for y in range(2 * N):
            h, bar_y = y % N, y >= N
            barred = bar_x != bar_y
            if omega[g][h] == -1:
                barred = not barred
            row.append(table[g][h] + (N if barred else 0))
        out.append(row)
    return out


# ---- IL: irrep labels ----------------------------------------------------------


def _record_symbol(record: dict[str, Any]) -> Symbol:
    """The symbol the structured fields of a record spell."""
    return Symbol(
        record["letter"],
        record.get("subscript"),
        record.get("suffix"),
        record.get("component"),
    )


def check_il1(doc: Document) -> list[Failure]:
    """Exactly one canonical symbol per irrep, grammar-conformant, pairwise distinct,
    and spelled by the record's letter, subscript, suffix and component."""
    failures = []
    symbols = []
    for record in _irrep_records(doc):
        symbol = record["symbol"]
        symbols.append(symbol)
        if not SYMBOL.match(symbol):
            failures.append(
                Failure("IL-1", f"symbol '{symbol}' does not match the grammar")
            )
        elif _record_symbol(record).text != symbol:
            failures.append(
                Failure(
                    "IL-1",
                    f"the fields of irrep {record['index']} do not spell '{symbol}'",
                )
            )
    for symbol, count in Counter(symbols).items():
        if count > 1:
            failures.append(Failure("IL-1", f"symbol '{symbol}' appears {count} times"))
    return failures


def _turn(record: dict[str, Any]) -> Fraction | None:
    """The turn of the rotation part of a principal-axis element, or ``None``."""
    kind = record["kind"]
    if kind in ("identity", "inversion"):
        return Fraction(0)
    if record.get("axis_class") not in (None, "principal"):
        return None
    if kind == "reflection":
        return Fraction(1, 2)  # sigma_h = i C2
    if record["order"] == 2 and kind == "rotation":
        return Fraction(1, 2)
    turn = Fraction(record["power"], record["order"])
    return turn if kind == "rotation" else (turn + Fraction(1, 2)) % 1


def _named_from_document(doc: Document) -> NamedElements | None:
    """The named elements of the irrep-labels chapter, located in the element records."""
    group = doc["group"]
    rotation_group = group["rotation_group"]
    items = [record for record in _items(doc) if not record.get("barred")]
    inversion = next((r["index"] for r in items if r["kind"] == "inversion"), None)
    if rotation_group in POLYHEDRAL_PRINCIPAL_ORDER:
        # the polyhedral principal element: the first element of the principal order
        order = POLYHEDRAL_PRINCIPAL_ORDER[rotation_group]
        principal = next((r["index"] for r in items if r["order"] == order), None)
        if principal is None:
            return None
        return NamedElements(
            principal=principal,
            proper_principal=principal,
            order=order,
            proper_order=order,
            secondary=None,
            double_prime=None,
            inversion=inversion,
            horizontal=None,
            d2=False,
            polyhedral=rotation_group,
        )
    N = int(rotation_group[1:])
    family, n = group["family"], group["n"]
    odd_h = family in ("Cnh", "Dnh") and n % 2 == 1
    d2 = family in ("Dn", "Dnh") and n == 2

    def principal_with_turn(turn: Fraction) -> int | None:
        candidates = [record["index"] for record in items if _turn(record) == turn % 1]
        proper = [g for g in candidates if items[g]["parity"] == 1]
        chosen = proper or candidates
        return chosen[0] if chosen else None

    def labelled(*labels: str) -> int | None:
        return next((r["index"] for r in items if r["label"] in labels), None)

    principal = principal_with_turn(Fraction(1, N))
    proper_principal = principal_with_turn(Fraction(2, N)) if odd_h else principal
    if principal is None or proper_principal is None:
        return None
    return NamedElements(
        principal=principal,
        proper_principal=proper_principal,
        order=N,
        proper_order=n if odd_h else N,
        secondary=labelled("C2'(1)", "sigma_v(1)"),
        double_prime=labelled("C2''(1)") if d2 else None,
        inversion=inversion,
        horizontal=labelled("sigma_h") if odd_h else None,
        d2=d2,
    )


def _row_on_elements(doc: Document, row: list[Value]) -> list[Value]:
    """A character row spread over the elements of ``G`` (the lifted elements)."""
    N = doc["group"]["order"]
    column_of: dict[int, int] = {}
    for column, cls in enumerate(_columns(doc)):
        for member in cls["members"]:
            if member < N:
                column_of[member] = column
    return [row[column_of[g]] for g in range(N)]


def check_il2(doc: Document) -> list[Failure]:
    """The symbol is a function of the character: recomputed from the character table
    by the rules of the chapter, it is reproduced."""
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("IL-2")
    named = _named_from_document(doc)
    if named is None:
        return [
            Failure("IL-2", "the named elements are missing from the element records")
        ]
    F = _field(doc)
    failures = []
    for record, row in zip(_irrep_records(doc), rows, strict=True):
        try:
            symbol = derive_symbol(
                F, named, _row_on_elements(doc, row), record["projective"]
            )
        except (ValueError, KeyError) as exc:
            failures.append(
                Failure(
                    "IL-2",
                    f"no symbol follows from the character of irrep {record['index']}: {exc}",
                )
            )
            continue
        if symbol.text != record["symbol"]:
            failures.append(
                Failure(
                    "IL-2",
                    f"irrep {record['index']} is named '{record['symbol']}' but its character gives '{symbol.text}'",
                )
            )
    return failures


def check_il3(doc: Document) -> list[Failure]:
    """``dimension = chi(E)``; ``frobenius_schur = (1/|G|) sum omega(g, g) chi(g^2)``;
    ``partner`` is present exactly when the indicator is ``0``, and names the
    conjugate irrep when that row is present."""
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("IL-3")
    F = _field(doc)
    table, omega = _base_cayley(doc), _factor(doc)
    N = len(table)
    records = _irrep_records(doc)
    failures = []
    conjugates = {tuple(F.conjugate(x) for x in row): k for k, row in enumerate(rows)}
    for k, (record, row) in enumerate(zip(records, rows, strict=True)):
        chi = _row_on_elements(doc, row)
        if chi[0] != F.rational(record["dimension"]):
            failures.append(
                Failure(
                    "IL-3",
                    f"irrep {record['index']} declares dimension {record['dimension']}",
                )
            )
        total = F.total(
            F.neg(chi[table[g][g]])
            if record["projective"] and omega[g][g] == -1
            else chi[table[g][g]]
            for g in range(N)
        )
        indicator = F.scale(total, Fraction(1, N))
        if indicator != F.rational(record["frobenius_schur"]):
            failures.append(
                Failure(
                    "IL-3",
                    f"irrep {record['index']} declares frobenius_schur {record['frobenius_schur']}",
                )
            )
        has_partner = "partner" in record
        if has_partner != (record["frobenius_schur"] == 0):
            failures.append(
                Failure(
                    "IL-3",
                    f"irrep {record['index']} has a partner iff its indicator is 0",
                )
            )
        if has_partner:
            partner = conjugates.get(tuple(row))
            if record["partner"] == record["index"]:
                failures.append(
                    Failure("IL-3", f"irrep {record['index']} is its own partner")
                )
            elif partner is not None and records[partner]["index"] != record["partner"]:
                failures.append(
                    Failure(
                        "IL-3",
                        f"the partner of irrep {record['index']} is not its complex conjugate",
                    )
                )
    return failures


def check_il4(doc: Document) -> list[Failure]:
    """The rows follow the canonical irrep order and carry increasing global indices."""
    records = _irrep_records(doc)
    failures = []
    keys = [_record_symbol(r).sort_key(r["projective"]) for r in records]
    if keys != sorted(keys):
        failures.append(Failure("IL-4", "the irreps are not in canonical irrep order"))
    indices = [r["index"] for r in records]
    if indices != sorted(set(indices)):
        failures.append(Failure("IL-4", "irrep indices are not strictly increasing"))
    return failures


# ---- CH: characters -----------------------------------------------------------------


def check_ch1(doc: Document) -> list[Failure]:
    """Row orthogonality over the column group: ``sum_C |C| chi_i(C) chi_j(C)* = |group| delta_ij``."""
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("CH-1")
    F = _field(doc)
    sizes = [cls["size"] for cls in _columns(doc)]
    order = _column_group_order(doc)
    records = _irrep_records(doc)
    failures = []
    conjugated = [[F.conjugate(x) for x in row] for row in rows]
    for i, row in enumerate(rows):
        for j in range(i, len(rows)):
            inner = F.total(
                F.scale(F.mul(x, y), size)
                for x, y, size in zip(row, conjugated[j], sizes, strict=True)
            )
            expected = F.rational(order if i == j else 0)
            if inner != expected:
                failures.append(
                    Failure(
                        "CH-1",
                        f"rows {records[i]['index']} and {records[j]['index']} are not orthogonal",
                    )
                )
    return failures


def check_ch2(doc: Document) -> list[Failure]:
    """Column orthogonality: ``sum_i chi_i(C) chi_i(C')* = |centraliser| delta_CC'`` over the
    complete set of irreps of the column group. When only one block of irreps is
    present, its own contribution is checked: the ordinary block contributes
    ``|Z_G(g)| delta_base``, the projective block ``|Z_G-dagger(C)| delta_CC' - |Z_G(g)| delta_base``."""
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("CH-2")
    F = _field(doc)
    columns = _columns(doc)
    order = doc["group"]["order"]
    double = doc["characters"]["class_group"] == "G_double"
    ordinary, projective = _blocks_present(doc)
    failures = []
    conjugated = [[F.conjugate(x) for x in row] for row in rows]
    for c, cls in enumerate(columns):
        for c2 in range(c, len(columns)):
            other = columns[c2]
            same_base = (
                cls.get("base_class", c) == other.get("base_class", c2)
                if double
                else c == c2
            )
            centraliser_g = Fraction(order, _base_size(cls))
            expected = Fraction(0)
            if ordinary and same_base:
                expected += centraliser_g
            if projective:
                if c == c2:
                    expected += Fraction(2 * order if double else order, cls["size"])
                if same_base:
                    expected -= centraliser_g
            if not double and projective and not ordinary:
                expected = Fraction(
                    0
                )  # projective rows over G classes: not a class function
            inner = F.total(
                F.mul(row[c], conj[c2])
                for row, conj in zip(rows, conjugated, strict=True)
            )
            if inner != F.rational(expected):
                failures.append(
                    Failure(
                        "CH-2",
                        f"columns {c} ('{cls['label']}') and {c2} ('{other['label']}') are not orthogonal",
                    )
                )
    return failures


def check_ch3(doc: Document) -> list[Failure]:
    """``sum d^2 = |G|`` over each complete block of irreps; the number of ordinary
    irreps equals the class count of ``G``."""
    records = _irrep_records(doc)
    order = doc["group"]["order"]
    ordinary, projective = _blocks_present(doc)
    failures = []
    for present, flag, name in (
        (ordinary, False, "ordinary"),
        (projective, True, "projective"),
    ):
        if not present:
            continue
        total = sum(r["dimension"] ** 2 for r in records if r["projective"] is flag)
        if total != order:
            failures.append(
                Failure(
                    "CH-3", f"sum of d^2 over the {name} irreps is {total}, not {order}"
                )
            )
    if ordinary:
        columns = _columns(doc)
        if doc["characters"]["class_group"] == "G_double":
            class_count = len({cls["base_class"] for cls in columns})
        else:
            class_count = len(columns)
        count = sum(1 for r in records if not r["projective"])
        if count != class_count:
            failures.append(
                Failure(
                    "CH-3", f"{count} ordinary irreps for {class_count} classes of G"
                )
            )
    return failures


def check_ch4(doc: Document) -> list[Failure]:
    """Characters are class functions of ``G-dagger``: on a split pair the barred class
    repeats an ordinary character and negates a projective one; projective
    characters vanish on merged classes."""
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("CH-4")
    if doc["characters"]["class_group"] != "G_double":
        return []
    F = _field(doc)
    columns = _columns(doc)
    lifted_of_base = {
        cls["base_class"]: c
        for c, cls in enumerate(columns)
        if not cls.get("barred", False)
    }
    failures = []
    for record, row in zip(_irrep_records(doc), rows, strict=True):
        for c, cls in enumerate(columns):
            if not cls["split"]:
                if record["projective"] and not F.is_zero(row[c]):
                    failures.append(
                        Failure(
                            "CH-4",
                            f"irrep {record['index']} is non-zero on the merged class '{cls['label']}'",
                        )
                    )
            elif cls["barred"]:
                lifted = row[lifted_of_base[cls["base_class"]]]
                expected = F.neg(lifted) if record["projective"] else lifted
                if row[c] != expected:
                    failures.append(
                        Failure(
                            "CH-4",
                            f"irrep {record['index']} on '{cls['label']}' is not the signed value of the lifted class",
                        )
                    )
    return failures


def check_ch5(doc: Document) -> list[Failure]:
    """The Frobenius-Schur indicator is ``+-1`` iff the character is real and ``0`` iff
    it is not; over all irreps of ``G-dagger``, ``sum nu_i d_i`` counts the
    elements of ``G-dagger`` that square to the identity."""
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("CH-5")
    F = _field(doc)
    records = _irrep_records(doc)
    failures = []
    for record, row in zip(records, rows, strict=True):
        real = all(F.is_real(x) for x in row)
        if real != (record["frobenius_schur"] != 0):
            failures.append(
                Failure(
                    "CH-5",
                    f"irrep {record['index']} has indicator {record['frobenius_schur']} but is {'real' if real else 'complex'}",
                )
            )
    ordinary, projective = _blocks_present(doc)
    if ordinary and projective and "cayley" in doc and "factor" in doc:
        table, omega = _base_cayley(doc), _factor(doc)
        involutions = 2 * sum(
            1 for g in range(len(table)) if table[g][g] == 0 and omega[g][g] == 1
        )
        total = sum(r["frobenius_schur"] * r["dimension"] for r in records)
        if total != involutions:
            failures.append(
                Failure(
                    "CH-5",
                    f"sum of nu_i d_i is {total}, G-dagger has {involutions} elements squaring to E",
                )
            )
    return failures


def check_ch6(doc: Document) -> list[Failure]:
    """The characters are those of the Burnside class-matrix method: each row gives a
    homomorphism of the class algebra, ``w_i(C_j) w_i(C_k) = sum_l c_jkl w_i(C_l)``
    with ``w_i(C) = |C| chi_i(C) / d_i`` and the class multiplication constants
    ``c_jkl`` read off the Cayley table (of ``G-dagger``, built from the table of
    ``G`` and the factor system), and distinct rows are distinct homomorphisms."""
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("CH-6")
    F = _field(doc)
    columns = _columns(doc)
    double = doc["characters"]["class_group"] == "G_double"
    table = _double_cayley(doc) if double else _base_cayley(doc)
    members = [list(cls["members"]) for cls in columns]
    class_of = {g: c for c, cls in enumerate(members) for g in cls}
    n = len(columns)
    constants: list[list[dict[int, Fraction]]] = [
        [{} for _ in range(n)] for _ in range(n)
    ]
    for j in range(n):
        for k in range(j, n):
            counter: Counter[int] = Counter()
            for x in members[j]:
                row_x = table[x]
                for y in members[k]:
                    counter[class_of[row_x[y]]] += 1
            constants[j][k] = {
                l: Fraction(count, len(members[l])) for l, count in counter.items()
            }
    records = _irrep_records(doc)
    failures = []
    for record, row in zip(records, rows, strict=True):
        d = record["dimension"]
        central = [F.scale(x, Fraction(len(members[c]), d)) for c, x in enumerate(row)]
        for j in range(n):
            for k in range(j, n):
                left = F.mul(central[j], central[k])
                right = F.total(
                    F.scale(central[l], c) for l, c in constants[j][k].items()
                )
                if left != right:
                    failures.append(
                        Failure(
                            "CH-6",
                            f"irrep {record['index']} is not a class-algebra homomorphism at classes {j}, {k}",
                        )
                    )
                    break
            else:
                continue
            break
    if len({tuple(row) for row in rows}) != len(rows):
        failures.append(Failure("CH-6", "two rows carry the same character"))
    return failures


# ---- MX: irrep matrices ------------------------------------------------------------

# the weights of I5/2 of I-dagger on s(C5), the one irrep whose canonical basis is
# ordered by weight, not by reduced exponent (docs/spec/basis-conventions.md,
# "Canonical basis: the principle")
I52_WEIGHTS = (5, 3, 1, -1, -3, -5)


def _matrix_entries(doc: Document) -> list[tuple[dict[str, Any], list[Matrix]]] | None:
    """``(record, D(g) for g in G)`` per irrep of the matrices table as field values,
    or ``None`` when an entry is not a well-formed exact value (``MX-7`` reports it)."""
    F = _field(doc)
    out = []
    try:
        for entry in doc["matrices"]["irreps"]:
            if "exact" not in entry:
                return None
            matrices = [
                [[F.from_exact(pairs) for pairs in row] for row in matrix]
                for matrix in entry["exact"]
            ]
            out.append((entry, matrices))
    except (ValueError, TypeError):
        return None
    return out


def _no_matrices(check: str) -> list[Failure]:
    return [
        Failure(check, "the matrices table carries no well-formed exact form to check")
    ]


def _projective_sign(doc: Document, projective: bool) -> Callable[[int, int], int]:
    """``omega(g, h)`` for a projective irrep, ``1`` for an ordinary one."""
    if not projective:
        return lambda g, h: 1
    omega = _factor(doc)
    return lambda g, h: omega[g][h]


def _generating_set(table: list[list[int]]) -> list[int]:
    """A generating set of the group by greedy choice, elements of the highest
    order first (two or three generators for every point group)."""
    N = len(table)
    orders = _element_orders(table)
    generators: list[int] = []
    reached = {0}
    for g in sorted(range(1, N), key=lambda g: (-orders[g], g)):
        if g in reached:
            continue
        generators.append(g)
        frontier = list(reached)
        while frontier:
            following = []
            for x in frontier:
                for s in generators:
                    y = table[x][s]
                    if y not in reached:
                        reached.add(y)
                        following.append(y)
            frontier = following
    return generators


def _mx_matrix_product_check(
    doc: Document, check: str, pairs: Callable[[int], Iterable[tuple[int, int]]]
) -> list[Failure]:
    """``D(g) D(h) = omega(g, h) D(gh)`` over the pairs a function yields for ``|G|``."""
    entries = _matrix_entries(doc)
    if entries is None:
        return _no_matrices(check)
    F = _field(doc)
    table = _base_cayley(doc)
    N = len(table)
    failures = []
    for record, matrices in entries:
        if len(matrices) != N:
            continue  # MX-7 reports the shape
        sign = _projective_sign(doc, record["projective"])
        for g, h in pairs(N):
            product = matrix_mul(F, matrices[g], matrices[h])
            expected = matrices[table[g][h]]
            if sign(g, h) == -1:
                expected = [[F.neg(x) for x in row] for row in expected]
            if product != expected:
                failures.append(
                    Failure(
                        check,
                        f"irrep {record['index']}: D({g}) D({h}) is not omega({g}, {h}) D({table[g][h]})",
                    )
                )
                break
    return failures


def check_mx1(doc: Document) -> list[Failure]:
    """Homomorphism by closure: ``D(E)`` is the identity and ``D(g) D(s) = omega(g, s)
    D(gs)`` for every element ``g`` and every generator ``s`` of a generating set, so
    that the matrices are the closure of the generator images along the words."""
    entries = _matrix_entries(doc)
    if entries is None:
        return _no_matrices("MX-1")
    F = _field(doc)
    failures = []
    for record, matrices in entries:
        if matrices and matrices[0] != matrix_identity(F, len(matrices[0])):
            failures.append(
                Failure("MX-1", f"irrep {record['index']}: D(E) is not the identity")
            )
    generators = _generating_set(_base_cayley(doc))
    failures.extend(
        _mx_matrix_product_check(
            doc, "MX-1", lambda N: ((g, s) for g in range(N) for s in generators)
        )
    )
    return failures


def check_mx2(doc: Document) -> list[Failure]:
    """``D(g) D(h) = omega(g, h) D(gh)`` for all pairs and all irreps."""
    return _mx_matrix_product_check(
        doc, "MX-2", lambda N: ((g, h) for g in range(N) for h in range(N))
    )


def check_mx3(doc: Document) -> list[Failure]:
    """Every matrix is unitary under the field's conjugation ``zeta -> zeta^-1``."""
    entries = _matrix_entries(doc)
    if entries is None:
        return _no_matrices("MX-3")
    F = _field(doc)
    failures = []
    for record, matrices in entries:
        for g, matrix in enumerate(matrices):
            d = len(matrix)
            adjoint = [[F.conjugate(matrix[j][i]) for j in range(d)] for i in range(d)]
            if matrix_mul(F, adjoint, matrix) != matrix_identity(F, d):
                failures.append(
                    Failure("MX-3", f"irrep {record['index']}: D({g}) is not unitary")
                )
                break
    return failures


def check_mx4(doc: Document) -> list[Failure]:
    """The trace of every matrix equals the character of the irrep of the same index."""
    entries = _matrix_entries(doc)
    if entries is None:
        return _no_matrices("MX-4")
    rows = _character_rows(doc)
    if rows is None:
        return _no_exact("MX-4")
    F = _field(doc)
    by_index = {
        record["index"]: _row_on_elements(doc, row)
        for record, row in zip(_irrep_records(doc), rows, strict=True)
    }
    failures = []
    for record, matrices in entries:
        characters = by_index.get(record["index"])
        if characters is None:
            failures.append(
                Failure(
                    "MX-4",
                    f"irrep {record['index']} has matrices but no character row",
                )
            )
            continue
        for g, (matrix, character) in enumerate(
            zip(matrices, characters, strict=False)
        ):
            if matrix_trace(F, matrix) != character:
                failures.append(
                    Failure(
                        "MX-4",
                        f"irrep {record['index']}: the trace of D({g}) is not the character",
                    )
                )
                break
    return failures


def _reduced_exponent(F: CyclotomicField, x: Value, period: int) -> int | None:
    """``e`` in ``(-period/2, period/2]`` with ``x = zeta_period^e``, or ``None``."""
    for e in range(-((period - 1) // 2), period // 2 + 1):
        if x == F.root_of_unity(period, e):
            return e
    return None


def check_mx5(doc: Document) -> list[Failure]:
    """Canonical basis: ``D(c)`` of the principal element is diagonal, its eigenvalues
    ``zeta_N^e`` (``zeta_2N^e`` for a projective irrep) listed with the reduced
    exponents descending; ``I5/2`` of ``I-dagger`` is ordered by weight instead."""
    if doc["basis_convention"]["name"] != "canonical":
        return []
    entries = _matrix_entries(doc)
    if entries is None:
        return _no_matrices("MX-5")
    named = _named_from_document(doc)
    if named is None:
        return [Failure("MX-5", "the principal element is not among the elements")]
    F = _field(doc)
    c, order = named.principal, named.order
    failures = []
    for record, matrices in entries:
        if c >= len(matrices):
            continue  # MX-7 reports the shape
        matrix = matrices[c]
        d = len(matrix)
        if d < 2:
            continue
        if any(
            not F.is_zero(matrix[i][j]) for i in range(d) for j in range(d) if i != j
        ):
            failures.append(
                Failure("MX-5", f"irrep {record['index']}: D(c) is not diagonal")
            )
            continue
        period = 2 * order if record["projective"] else order
        diagonal = [matrix[i][i] for i in range(d)]
        if record["projective"] and record["letter"] == "I":
            expected = [F.root_of_unity(period, w) for w in I52_WEIGHTS]
            if diagonal != expected:
                failures.append(
                    Failure(
                        "MX-5",
                        f"irrep {record['index']}: the eigenvalues of D(c) are not in weight order",
                    )
                )
            continue
        exponents = [_reduced_exponent(F, x, period) for x in diagonal]
        reduced = [e for e in exponents if e is not None]
        if len(reduced) != len(exponents):
            failures.append(
                Failure(
                    "MX-5",
                    f"irrep {record['index']}: an eigenvalue of D(c) is not a {period}-th root of unity",
                )
            )
        elif any(a <= b for a, b in pairwise(reduced)):
            failures.append(
                Failure(
                    "MX-5",
                    f"irrep {record['index']}: the exponents of D(c) are {reduced}, not descending",
                )
            )
    return failures


def check_mx7(doc: Document) -> list[Failure]:
    """Every irrep carries ``|G|`` matrices of shape ``dimension x dimension`` whose
    entries are well-formed exact values of the declared field."""
    N = doc["group"]["order"]
    degree = doc["field"]["degree"]
    failures = []
    for entry in doc["matrices"]["irreps"]:
        name = f"irrep {entry['index']}"
        d = entry["dimension"]
        for form in ("exact", "numeric", "pretty"):
            if form not in entry:
                continue
            matrices = entry[form]
            if len(matrices) != N:
                failures.append(
                    Failure(
                        "MX-7",
                        f"{name}: the {form} form has {len(matrices)} matrices, not {N}",
                    )
                )
                continue
            if any(
                len(matrix) != d or any(len(row) != d for row in matrix)
                for matrix in matrices
            ):
                failures.append(
                    Failure("MX-7", f"{name}: a {form} matrix is not {d} x {d}")
                )
        if "exact" not in entry:
            continue
        for value in _leaves(entry["exact"], 3):
            if (
                not isinstance(value, list)
                or len(value) != degree
                or any(
                    not isinstance(pair, list)
                    or len(pair) != 2
                    or not isinstance(pair[0], int)
                    or not isinstance(pair[1], int)
                    or pair[1] <= 0
                    or math.gcd(pair[0], pair[1]) != 1
                    for pair in value
                )
            ):
                failures.append(
                    Failure(
                        "MX-7",
                        f"{name}: an entry is not a reduced coefficient vector of degree {degree}",
                    )
                )
                break
    return failures


def check_mx8(doc: Document) -> list[Failure]:
    """Irreducible and pairwise inequivalent: ``(1/|G|) sum_g |chi(g)|^2 = 1`` for the
    trace of every irrep, and ``sum_g chi_i(g) chi_j(g)* = 0`` for two irreps of the
    same kind (an ordinary and a projective irrep are orthogonal over ``G-dagger``
    by the sign of the barred elements alone)."""
    entries = _matrix_entries(doc)
    if entries is None:
        return _no_matrices("MX-8")
    F = _field(doc)
    N = doc["group"]["order"]
    traces = [
        (record, [matrix_trace(F, matrix) for matrix in matrices])
        for record, matrices in entries
        if len(matrices) == N
    ]
    failures = []
    for i, (record, chi) in enumerate(traces):
        conjugated = [F.conjugate(x) for x in chi]
        for j in range(i, len(traces)):
            other, psi = traces[j]
            if i != j and other["projective"] != record["projective"]:
                continue
            inner = F.total(F.mul(x, y) for x, y in zip(psi, conjugated, strict=True))
            if inner != F.rational(N if i == j else 0):
                failures.append(
                    Failure(
                        "MX-8",
                        f"irrep {record['index']} is not irreducible"
                        if i == j
                        else f"irreps {record['index']} and {other['index']} are not inequivalent",
                    )
                )
    return failures


# ---- VL: values -------------------------------------------------------------------------


def _value_tables(doc: Document) -> list[tuple[str, dict[str, Any]]]:
    """The tables that carry values: ``characters`` and every irrep of ``matrices``."""
    out = []
    if "characters" in doc:
        out.append(("characters", doc["characters"]))
    if "matrices" in doc:
        for irrep in doc["matrices"]["irreps"]:
            out.append((f"matrices irrep {irrep['index']}", irrep))
    return out


def _leaves(value: Any, depth: int) -> Iterable[Any]:
    """The entries ``depth`` levels down in a nested array."""
    if depth == 0:
        yield value
        return
    for item in value:
        yield from _leaves(item, depth - 1)


def _value_depth(name: str) -> int:
    """Nesting of one value inside a table: ``[row][column]`` for characters,
    ``[element][row][column]`` for a matrix table."""
    return 2 if name == "characters" else 3


def check_vl1(doc: Document) -> list[Failure]:
    """Every exact value has ``phi(M)`` reduced integer pairs (``q > 0``, ``gcd = 1``)."""
    degree = doc["field"]["degree"]
    failures = []
    for name, table in _value_tables(doc):
        if "exact" not in table:
            continue
        for value in _leaves(table["exact"], _value_depth(name)):
            if len(value) != degree:
                failures.append(
                    Failure(
                        "VL-1",
                        f"an exact value of {name} has {len(value)} coefficients, not {degree}",
                    )
                )
                break
            if any(
                not isinstance(p, int)
                or not isinstance(q, int)
                or q <= 0
                or math.gcd(p, q) != 1
                for p, q in value
            ):
                failures.append(
                    Failure("VL-1", f"an exact value of {name} is not reduced")
                )
                break
    return failures


def check_vl2(doc: Document) -> list[Failure]:
    """Every numeric value is the correctly rounded double of its exact value, and every
    pretty value parses back to it."""
    F = _field(doc)
    failures = []
    for name, table in _value_tables(doc):
        if "exact" not in table:
            continue
        depth = _value_depth(name)
        exact = list(_leaves(table["exact"], depth))
        try:
            [F.from_exact(pairs) for pairs in exact]
        except (ValueError, TypeError):
            continue  # VL-1 reports the malformed value
        numeric_leaves = (
            list(_leaves(table["numeric"], depth)) if "numeric" in table else exact
        )
        pretty_leaves = (
            list(_leaves(table["pretty"], depth)) if "pretty" in table else exact
        )
        if len(numeric_leaves) != len(exact) or len(pretty_leaves) != len(exact):
            failures.append(
                Failure("VL-2", f"the value forms of {name} differ in shape")
            )
            continue
        if "numeric" in table:
            for pairs, approx in zip(exact, numeric_leaves, strict=True):
                expected = numeric(F, F.from_exact(pairs))
                if any(
                    a != b or math.copysign(1.0, a) != math.copysign(1.0, b)
                    for a, b in zip(expected, approx, strict=True)
                ):
                    failures.append(
                        Failure(
                            "VL-2",
                            f"a numeric value of {name} is {approx}, not {expected}",
                        )
                    )
                    break
        if "pretty" in table:
            for pairs, text in zip(exact, pretty_leaves, strict=True):
                try:
                    parsed = parse(F, text)
                except ValueSyntaxError as exc:
                    failures.append(
                        Failure(
                            "VL-2", f"a pretty value of {name} does not parse: {exc}"
                        )
                    )
                    break
                if parsed != F.from_exact(pairs):
                    failures.append(
                        Failure(
                            "VL-2",
                            f"the pretty value '{text}' of {name} is not its exact value",
                        )
                    )
                    break
    return failures


def check_vl3(doc: Document) -> list[Failure]:
    """The declared conductor is minimal: no maximal proper divisor ``M / p`` of ``M``
    contains every value of the document (checked on the matrices, which fix the
    field, and on the characters)."""
    F = _field(doc)
    M = F.conductor
    values: list[Value] = []
    for name, table in _value_tables(doc):
        if "exact" not in table:
            continue
        try:
            values.extend(
                F.from_exact(pairs)
                for pairs in _leaves(table["exact"], _value_depth(name))
            )
        except (ValueError, TypeError):
            continue  # VL-1 and MX-7 report the malformed value
    failures = []
    for p in range(2, M + 1):
        if M % p or any(p % q == 0 for q in range(2, p)):
            continue
        smaller = cyclotomic_field(M // p)
        if all(F.descend(x, smaller) is not None for x in values):
            failures.append(
                Failure(
                    "VL-3",
                    f"every value lies in Q(zeta_{M // p}); the conductor {M} is not minimal",
                )
            )
    return failures


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
    "DG-1": (check_dg1, ("elements", "factor")),
    "DG-2": (check_dg2, ("cayley", "factor")),
    "DG-3": (check_dg3, ("cayley", "factor")),
    "DG-4": (check_dg4, ("elements", "cayley")),
    "DG-5": (check_dg5, ("elements", "cayley", "factor")),
    "DG-6": (check_dg6, ("elements",)),
    "IL-1": (check_il1, ("characters",)),
    "IL-2": (check_il2, ("elements", "characters")),
    "IL-3": (check_il3, ("elements", "cayley", "factor", "characters")),
    "IL-4": (check_il4, ("characters",)),
    "CH-1": (check_ch1, ("characters",)),
    "CH-2": (check_ch2, ("characters",)),
    "CH-3": (check_ch3, ("characters",)),
    "CH-4": (check_ch4, ("characters",)),
    "CH-5": (check_ch5, ("characters",)),
    "CH-6": (check_ch6, ("elements", "cayley", "factor", "characters")),
    "MX-1": (check_mx1, ("cayley", "factor", "matrices")),
    "MX-2": (check_mx2, ("cayley", "factor", "matrices")),
    "MX-3": (check_mx3, ("matrices",)),
    "MX-4": (check_mx4, ("elements", "characters", "matrices")),
    "MX-5": (check_mx5, ("elements", "matrices")),
    "MX-7": (check_mx7, ("matrices",)),
    "MX-8": (check_mx8, ("matrices",)),
    "VL-1": (check_vl1, ("characters",)),
    "VL-2": (check_vl2, ("characters",)),
    "VL-3": (check_vl3, ("matrices",)),
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
