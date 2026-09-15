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
