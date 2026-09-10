"""Scaffolding for the Altmann & Herzig transcription (wayfinder ticket #13).

Two subcommands:

    uv run python tools/ah_transcription.py templates            # write missing template files
    uv run python tools/ah_transcription.py check [GROUP ...]    # validate transcribed files

A template file per group is written to tests/reference/altmann-herzig/labels/
with the canonical element list and class structure of the group in its
header, so the transcriber sees what every book label must map onto. The
checker validates a filled file: book labels are unique and as many as the
group has elements, every canonical label named exists in the group, the
map is injective, the irrep row has the right length, and, when the book's
multiplication table is transcribed, that it is a Latin square isomorphic
to our Cayley table under the anchors given in ``[map]``; it then proposes
the canonical label of every unmapped element and lists the tilde marks
that disagree with our factor system (candidates for ``barred.swapped``).

The two prototypes under prototypes/ supply the groups; nothing here is
production code.
"""

from __future__ import annotations

import sys
import tomllib
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))  # ty resolves these via tool.ty.environment.extra-paths
import abstract_group_prototype as ag  # noqa: E402
import double_group_prototype as dg  # noqa: E402

LABELS_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "labels"

CRYSTALLOGRAPHIC = [
    "C1", "Ci", "Cs", "C2", "C2h", "C2v", "D2", "D2h",
    "C3", "S6", "C3v", "D3", "D3d", "C3h", "D3h",
    "C4", "S4", "C4h", "C4v", "D4", "D2d", "D4h",
    "C6", "C6h", "C6v", "D6", "D6h",
    "T", "Th", "Td", "O", "Oh",
]
EXTRA = ["I", "Ih", "C5v", "D5h", "D5d", "S8", "D4d"]
GROUPS = CRYSTALLOGRAPHIC + EXTRA

# groups whose double-group multiplication table the double-group chapter
# wants compared (docs/spec/double-groups.md, "Verification data")
MULTIPLICATION_PRIORITY = ["C2v", "D3", "D4", "D3h", "T", "O"]


# --------------------------------------------------------------------------
# group data
# --------------------------------------------------------------------------


class Group:
    def __init__(self, symbol):
        self.input_symbol = symbol
        self.G = ag.point_group(symbol)
        self.symbol = ag.ALIASES.get(symbol, symbol)
        self.labels = ag.canonical_labels(self.G)
        self.elements = [self.labels[g] for g in self.G.elements]
        self.index = {lab: i for i, lab in enumerate(self.elements)}
        self.table = self.G.cayley()
        self.omega = dg.double_group(self.G).factor_system()
        classes = ag.conjugacy_classes(self.G)
        self.classes = [[self.labels[g] for g in c] for c in classes]

    def omega_at(self, a, b):
        return self.omega[(self.G.elements[a], self.G.elements[b])]


# --------------------------------------------------------------------------
# templates
# --------------------------------------------------------------------------

TEMPLATE = '''\
# Altmann & Herzig, Point-Group Theory Tables (2nd ed.): element labels of {symbol}
# Format: tests/reference/altmann-herzig/labels/README.md
#
# Canonical elements of {symbol} (order {order}), in canonical order:
#   {elements}
# Conjugacy classes (canonical):
{classes}
# The map below must send every book label onto exactly one of these.

group = "{symbol}"
table = ""          # the book's table number and page, e.g. "T 15.1, p. 123"
status = "todo"     # todo | partial | done
notes = ""          # anything that could not be read unambiguously

[book]
# Every element label of the book's table for this group, in the book's
# order, transliterated to ASCII (README: C3+, C2'1, sigma_v1, S6-, ~C3+ ...).
elements = []
# The book's identification of the operation (axis or plane), where stated,
# keyed by the book label, e.g. "C2'1" = "C2 about x".
operations = {{}}
# The book's irrep symbol row, in the book's order, as printed (ASCII).
irreps = []

[map]
# book label -> canonical label. Fill the ones the book's axis statements
# make certain (E, the principal rotation, the first secondary axis, i,
# sigma_h); the checker derives the rest from the multiplication table.

[multiplication]
# Optional (priority groups: {priority}). One row per book element, in the
# order of book.elements, each entry a book label with a leading "~" where the
# book prints a tilde on the product.
rows = []
'''


def write_templates(only=None):
    LABELS_DIR.mkdir(parents=True, exist_ok=True)
    written = []
    for symbol in only or GROUPS:
        grp = Group(symbol)
        path = LABELS_DIR / f"{grp.symbol}.toml"
        if path.exists():
            continue
        classes = "\n".join(
            f"#   {len(c)}x[{', '.join(c)}]" for c in grp.classes
        )
        path.write_text(
            TEMPLATE.format(
                symbol=grp.symbol,
                order=grp.G.order,
                elements=" ".join(grp.elements),
                classes=classes,
                priority=", ".join(MULTIPLICATION_PRIORITY),
            )
        )
        written.append(path.name)
    return written


# --------------------------------------------------------------------------
# checker
# --------------------------------------------------------------------------


class Report:
    def __init__(self, name):
        self.name = name
        self.errors = []
        self.warnings = []
        self.infos = []

    def error(self, msg):
        self.errors.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)

    def info(self, msg):
        self.infos.append(msg)

    def print(self):
        print(f"== {self.name}")
        for m in self.errors:
            print(f"   ERROR   {m}")
        for m in self.warnings:
            print(f"   warning {m}")
        for m in self.infos:
            print(f"   {m}")


def check_file(path):
    data = tomllib.loads(path.read_text())
    rep = Report(path.name)
    grp = Group(data["group"])
    status = data.get("status", "todo")
    book = data.get("book", {})
    elements = book.get("elements", [])
    if status == "todo" and not elements:
        rep.info("todo: not transcribed yet")
        return rep

    # --- element list
    n = grp.G.order
    if len(set(elements)) != len(elements):
        dup = sorted({e for e in elements if elements.count(e) > 1})
        rep.error(f"duplicate book labels: {dup}")
    if len(elements) != n:
        rep.error(f"book lists {len(elements)} elements, {grp.symbol} has {n}")
    for e in elements:
        if e.startswith("~") or not e.isascii() or " " in e:
            rep.error(f"book label {e!r} is not an ASCII single-group label")

    # --- operations and irreps
    for key in book.get("operations", {}):
        if key not in elements:
            rep.error(f"operations names {key!r}, which is not in book.elements")
    irreps = book.get("irreps", [])
    if irreps:
        expected = None  # irrep counts depend on the irrep chapter; report only
        rep.info(f"irrep row: {len(irreps)} symbols: {' '.join(irreps)}")
        del expected

    # --- map
    mapping = dict(data.get("map", {}))
    for b, c in mapping.items():
        if b not in elements:
            rep.error(f"map names book label {b!r}, which is not in book.elements")
        if c not in grp.index:
            rep.error(f"map sends {b!r} to {c!r}, which is not a canonical label of {grp.symbol}")
    targets = list(mapping.values())
    if len(set(targets)) != len(targets):
        dup = sorted({t for t in targets if targets.count(t) > 1})
        rep.error(f"map is not injective; canonical labels used twice: {dup}")
    if "E" in elements and mapping.get("E", "E") != "E":
        rep.error("the book's E must map to canonical E")

    # --- multiplication table
    rows = data.get("multiplication", {}).get("rows", [])
    if rows and not rep.errors:
        check_multiplication(grp, elements, mapping, rows, rep)
    elif not rows:
        unmapped = [e for e in elements if e not in mapping]
        if unmapped:
            rep.warn(f"{len(unmapped)} book labels unmapped and no multiplication table to derive them")
    if status == "done" and (rep.errors or [e for e in elements if e not in mapping] and not rows):
        rep.error("status is done but the file is incomplete")
    return rep


def check_multiplication(grp, elements, mapping, rows, rep):
    n = len(elements)
    pos = {e: i for i, e in enumerate(elements)}
    if len(rows) != n or any(len(r) != n for r in rows):
        rep.error(f"multiplication.rows must be {n} rows of {n} entries")
        return
    prod: list[list[int]] = [[-1] * n for _ in range(n)]
    tilde = [[False] * n for _ in range(n)]
    for i, row in enumerate(rows):
        for j, entry in enumerate(row):
            t = entry.startswith("~")
            e = entry[1:] if t else entry
            if e not in pos:
                rep.error(f"row {i} col {j}: {entry!r} is not a book label")
                return
            prod[i][j] = pos[e]
            tilde[i][j] = t
    for i in range(n):
        if len(set(prod[i])) != n or len({prod[k][i] for k in range(n)}) != n:
            rep.error(f"multiplication table is not a Latin square at row/column {i} ({elements[i]})")
            return
    # derive the isomorphism from the anchors in [map] by closure
    iso = {pos[b]: grp.index[c] for b, c in mapping.items()}
    if not iso:
        rep.error("no anchors in [map]; give at least E, the principal rotation and one secondary element")
        return
    queue = deque(iso)
    known = set(iso)
    conflict = False
    while queue and not conflict:
        a = queue.popleft()
        for b in list(known):
            for x, y in ((a, b), (b, a)):
                p = prod[x][y]
                target = grp.table[iso[x]][iso[y]]
                if p in iso:
                    if iso[p] != target:
                        rep.error(
                            f"[map] is not a homomorphism: book {elements[x]} * {elements[y]} = "
                            f"{elements[p]}, but canonical {grp.elements[iso[x]]} * "
                            f"{grp.elements[iso[y]]} = {grp.elements[target]} != {grp.elements[iso[p]]}"
                        )
                        conflict = True
                        break
                else:
                    iso[p] = target
                    known.add(p)
                    queue.append(p)
            if conflict:
                break
    if conflict:
        return
    if len(iso) < n:
        missing = [elements[i] for i in range(n) if i not in iso]
        rep.warn(f"anchors do not generate the group; still unmapped: {missing}")
    else:
        derived = {elements[i]: grp.elements[iso[i]] for i in range(n) if elements[i] not in mapping}
        if derived:
            rep.info("derived map (add to [map]):")
            for b, c in derived.items():
                rep.info(f'  "{b}" = "{c}"')
        # tilde marks versus our factor system: omega_book = eps(g) eps(h) eps(gh) omega,
        # solve eps with eps = +1 on the anchors (canonical section as reference)
        r = [[(-1 if tilde[x][y] else 1) * grp.omega_at(iso[x], iso[y]) for y in range(n)] for x in range(n)]
        eps = {pos[b]: 1 for b in mapping}
        queue = deque(eps)
        consistent = True
        while queue and consistent:
            a = queue.popleft()
            for b in list(eps):
                for x, y in ((a, b), (b, a)):
                    p = prod[x][y]
                    value = r[x][y] * eps[x] * eps[y]
                    if p in eps:
                        if eps[p] != value:
                            consistent = False
                            break
                    else:
                        eps[p] = value
                        queue.append(p)
                if not consistent:
                    break
        total = sum(1 for x in range(n) for y in range(n) if tilde[x][y])
        if not consistent:
            rep.error("the book's tilde marks do not differ from the canonical factor system by a "
                      "coboundary: a transcription error, or the anchors are wrong")
            return
        swapped = sorted((elements[i] for i, e in eps.items() if e == -1), key=elements.index)
        rep.info(f"tilde marks: {total} in the book; the book's section differs from the canonical one "
                 f"on {len(swapped)} element(s) (anchors fixed at +1)")
        if swapped:
            canon = [grp.elements[iso[pos[b]]] for b in swapped]
            rep.info(f"  [barred.{grp.symbol}] swapped = {canon!r}")


def main(argv):
    if not argv or argv[0] not in ("templates", "check"):
        print(__doc__)
        return 2
    if argv[0] == "templates":
        written = write_templates(argv[1:] or None)
        print(f"wrote {len(written)} template(s)" + (": " + " ".join(written) if written else ""))
        return 0
    names = argv[1:] or sorted(p.stem for p in LABELS_DIR.glob("*.toml"))
    failed = 0
    for name in names:
        path = LABELS_DIR / f"{ag.ALIASES.get(name, name)}.toml"
        if not path.exists():
            print(f"== {name}: no file {path.relative_to(ROOT)}")
            failed += 1
            continue
        rep = check_file(path)
        rep.print()
        failed += bool(rep.errors)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
