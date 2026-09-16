"""The ``altmann-herzig`` label convention: generation and check (wayfinder ticket #14).

    uv run python tools/ah_convention.py check       # rules + data reproduce every transcribed file
    uv run python tools/ah_convention.py generate    # (re)write the shipped convention file

The convention file (``src/pointgroup_tables/data/conventions/labels/altmann-herzig.toml``,
spec: ``docs/spec/label-conventions.md``) is derived data: the secondary-direction
table below is read from the book's parameter tables, the polyhedral overrides
and barred lists come from the transcriptions under
``tests/reference/altmann-herzig/labels/``. ``check`` applies the rules of the
spec to the canonical labels of every transcribed group and compares, element
by element, with the book's labels; it also derives the barred (swapped-lift)
list of every axial group from the pole table alone and compares it with the
list the transcription checker solves from the book's tilde marks.

The two prototypes under prototypes/ supply the groups; nothing here is
production code.
"""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
sys.path.insert(0, str(ROOT / "tools"))
import abstract_group_prototype as ag
import ah_book as ab
import ah_transcription as at

LABELS_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "labels"
CHARACTERS_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "characters"
OUT_FILE = (
    ROOT
    / "src"
    / "pointgroup_tables"
    / "data"
    / "conventions"
    / "labels"
    / "altmann-herzig.toml"
)
SOURCE = "S. L. Altmann and P. Herzig, Point-Group Theory Tables, 2nd ed. (Wien: Clarendon, 2011)"

# The book's numbering of the secondary directions of the master group D_Nh:
# the poles of C2'1, C2'2, ... and of C2''1, C2''2, ... as integers in units of
# 180/N degrees (0 <= pole < 2N).  N = 2..10 from the parameter tables T 31.1,
# T 33.1, T 35.1, T 37.1, T 39.1 (transcribed in tools/ah_book.py); N = 12..20
# from the tables of D6d (T 45.1), D7h (T 36.1), D8d (T 47.1), D9h (T 38.1),
# D10d (T 49.1), whose C2' and sigma_v / sigma_d rows are the C2' and C2''
# directions of D_Nh.  The book tabulates no group whose master is beyond D20h.
SECONDARY: dict[int, tuple[list[int], list[int]]] = {
    2: ([0], [1]),
    4: ([0, 2], [1, 3]),
    6: ([0, 4, 8], [3, 7, 11]),
    8: ([0, 4, 2, 6], [1, 5, 3, 7]),
    10: ([0, 4, 8, 12, 16], [5, 9, 13, 17, 1]),
    12: ([0, 8, 16, 6, 14, 22], [3, 11, 19, 9, 17, 1]),
    14: ([0, 4, 8, 12, 16, 20, 24], [7, 11, 15, 19, 23, 27, 31]),
    16: ([0, 8, 4, 12, 2, 10, 6, 14], [1, 9, 5, 13, 3, 11, 7, 15]),
    18: ([0, 12, 24, 4, 16, 28, 8, 20, 32], [9, 21, 33, 13, 25, 1, 17, 29, 5]),
    20: (
        [0, 8, 16, 24, 32, 10, 18, 26, 34, 2],
        [5, 13, 21, 29, 37, 15, 23, 31, 39, 7],
    ),
}

# Cartesian names of the D2h family (the book's subgroup table T 31.0).
OVERRIDES: dict[str, dict[str, str]] = {
    "D2": {"C2": "C2z", "C2'(1)": "C2x", "C2''(1)": "C2y"},
    "D2h": {
        "C2": "C2z",
        "sigma_h": "sigma_z",
        "C2'(1)": "C2x",
        "C2''(1)": "C2y",
        "sigma_d(1)": "sigma_x",
        "sigma_v(1)": "sigma_y",
    },
    "C2v": {"sigma_v(1)": "sigma_y", "sigma_d(1)": "sigma_x"},
}

POLYHEDRAL = ("T", "Th", "Td", "O", "Oh", "I", "Ih")
BOOK_NAME = {"C1h": "Cs", "S2": "Ci"}  # file stem -> name in tools/ah_book.py


def plusminus(label: str) -> str:
    """C3^1 -> C3+, C3^-1 -> C3-, C8^3 -> C8^3+, S10^-3 -> S10^3-; other labels unchanged."""
    m = re.fullmatch(r"([CS])(\d+)\^(-?)(\d+)", label)
    if not m:
        return label
    kind, order, sign, power = m.groups()
    sgn = "-" if sign else "+"
    return f"{kind}{order}{sgn}" if power == "1" else f"{kind}{order}^{power}{sgn}"


class Axial:
    """Reference labels and reference section of an axial group under the convention's rules."""

    def __init__(self, symbol: str):
        self.grp = at.Group(symbol)
        self.G = self.grp.G
        self.fam, self.n = ag.parse(self.G.symbol)
        R = self.G.rotation_group
        self.N = int(R.name[1:])
        self.dihedral = R.name.startswith("D")
        self.M = self.N if self.N % 2 == 0 else 2 * self.N  # the master D_Mh
        # Cnv with n = 2 (mod 4): canonical direction 0 is the master's C2''1 direction
        self.twist = self.fam == "Cv" and self.n % 4 == 2

    @property
    def covered(self) -> bool:
        return not self.dihedral or self.M in SECONDARY

    def master_direction(self, k: int) -> int:
        if self.twist:
            k = (k + self.N // 2) % self.N
        return k * (self.M // self.N)

    def lookup(self, k: int) -> tuple[str, int, int]:
        """(primes, book index, book pole) of the master's secondary axis in direction k."""
        primed, double = SECONDARY[self.M]
        u = self.master_direction(k)
        for primes, poles in (("'", primed), ("''", double)):
            for i, pole in enumerate(poles):
                if pole % self.M == u:
                    return primes, i + 1, pole
        raise AssertionError((self.G.symbol, k, u))

    def reference_label(self, g) -> str:
        canon = self.grp.labels[g]
        override = OVERRIDES.get(self.G.symbol, {})
        if canon in override:
            return override[canon]
        r, p = g
        if not self.dihedral or r[1] == 0 or not self.covered:
            return plusminus(canon)
        primes, index, _ = self.lookup(r[0])
        if p == +1:
            return f"C2{primes}{index}"
        return canon.split("(")[0] + str(
            index
        )  # sigma_v / sigma_d: the canonical letter

    def swapped(self) -> list[str]:
        """Canonical labels of the elements whose book lift is our barred lift."""
        out = []
        if not self.dihedral or not self.covered:
            return out
        for g in self.G.elements:
            r, _ = g
            if r[1] == 0:
                continue
            k = r[0]
            ours = 180.0 * (k if 2 * k <= self.N else k - self.N) / self.N
            if self.twist:
                ours += 90.0
            _, _, pole = self.lookup(k)
            if abs(((180.0 * pole / self.M - ours) % 360.0) - 180.0) < 1e-6:
                out.append(self.grp.labels[g])
        return sorted(out)


def transcribed(
    name: str, data: dict, grp: at.Group
) -> tuple[dict[str, str], list[str]]:
    """canonical -> book label from the file's [map], and the swapped list solved from the tilde marks."""
    book, _ = ab.book_group(name)
    iso, anchors = ab.find_map(name, book, grp)
    for b, c in data["map"].items():
        assert grp.elements[iso[book.index[b]]] == c, (name, b, c)
    eps = ab.solve_epsilon(book, grp, iso, list(anchors))
    assert eps is not None, name
    swapped = sorted(grp.elements[iso[b]] for b, v in eps.items() if v == -1)
    return {c: b for b, c in data["map"].items()}, swapped


def check() -> tuple[dict[str, dict[str, str]], dict[str, list[str]], int]:
    failures = 0
    overrides: dict[str, dict[str, str]] = {}
    barred: dict[str, list[str]] = {}
    for path in sorted(LABELS_DIR.glob("*.toml")):
        data = tomllib.load(path.open("rb"))
        name = BOOK_NAME.get(path.stem, path.stem)
        symbol = data["group"]
        grp = at.Group(symbol)
        want, want_swapped = transcribed(name, data, grp)
        if symbol in POLYHEDRAL:
            overrides[symbol] = want
            barred[symbol] = want_swapped
            print(
                f"{path.stem:5s} polyhedral: {len(want)} override entries, {len(want_swapped)} swapped (data)"
            )
            continue
        ax = Axial(symbol)
        got = {ax.grp.labels[g]: ax.reference_label(g) for g in ax.G.elements}
        bad = {c: (got[c], want[c]) for c in want if got[c] != want[c]}
        got_swapped = ax.swapped()
        ok = not bad and got_swapped == want_swapped
        failures += not ok
        print(
            f"{path.stem:5s} {'OK ' if ok else 'BAD'} master D{ax.M}h"
            + (" twist" if ax.twist else "")
            + f" swapped={got_swapped}"
        )
        if bad:
            print("      label mismatches (canonical: got, book):", bad)
        if got_swapped != want_swapped:
            print(
                "      swapped mismatch, the transcription checker says:", want_swapped
            )
    print("FAILURES:", failures)
    return overrides, barred, failures


# The book's spelling of a few class labels in the character transcription differs
# from the labels transcription (whose README fixes sigma_x) and from the reference
# class label of the convention (15sigma of Ih, docs/spec/conjugacy-classes.md).
CLASS_SPELLING = {
    "sigmax": "sigma_x",
    "sigmay": "sigma_y",
    "sigmaz": "sigma_z",
    "15sigma_h": "15sigma",
}
IRREP_RULES = {
    "component_style": "prefix",
    "fraction_style": "plain",
    "projective_parity_separator": ",",
    "letters": {"G": "F"},
    "projective_letters": {"G": "F", "B": "A"},
}


def irrep_overrides() -> dict[str, dict[str, str]]:
    """``[irreps.overrides]``: the book's symbol of every irrep whose spelled canonical
    symbol differs from it, decided by the characters on the labelled classes of
    the character transcriptions (docs/spec/irrep-labels.md: "any disagreement
    becomes an [irreps.overrides] entry"). The element and class labels of the
    shipped convention identify the columns; its irrep rules play no part here."""
    import math

    from pointgroup_tables.document import build_model, model_document
    from pointgroup_tables.field import cyclotomic_field
    from pointgroup_tables.reference import reference_symbol
    from pointgroup_tables.symbols import parse_group
    from pointgroup_tables.values import parse

    out: dict[str, dict[str, str]] = {}
    for path in sorted(CHARACTERS_DIR.glob("*.toml")):
        data = tomllib.load(path.open("rb"))
        symbol = data["group"]
        document = model_document(symbol, label_convention="altmann-herzig")
        table = document["characters"]
        F_ours = cyclotomic_field(document["field"]["conductor"])
        F_book = cyclotomic_field(data["conductor"])
        F = cyclotomic_field(math.lcm(F_ours.conductor, F_book.conductor))
        labels = [c["label"] for c in data["classes"]]
        rows = [
            (
                row["label"],
                [F_book.lift(parse(F_book, t), F) for t in row["characters"]],
                row["projective"],
            )
            for row in data["irreps"]
        ]
        for entry in data.get("known_discrepancy", []):
            a, b = (labels.index(x) for x in entry["swap_projective"])
            rows = [
                (label, _swapped(values, a, b) if projective else values, projective)
                for label, values, projective in rows
            ]
        columns = []
        for label in labels:
            label = CLASS_SPELLING.get(label, label)
            matches = [
                c
                for c, column in enumerate(table["classes"])
                if column["reference_label"] == label
                or column["reference_label"].startswith(label + "+")
            ]
            assert len(matches) == 1, (symbol, label)
            columns.append(matches[0])
        irreps = build_model(parse_group(symbol)).irreps
        ours = {
            (
                tuple(
                    F_ours.lift(F_ours.from_exact(table["exact"][k][c]), F)
                    for c in columns
                ),
                record["projective"],
            ): reference_symbol(irreps[k], IRREP_RULES)
            for k, record in enumerate(table["irreps"])
        }
        canonical = {
            reference_symbol(irrep, IRREP_RULES): irrep.symbol.text for irrep in irreps
        }
        for label, values, projective in rows:
            spelled = ours[(tuple(values), projective)]
            if spelled != label:
                out.setdefault(symbol, {})[canonical[spelled]] = label
    return out


def _swapped(values: list, a: int, b: int) -> list:
    values = list(values)
    values[a], values[b] = values[b], values[a]
    return values


def toml_list(items) -> str:
    return (
        "["
        + ", ".join(ab.toml_str(x) if isinstance(x, str) else str(x) for x in items)
        + "]"
    )


def emit(
    overrides: dict[str, dict[str, str]],
    barred: dict[str, list[str]],
    irreps: dict[str, dict[str, str]],
) -> str:
    lines = [
        "# The altmann-herzig label convention (docs/spec/label-conventions.md).",
        "# Generated by tools/ah_convention.py from the book's parameter tables and the",
        "# transcriptions under tests/reference/altmann-herzig/labels/; do not edit by hand.",
        "",
        'name = "altmann-herzig"',
        f"source = {ab.toml_str(SOURCE)}",
        'exponent_style = "plusminus"',
        'barred_format = "~{label}"',
        "",
        "# The book's range (docs/spec/output-schema.md, coverage): every axial family up to",
        "# n = 10 and all seven polyhedral groups are verified against the print; beyond,",
        "# the rules extrapolate.",
        "[coverage]",
        "axial_max_n = 10",
        "polyhedral = true",
        "",
        "# The book's numbering of the secondary directions of the master group D_Nh:",
        "# poles of C2'1, C2'2, ... and of C2''1, C2''2, ... in units of 180/N degrees.",
        "[secondary]",
    ]
    for N, (primed, double) in SECONDARY.items():
        lines.append(f'"{N}" = {{ primed = {primed}, double_primed = {double} }}')
    lines.append("")
    for symbol, table in OVERRIDES.items():
        lines.append(f"[overrides.{symbol}]")
        lines += [f"{ab.toml_str(c)} = {ab.toml_str(b)}" for c, b in table.items()]
        lines.append("")
    for symbol in POLYHEDRAL:
        grp = at.Group(symbol)
        lines.append(f"[overrides.{symbol}]")
        lines += [
            f"{ab.toml_str(c)} = {ab.toml_str(overrides[symbol][c])}"
            for c in grp.elements
        ]
        lines.append("")
    lines.append(
        "# Elements whose book lift is our barred lift; axial groups derive theirs from [secondary]."
    )
    for symbol in POLYHEDRAL:
        if barred[symbol]:
            lines += [
                f"[barred.{symbol}]",
                f"swapped = {toml_list(barred[symbol])}",
                "",
            ]
    lines += [
        "[irreps]",
        'component_style = "prefix"',
        'fraction_style = "plain"',
        'projective_parity_separator = ","',
        'letters = { G = "F" }',
        'projective_letters = { G = "F", B = "A" }',
        "",
        "# The book's symbol where it differs from the spelled canonical one: the members",
        "# of the complex pairs, the E_k of D4d and D5h, and the F of T and Th, decided",
        "# by the characters of the transcribed tables (tools/ah_convention.py).",
    ]
    for symbol, table in irreps.items():
        lines.append(f"[irreps.overrides.{symbol}]")
        lines += [f"{ab.toml_str(c)} = {ab.toml_str(b)}" for c, b in table.items()]
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    command = argv[0] if argv else "check"
    overrides, barred, failures = check()
    if command == "generate":
        if failures:
            print("not generated: the rules do not reproduce the transcriptions")
            return 1
        OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
        OUT_FILE.write_text(emit(overrides, barred, irrep_overrides()))
        print("wrote", OUT_FILE.relative_to(ROOT))
    elif command != "check":
        print(__doc__)
        return 2
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
