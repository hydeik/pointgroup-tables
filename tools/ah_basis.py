"""Build and verify the ``altmann-herzig`` basis convention (wayfinder ticket #18).

The convention file ``src/pointgroup_tables/data/conventions/basis/altmann-herzig.toml``
is generated from the transcriptions of ticket #17
(``tests/reference/altmann-herzig/matrices/{T,O,I}.toml``), the element-label
maps of ticket #13 (``tests/reference/altmann-herzig/labels/<Group>.toml``)
and the ``swapped`` lists of the label convention (ticket #14):

- each generator the book uses is named by its canonical label; where the
  book's untilded operation is the barred lift of the canonical element
  (``swapped``), the generator is written with the ``~`` prefix so that the
  book's matrix is kept verbatim;
- irreps are keyed by their canonical Mulliken symbol (the ``[canonical]``
  map of the transcription), one-dimensional irreps omitted;
- the conductor is the transcription's (24 for 2T and 2O, 120 for 2I);
- the ``[axial]`` table selects the spin-harmonic rule anchored by the
  ``altmann-herzig`` label convention (ticket #21; the rule itself is
  ``tools/ah_axial_rule.py``, verified there against the axial samples).

    uv run python tools/ah_basis.py check       # verify the convention file against the group
    uv run python tools/ah_basis.py generate    # (re)write the convention file, then verify it

Verification, exact in ``Q(zeta_M)``: the named generators are lifted by the
canonical section (negated for ``~``) and the matrices are closed by
breadth-first products over the binary cover ``R-dagger``; every element
must be reached by consistent products (this is ``D(a) D(b) = omega(a, b)
D(ab)`` for the spin factor system), every generator image must be unitary,
and the trace of every matrix must equal the character of the irrep the key
names (the Burnside characters of the schema prototype).
"""

from __future__ import annotations

import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
sys.path.insert(0, str(ROOT / "tools"))
import canonical_basis_prototype as cb

REF_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "matrices"
LABEL_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "labels"
LABEL_CONVENTION = (
    ROOT
    / "src"
    / "pointgroup_tables"
    / "data"
    / "conventions"
    / "labels"
    / "altmann-herzig.toml"
)
OUT = (
    ROOT
    / "src"
    / "pointgroup_tables"
    / "data"
    / "conventions"
    / "basis"
    / "altmann-herzig.toml"
)

SOURCE = "S. L. Altmann and P. Herzig, Point-Group Theory Tables, 2nd ed. (Wien: Clarendon, 2011)"
COVER = {"T": "2T", "O": "2O", "I": "2I"}


# --------------------------------------------------------------------------
# Building the convention from the reference data
# --------------------------------------------------------------------------


def swapped_labels(group):
    conv = tomllib.loads(LABEL_CONVENTION.read_text())
    return set(conv.get("barred", {}).get(group, {}).get("swapped", []))


def build(group):
    """The ``[polyhedral.<R-dagger>]`` block of the convention, from the references."""
    ref = tomllib.loads((REF_DIR / f"{group}.toml").read_text())
    labels = tomllib.loads((LABEL_DIR / f"{group}.toml").read_text())
    swapped = swapped_labels(group)
    elements = ref["book"]["elements"]
    generators = []
    for book_label in ref["book"]["generators"]:
        canonical = labels["map"][book_label]
        generators.append(("~" if canonical in swapped else "") + canonical)
    irreps = {}
    for book_symbol, canonical_symbol in ref["canonical"].items():
        mats = ref["irreps"][book_symbol]["matrices"]
        irreps[canonical_symbol] = {
            "book_symbol": book_symbol,
            "matrices": [mats[elements.index(g)] for g in ref["book"]["generators"]],
        }
    return {
        "conductor": ref["conductor"],
        "generators": generators,
        "book_generators": ref["book"]["generators"],
        "table": ref["table"],
        "irreps": irreps,
    }


def toml_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def key(symbol):
    return symbol if symbol.isalnum() else toml_str(symbol)


def render(blocks):
    lines = [
        "# Altmann & Herzig basis convention: the book's polyhedral irrep matrices as generator",
        "# images (docs/spec/basis-conventions.md). Generated and verified by tools/ah_basis.py from",
        "# tests/reference/altmann-herzig/matrices/{T,O,I}.toml (ticket #17, including the two",
        "# corrected misprints of T 74.7b recorded there as [[known_discrepancy]]) and the element",
        "# maps of tests/reference/altmann-herzig/labels/ (ticket #13). Decided in ticket #18.",
        "#",
        "# The axial families follow the spin-harmonic rule (docs/spec/basis-conventions.md, ticket",
        "# #21): the book's T n.6 bases <|j m>, +-|j -m>| in closed form, anchored by the altmann-herzig",
        "# label convention (alpha = pi/2 for the twisted Cnv, 0 otherwise). Verified against the",
        "# transcribed axial samples by tools/ah_axial_rule.py; the book covers n <= 10, larger n",
        "# extrapolates the same formulas.",
        "",
        'name = "altmann-herzig"',
        f"source = {toml_str(SOURCE)}",
        "",
        "[axial]",
        'rule = "spin-harmonic"',
        'anchors = "altmann-herzig"',
    ]
    for group, block in blocks.items():
        cover = COVER[group]
        lines.append("")
        lines.append(f"[polyhedral.{cover}]  # {block['table']}")
        lines.append(f"conductor = {block['conductor']}")
        gens = ", ".join(toml_str(g) for g in block["generators"])
        book = ", ".join(block["book_generators"])
        lines.append(
            f"generators = [{gens}]  # the book's {book}; ~ marks a generator whose book lift is the barred one"
        )
        for symbol, entry in block["irreps"].items():
            lines.append("")
            note = (
                ""
                if entry["book_symbol"] == symbol
                else f"  # the book's {entry['book_symbol']}"
            )
            lines.append(f"[polyhedral.{cover}.irreps.{key(symbol)}]{note}")
            lines.append("matrices = [")
            for A in entry["matrices"]:
                lines.append("  [")
                for row in A:
                    lines.append("    [" + ", ".join(toml_str(a) for a in row) + "],")
                lines.append("  ],")
            lines.append("]")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Verifying a convention file against the abstract groups
# --------------------------------------------------------------------------


def lift(P, label):
    """Quaternion of the canonical element named by ``label``, lifted by the canonical section
    (the barred lift for a ``~`` prefix)."""
    barred = label.startswith("~")
    name = label[1:] if barred else label
    g = next(g for g in P.G.elements if P.labels[g] == name)
    q = P.cov.section(g[0])
    return tuple((-a, -b) for a, b in q) if barred else q


def closure(P, F, gens, mats):
    """Breadth-first closure over R-dagger from the generator images; None on an inconsistency."""
    one = next(q for q in P.Q if all(b == 0 for _, b in q) and q[0][0] == 1)
    d = len(mats[0])
    D = {one: cb.midentity(F, d)}
    frontier = [one]
    while frontier:
        nxt = []
        for q in frontier:
            for gq, m in zip(gens, mats):
                p = P.cov.mul(q, gq)
                prod = cb.mmul(F, D[q], m)
                if p in D:
                    if not cb.meq(D[p], prod):
                        return None
                else:
                    D[p] = prod
                    nxt.append(p)
        frontier = nxt
    return D


def verify(path=OUT):
    data = tomllib.loads(Path(path).read_text())
    ok = True
    for kind, cover in COVER.items():
        block = data["polyhedral"][cover]
        M = block["conductor"]
        P = cb.Poly(
            kind, cb.BASE_CONDUCTOR[kind]
        )  # characters live in the canonical field
        F = cb.field(M)
        assert M % P.F.conductor == 0, (M, P.F.conductor)
        gens = [lift(P, label) for label in block["generators"]]
        print(f"{cover}: conductor {M}, generators {block['generators']}")
        for symbol, entry in block["irreps"].items():
            mats = [
                [[cb.parse_value(F, a) for a in row] for row in A]
                for A in entry["matrices"]
            ]
            for A in mats:
                if not cb.is_unitary(F, A):
                    print(f"   FAIL {symbol}: generator image not unitary")
                    ok = False
            D = closure(P, F, gens, mats)
            if D is None or len(D) != len(P.Q):
                print(f"   FAIL {symbol}: closure inconsistent or incomplete")
                ok = False
                continue
            chi, ir = P.character(symbol)
            bad = [q for q in P.Q if cb.mtrace(F, D[q]) != cb.lift(P.F, F, chi(q))]
            if bad:
                print(
                    f"   FAIL {symbol}: trace differs from the character on {len(bad)} elements"
                )
                ok = False
            else:
                print(
                    f"   ok   {symbol}: {len(D)} matrices, traces are the character (dimension {ir.dimension})"
                )
    return ok


def main(argv):
    if len(argv) < 2 or argv[1] not in ("check", "generate"):
        print(__doc__)
        return 2
    if argv[1] == "generate":
        blocks = {group: build(group) for group in COVER}
        OUT.write_text(render(blocks))
        print(f"wrote {OUT}")
    return 0 if verify() else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
