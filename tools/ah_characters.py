"""Transcription of the Altmann & Herzig character tables and sample matrices
(wayfinder ticket #19).

Two sources feed one set of reference files:

* the **typeset** character tables of the master groups and of the groups that
  own a page of their own (41 groups; four of them, ``S18``, ``S20``, ``C9h``,
  ``C10h``, continue on a second page listed in ``CONTINUED``) are parsed
  from the OCR text of the PDF (pymupdf ``get_text`` of the table T n.4);
* the **scanned** group pages (32 groups) are read from the page images and
  typed into ``BOOK`` below, in the book's own notation transliterated to
  ASCII (``e`` = epsilon, ``d`` = delta, ``h`` = eta, ``t`` = theta, ``2c5^2`` =
  2 cos(2 pi / 5), ``sqrt2``, ``ie*``); the generator matrices of ten axial
  groups (T n.7; six for ticket #19, ``D5``, ``D8``, ``D9``, ``D10`` for
  ticket #22) into ``MATRICES`` the same way.

Ticket #19 transcribed the 39 groups of the crystallographic range plus
``C5v``, ``D5h``, ``D5d``, ``S8``, ``D4d``; ticket #36 added the remaining 34
axial groups with ``n <= 10`` so that the label convention's irrep symbols
are decided by the book for every group in its declared range.

    uv run python tools/ah_characters.py generate [OCR_DIR]   # write tests/reference/altmann-herzig/{characters,matrices}/
    uv run python tools/ah_characters.py check                # structural comparison with the computed tables (REF-1)

``OCR_DIR`` holds one ``pNNN.txt`` per PDF page; without it only the scanned
groups (typed in ``BOOK``) are generated. The check rebuilds every group with
the prototypes, matches the book's classes to ours by size, element order and
parity and the book's irreps by their characters (projective rows up to one
sign per class, which is the book's tilde assignment), and compares the
book's reality column ``tau`` with the Frobenius-Schur indicator. Nothing here
is production code.
"""

from __future__ import annotations

import math
import re
import sys
import tomllib
from fractions import Fraction as Fr
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
sys.path.insert(0, str(ROOT / "tools"))
import abstract_group_prototype as ag
import canonical_basis_prototype as cb
import koster_bethe  # noqa: F401  (extends the schema prototype to the cyclic and T families)
import output_schema_prototype as osp

SOURCE = "S. L. Altmann and P. Herzig, Point-Group Theory Tables, 2nd ed. (Wien: Clarendon, 2011)"
CHAR_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "characters"
MAT_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "matrices"

# group -> (book table number, pdf page of T n.4, typeset?)   (book page = pdf page - 12)
PAGES = {
    "C1": (1, 120, False), "C2": (2, 122, False), "C3": (3, 124, False), "C4": (4, 126, False),
    "C6": (6, 131, False), "Ci": (11, 150, False), "Cs": (12, 152, False), "S4": (13, 156, False),
    "S6": (14, 158, False), "S8": (15, 162, True), "D2": (22, 206, False), "D3": (23, 208, False),
    "D4": (24, 211, False), "D6": (26, 219, False), "D2h": (31, 259, True), "D3h": (32, 262, False),
    "D4h": (33, 270, True), "D5h": (34, 276, True), "D6h": (35, 289, True), "D2d": (41, 378, False),
    "D3d": (42, 383, True), "D4d": (43, 387, False), "D5d": (44, 395, True), "C2v": (50, 494, False),
    "C3v": (51, 497, True), "C4v": (52, 501, False), "C5v": (53, 504, False), "C6v": (54, 509, False),
    "C2h": (60, 544, False), "C3h": (61, 546, False), "C4h": (62, 550, True), "C6h": (64, 558, True),
    "O": (69, 592, False), "T": (70, 602, False), "Oh": (71, 618, True), "Th": (72, 645, True),
    "Td": (73, 649, False), "I": (74, 664, True), "Ih": (75, 677, True),
    # ticket #36: the remaining axial groups with n <= 10
    "C5": (5, 128, False), "C7": (7, 134, False), "C8": (8, 138, True), "C9": (9, 141, True),
    "C10": (10, 145, True), "S10": (16, 165, True), "S12": (17, 169, True), "S14": (18, 174, True),
    "S16": (19, 179, True), "S18": (20, 186, True), "S20": (21, 194, True), "D5": (25, 215, False),
    "D7": (27, 225, False), "D8": (28, 232, False), "D9": (29, 240, True), "D10": (30, 248, True),
    "D7h": (36, 300, True), "D8h": (37, 320, True), "D9h": (38, 332, True), "D10h": (39, 362, True),
    "D6d": (45, 404, True), "D7d": (46, 419, True), "D8d": (47, 430, True), "D9d": (48, 454, True),
    "D10d": (49, 467, True), "C7v": (55, 513, False), "C8v": (56, 519, False), "C9v": (57, 523, True),
    "C10v": (58, 532, True), "C5h": (63, 554, True), "C7h": (65, 563, True), "C8h": (66, 569, True),
    "C9h": (67, 575, True), "C10h": (68, 583, True),
}  # fmt: skip

# group -> pdf page of "T n.4 Character table (cont.)": the remaining class columns of a
# table too wide for one page; the rows repeat their labels.
CONTINUED = {"S18": 187, "S20": 195, "C9h": 576, "C10h": 584}

# --------------------------------------------------------------------------
# Hand transcriptions of the scanned pages (read from the page images)
# --------------------------------------------------------------------------
#
# classes: the book's column headers, ASCII (sigma_h, C3+, S4-, C2', 2S8^3 ...)
# rows: (book irrep label as printed, characters, tau)
# eps / delta: the book's definition of epsilon / delta as the argument of exp(2 pi i / n)


def tbl(**fields: Any) -> dict[str, Any]:
    return fields


BOOK: dict[str, dict[str, Any]] = {
    "C4v": tbl(
        classes="E 2C4 C2 2sigma_v 2sigma_d",
        rows=[
            ("A1", "1 1 1 1 1", "a"),
            ("A2", "1 1 1 -1 -1", "a"),
            ("B1", "1 -1 1 1 -1", "a"),
            ("B2", "1 -1 1 -1 1", "a"),
            ("E", "2 0 -2 0 0", "a"),
            ("E1/2", "2 sqrt2 0 0 0", "c"),
            ("E3/2", "2 -sqrt2 0 0 0", "c"),
        ],
    ),
    "S4": tbl(
        classes="E S4- C2 S4+",
        eps=8,
        rows=[
            ("A", "1 1 1 1", "a"),
            ("B", "1 -1 1 -1", "a"),
            ("1E", "1 -i -1 i", "b"),
            ("2E", "1 i -1 -i", "b"),
            ("1E1/2", "1 e* -i e", "b"),
            ("2E1/2", "1 e i e*", "b"),
            ("1E3/2", "1 -e* -i -e", "b"),
            ("2E3/2", "1 -e i -e*", "b"),
        ],
    ),
    "T": tbl(
        classes="E 3C2 4C3+ 4C3-",
        eps=3,
        rows=[
            ("A", "1 1 1 1", "a"),
            ("1E", "1 1 e e*", "b"),
            ("2E", "1 1 e* e", "b"),
            ("T", "3 -1 0 0", "a"),
            ("E1/2", "2 0 1 1", "c"),
            ("1F3/2", "2 0 e e*", "b"),
            ("2F3/2", "2 0 e* e", "b"),
        ],
    ),
    "C6": tbl(
        classes="E C6+ C3+ C2 C3- C6-",
        eps=3,
        rows=[
            ("A", "1 1 1 1 1 1", "a"),
            ("B", "1 -1 1 -1 1 -1", "a"),
            ("1E1", "1 -e e* -1 e -e*", "b"),
            ("2E1", "1 -e* e -1 e* -e", "b"),
            ("1E2", "1 e e* 1 e e*", "b"),
            ("2E2", "1 e* e 1 e* e", "b"),
            ("1E1/2", "1 -ie -e* i -e ie*", "b"),
            ("2E1/2", "1 ie* -e -i -e* -ie", "b"),
            ("1E3/2", "1 -i -1 i -1 i", "b"),
            ("2E3/2", "1 i -1 -i -1 -i", "b"),
            ("1E5/2", "1 -ie* -e i -e* ie", "b"),
            ("2E5/2", "1 ie -e* -i -e -ie*", "b"),
        ],
    ),
    "D4d": tbl(
        classes="E 2C4 C2 4C2' 2S8^3 2S8 4sigma_d",
        rows=[
            ("A1", "1 1 1 1 1 1 1", "a"),
            ("A2", "1 1 1 -1 1 1 -1", "a"),
            ("B1", "1 1 1 1 -1 -1 -1", "a"),
            ("B2", "1 1 1 -1 -1 -1 1", "a"),
            ("E1", "2 0 -2 0 -sqrt2 sqrt2 0", "a"),
            ("E2", "2 -2 2 0 0 0 0", "a"),
            ("E3", "2 0 -2 0 sqrt2 -sqrt2 0", "a"),
            ("E1/2", "2 sqrt2 0 0 2c8 2c8^3 0", "c"),
            ("E3/2", "2 -sqrt2 0 0 2c8^3 -2c8 0", "c"),
            ("E5/2", "2 -sqrt2 0 0 -2c8^3 2c8 0", "c"),
            ("E7/2", "2 sqrt2 0 0 -2c8 -2c8^3 0", "c"),
        ],
    ),
}

BOOK.update(
    {
        "C1": tbl(classes="E", rows=[("A", "1", "a"), ("A1/2", "1", "a")]),
        "C2": tbl(
            classes="E C2",
            rows=[
                ("A", "1 1", "a"),
                ("B", "1 -1", "a"),
                ("1E1/2", "1 i", "b"),
                ("2E1/2", "1 -i", "b"),
            ],
        ),
        "C3": tbl(
            classes="E C3+ C3-",
            eps=3,
            rows=[
                ("A", "1 1 1", "a"),
                ("1E", "1 e* e", "b"),
                ("2E", "1 e e*", "b"),
                ("1E1/2", "1 -e* -e", "b"),
                ("2E1/2", "1 -e -e*", "b"),
                ("A3/2", "1 -1 -1", "a"),
            ],
        ),
        "C4": tbl(
            classes="E C4+ C2 C4-",
            eps=8,
            rows=[
                ("A", "1 1 1 1", "a"),
                ("B", "1 -1 1 -1", "a"),
                ("1E", "1 -i -1 i", "b"),
                ("2E", "1 i -1 -i", "b"),
                ("1E1/2", "1 e* -i e", "b"),
                ("2E1/2", "1 e i e*", "b"),
                ("1E3/2", "1 -e* -i -e", "b"),
                ("2E3/2", "1 -e i -e*", "b"),
            ],
        ),
        "Ci": tbl(
            classes="E i",
            rows=[
                ("Ag", "1 1", "a"),
                ("Au", "1 -1", "a"),
                ("A1/2,g", "1 1", "a"),
                ("A1/2,u", "1 -1", "a"),
            ],
        ),
        "C2h": tbl(
            classes="E C2 i sigma_h",
            rows=[
                ("Ag", "1 1 1 1", "a"),
                ("Bg", "1 -1 1 -1", "a"),
                ("Au", "1 1 -1 -1", "a"),
                ("Bu", "1 -1 -1 1", "a"),
                ("1E1/2,g", "1 i 1 i", "b"),
                ("2E1/2,g", "1 -i 1 -i", "b"),
                ("1E1/2,u", "1 i -1 -i", "b"),
                ("2E1/2,u", "1 -i -1 i", "b"),
            ],
        ),
        "C2v": tbl(
            classes="E C2 sigmax sigmay",
            rows=[
                ("A1", "1 1 1 1", "a"),
                ("A2", "1 1 -1 -1", "a"),
                ("B1", "1 -1 -1 1", "a"),
                ("B2", "1 -1 1 -1", "a"),
                ("E1/2", "2 0 0 0", "c"),
            ],
        ),
        "C3h": tbl(
            classes="E C3+ C3- sigma_h S3+ S3-",
            eps=3,
            rows=[
                ("A'", "1 1 1 1 1 1", "a"),
                ("1E'", "1 e* e 1 e* e", "b"),
                ("2E'", "1 e e* 1 e e*", "b"),
                ("A''", "1 1 1 -1 -1 -1", "a"),
                ("1E''", "1 e* e -1 -e* -e", "b"),
                ("2E''", "1 e e* -1 -e -e*", "b"),
                ("1E1/2", "1 -e* -e i ie* -ie", "b"),
                ("2E1/2", "1 -e -e* -i -ie ie*", "b"),
                ("1E3/2", "1 -1 -1 i i -i", "b"),
                ("2E3/2", "1 -1 -1 -i -i i", "b"),
                ("1E5/2", "1 -e -e* i ie -ie*", "b"),
                ("2E5/2", "1 -e* -e -i -ie* ie", "b"),
            ],
        ),
    }
)

BOOK.update(
    {
        "Cs": tbl(
            classes="E sigma_h",
            rows=[
                ("A'", "1 1", "a"),
                ("A''", "1 -1", "a"),
                ("1E1/2", "1 i", "b"),
                ("2E1/2", "1 -i", "b"),
            ],
        ),
        "S6": tbl(
            classes="E C3+ C3- i S6- S6+",
            eps=3,
            rows=[
                ("Ag", "1 1 1 1 1 1", "a"),
                ("1Eg", "1 e* e 1 e* e", "b"),
                ("2Eg", "1 e e* 1 e e*", "b"),
                ("Au", "1 1 1 -1 -1 -1", "a"),
                ("1Eu", "1 e* e -1 -e* -e", "b"),
                ("2Eu", "1 e e* -1 -e -e*", "b"),
                ("1E1/2,g", "1 -e* -e 1 -e* -e", "b"),
                ("2E1/2,g", "1 -e -e* 1 -e -e*", "b"),
                ("A3/2,g", "1 -1 -1 1 -1 -1", "a"),
                ("1E1/2,u", "1 -e* -e -1 e* e", "b"),
                ("2E1/2,u", "1 -e -e* -1 e e*", "b"),
                ("A3/2,u", "1 -1 -1 -1 1 1", "a"),
            ],
        ),
        "D2": tbl(
            classes="E C2z C2x C2y",
            rows=[
                ("A", "1 1 1 1", "a"),
                ("B1", "1 1 -1 -1", "a"),
                ("B2", "1 -1 -1 1", "a"),
                ("B3", "1 -1 1 -1", "a"),
                ("E1/2", "2 0 0 0", "c"),
            ],
        ),
        "D3": tbl(
            classes="E 2C3 3C2'",
            rows=[
                ("A1", "1 1 1", "a"),
                ("A2", "1 1 -1", "a"),
                ("E", "2 -1 0", "a"),
                ("E1/2", "2 1 0", "c"),
                ("1E3/2", "1 -1 i", "b"),
                ("2E3/2", "1 -1 -i", "b"),
            ],
        ),
        "D4": tbl(
            classes="E 2C4 C2 2C2' 2C2''",
            rows=[
                ("A1", "1 1 1 1 1", "a"),
                ("A2", "1 1 1 -1 -1", "a"),
                ("B1", "1 -1 1 1 -1", "a"),
                ("B2", "1 -1 1 -1 1", "a"),
                ("E", "2 0 -2 0 0", "a"),
                ("E1/2", "2 sqrt2 0 0 0", "c"),
                ("E3/2", "2 -sqrt2 0 0 0", "c"),
            ],
        ),
        "D6": tbl(
            classes="E 2C6 2C3 C2 3C2' 3C2''",
            rows=[
                ("A1", "1 1 1 1 1 1", "a"),
                ("A2", "1 1 1 1 -1 -1", "a"),
                ("B1", "1 -1 1 -1 1 -1", "a"),
                ("B2", "1 -1 1 -1 -1 1", "a"),
                ("E1", "2 1 -1 -2 0 0", "a"),
                ("E2", "2 -1 -1 2 0 0", "a"),
                ("E1/2", "2 sqrt3 1 0 0 0", "c"),
                ("E3/2", "2 0 -2 0 0 0", "c"),
                ("E5/2", "2 -sqrt3 1 0 0 0", "c"),
            ],
        ),
    }
)

BOOK.update(
    {
        "D2d": tbl(
            classes="E C2 2C2' 2S4 2sigma_d",
            rows=[
                ("A1", "1 1 1 1 1", "a"),
                ("A2", "1 1 -1 1 -1", "a"),
                ("B1", "1 1 1 -1 -1", "a"),
                ("B2", "1 1 -1 -1 1", "a"),
                ("E", "2 -2 0 0 0", "a"),
                ("E1/2", "2 0 0 sqrt2 0", "c"),
                ("E3/2", "2 0 0 -sqrt2 0", "c"),
            ],
        ),
        "D3h": tbl(
            classes="E 2C3 3C2' sigma_h 2S3 3sigma_v",
            rows=[
                ("A1'", "1 1 1 1 1 1", "a"),
                ("A2'", "1 1 -1 1 1 -1", "a"),
                ("E'", "2 -1 0 2 -1 0", "a"),
                ("A1''", "1 1 1 -1 -1 -1", "a"),
                ("A2''", "1 1 -1 -1 -1 1", "a"),
                ("E''", "2 -1 0 -2 1 0", "a"),
                ("E1/2", "2 1 0 0 sqrt3 0", "c"),
                ("E3/2", "2 -2 0 0 0 0", "c"),
                ("E5/2", "2 1 0 0 -sqrt3 0", "c"),
            ],
        ),
        "C5v": tbl(
            classes="E 2C5 2C5^2 5sigma_v",
            rows=[
                ("A1", "1 1 1 1", "a"),
                ("A2", "1 1 1 -1", "a"),
                ("E1", "2 2c5^2 2c5^4 0", "a"),
                ("E2", "2 2c5^4 2c5^2 0", "a"),
                ("E1/2", "2 -2c5^4 2c5^2 0", "c"),
                ("E3/2", "2 -2c5^2 2c5^4 0", "c"),
                ("1E5/2", "1 -1 1 i", "b"),
                ("2E5/2", "1 -1 1 -i", "b"),
            ],
        ),
        "C6v": tbl(
            classes="E 2C6 2C3 C2 3sigma_d 3sigma_v",
            rows=[
                ("A1", "1 1 1 1 1 1", "a"),
                ("A2", "1 1 1 1 -1 -1", "a"),
                ("B1", "1 -1 1 -1 -1 1", "a"),
                ("B2", "1 -1 1 -1 1 -1", "a"),
                ("E1", "2 1 -1 -2 0 0", "a"),
                ("E2", "2 -1 -1 2 0 0", "a"),
                ("E1/2", "2 sqrt3 1 0 0 0", "c"),
                ("E3/2", "2 0 -2 0 0 0", "c"),
                ("E5/2", "2 -sqrt3 1 0 0 0", "c"),
            ],
        ),
        "O": tbl(
            classes="E 3C2 8C3 6C4 6C2'",
            rows=[
                ("A1", "1 1 1 1 1", "a"),
                ("A2", "1 1 1 -1 -1", "a"),
                ("E", "2 2 -1 0 0", "a"),
                ("T1", "3 -1 0 1 -1", "a"),
                ("T2", "3 -1 0 -1 1", "a"),
                ("E1/2", "2 0 1 sqrt2 0", "c"),
                ("E5/2", "2 0 1 -sqrt2 0", "c"),
                ("F3/2", "4 0 -1 0 0", "c"),
            ],
        ),
        "Td": tbl(
            classes="E 3C2 8C3 6S4 6sigma_d",
            rows=[
                ("A1", "1 1 1 1 1", "a"),
                ("A2", "1 1 1 -1 -1", "a"),
                ("E", "2 2 -1 0 0", "a"),
                ("T1", "3 -1 0 1 -1", "a"),
                ("T2", "3 -1 0 -1 1", "a"),
                ("E1/2", "2 0 1 sqrt2 0", "c"),
                ("E5/2", "2 0 1 -sqrt2 0", "c"),
                ("F3/2", "4 0 -1 0 0", "c"),
            ],
        ),
    }
)

# Documented disagreements between the book's table and its own parameter, basis and
# matrix tables, applied in order (tests/altmann_herzig.py applies the same entries).
# swap_projective: two class columns whose projective-row entries the book prints interchanged.
# swap_labels: two irrep symbols the book prints on each other's rows.
DISCREPANCIES: dict[str, list[dict[str, Any]]] = {
    "D5h": [
        tbl(
            swap_projective=("2S5", "2S5^2"),
            note=(
                "The projective rows (E1/2 ... E9/2) of T 34.4 carry the values of the 2S5 and 2S5^2 "
                "columns interchanged: by the book's own parameters (T 39.1: S5+ has phi = 3pi/5, S5^2+ "
                "has phi = pi/5) the spin-1/2 character of S5 is 2cos(3pi/10) = 2c10^3, and the master "
                "table T 39.4 (D10h) prints E1/2,g = 2c10^3 on 2S5 and 2c10 on 2S5^2, whereas T 34.4 prints "
                "2c10 on 2S5 and 2c10^3 on 2S5^2. The ordinary rows are consistent with the header. "
                "With the two columns interchanged in the projective rows the table matches the computed one."
            ),
        ),
        tbl(
            swap_labels=("E3/2", "E7/2"),
            note=(
                "The rows E3/2 and E7/2 of T 34.4 carry each other's symbol. With the 2S5 / 2S5^2 columns "
                "corrected, the printed E3/2 row has 2c10 on 2S5 and -2c10^3 on 2S5^2 and the printed E7/2 row "
                "-2c10 on 2S5 and 2c10^3 on 2S5^2, whereas the book's own matrices T 34.7 (p. 265) give the irrep "
                "it labels E3/2 the matrices S5+ = diag(i delta, -i delta*) and S5^2+ = diag(-i epsilon, i epsilon*) "
                "(delta = exp(2 pi i/5), epsilon = exp(4 pi i/5)), of traces -2c10 and 2c10^3, and the basis table "
                "T 34.6 spans that E3/2 by |3/2 3/2>, whose characters are the same; the two irreps differ only by "
                "the sign on the improper classes (E7/2 = E3/2 x A1''). T 34.6 and T 34.7 agree with the canonical "
                "symbols, so the symbols E3/2 and E7/2 of T 34.4 are interchanged (ticket #37)."
            ),
        ),
    ],
}

# Generator matrices of T n.7 for the axial samples, read from the page images.
# elements: book labels of the generators; irreps: book label -> one matrix per element,
# entries in the token grammar above (the book's overbar is a minus sign). The D_n tables
# (ticket #22) name their roots of unity in `symbols` (d = delta, e = epsilon, h = eta,
# t = theta of the table's footnote, as exp(2 pi i * num / den)).
MATRICES: dict[str, dict[str, Any]] = {
    "D5": tbl(
        table="T 25.7, p. 204",
        eps=5,
        symbols={"d": (1, 5), "e": (2, 5)},
        elements=["C5+", "C2'1"],
        irreps={
            "E1": [[["d*", "0"], ["0", "d"]], [["0", "-1"], ["-1", "0"]]],
            "E2": [[["e*", "0"], ["0", "e"]], [["0", "1"], ["1", "0"]]],
            "E1/2": [[["-e", "0"], ["0", "-e*"]], [["0", "-i"], ["-i", "0"]]],
            "E3/2": [[["-d", "0"], ["0", "-d*"]], [["0", "i"], ["i", "0"]]],
        },
    ),
    "D8": tbl(
        table="T 28.7, p. 222",
        eps=8,
        symbols={"d": (1, 16), "e": (1, 8)},
        elements=["C8+", "C2'1"],
        irreps={
            "E1": [[["e*", "0"], ["0", "e"]], [["0", "-1"], ["-1", "0"]]],
            "E2": [[["-i", "0"], ["0", "i"]], [["0", "1"], ["1", "0"]]],
            "E3": [[["-e*", "0"], ["0", "-e"]], [["0", "-1"], ["-1", "0"]]],
            "E1/2": [[["d*", "0"], ["0", "d"]], [["0", "-i"], ["-i", "0"]]],
            "E3/2": [[["id*", "0"], ["0", "-id"]], [["0", "-i"], ["-i", "0"]]],
            "E5/2": [[["-id*", "0"], ["0", "id"]], [["0", "-i"], ["-i", "0"]]],
            "E7/2": [[["-d*", "0"], ["0", "-d"]], [["0", "-i"], ["-i", "0"]]],
        },
    ),
    "D9": tbl(
        table="T 29.7, p. 229",
        eps=9,
        symbols={"d": (1, 9), "e": (2, 9), "h": (3, 9), "t": (4, 9)},
        elements=["C9+", "C2'1"],
        irreps={
            "E1": [[["d*", "0"], ["0", "d"]], [["0", "-1"], ["-1", "0"]]],
            "E2": [[["e", "0"], ["0", "e*"]], [["0", "-1"], ["-1", "0"]]],
            "E3": [[["h*", "0"], ["0", "h"]], [["0", "-1"], ["-1", "0"]]],
            "E4": [[["t*", "0"], ["0", "t"]], [["0", "-1"], ["-1", "0"]]],
            "E1/2": [[["-t", "0"], ["0", "-t*"]], [["0", "-i"], ["-i", "0"]]],
            "E3/2": [[["-h", "0"], ["0", "-h*"]], [["0", "i"], ["i", "0"]]],
            "E5/2": [[["-e*", "0"], ["0", "-e"]], [["0", "-i"], ["-i", "0"]]],
            "E7/2": [[["-d", "0"], ["0", "-d*"]], [["0", "-i"], ["-i", "0"]]],
        },
    ),
    "D10": tbl(
        table="T 30.7, pp. 237-238",
        eps=5,
        symbols={"d": (1, 5), "e": (2, 5)},
        elements=["C10+", "C2'1"],
        irreps={
            "E1": [[["-e", "0"], ["0", "-e*"]], [["0", "-1"], ["-1", "0"]]],
            "E2": [[["d*", "0"], ["0", "d"]], [["0", "1"], ["1", "0"]]],
            "E3": [[["-d*", "0"], ["0", "-d"]], [["0", "1"], ["1", "0"]]],
            "E4": [[["e", "0"], ["0", "e*"]], [["0", "-1"], ["-1", "0"]]],
            "E1/2": [[["-id", "0"], ["0", "id*"]], [["0", "-i"], ["-i", "0"]]],
            "E3/2": [[["ie*", "0"], ["0", "-ie"]], [["0", "i"], ["i", "0"]]],
            "E5/2": [[["-i", "0"], ["0", "i"]], [["0", "-i"], ["-i", "0"]]],
            "E7/2": [[["-ie*", "0"], ["0", "ie"]], [["0", "i"], ["i", "0"]]],
            "E9/2": [[["id", "0"], ["0", "-id*"]], [["0", "-i"], ["-i", "0"]]],
        },
    ),
    "D3": tbl(
        table="T 23.7, p. 197",
        eps=3,
        elements=["C3+", "C2'1"],
        irreps={
            "E": [[["e*", "0"], ["0", "e"]], [["0", "-1"], ["-1", "0"]]],
            "E1/2": [[["-e", "0"], ["0", "-e*"]], [["0", "-i"], ["-i", "0"]]],
        },
    ),
    "C4v": tbl(
        table="T 52.7, p. 490",
        eps=8,
        elements=["C4+", "sigma_v1"],
        irreps={
            "E": [[["-i", "0"], ["0", "i"]], [["0", "-1"], ["-1", "0"]]],
            "E1/2": [[["e*", "0"], ["0", "e"]], [["0", "-i"], ["-i", "0"]]],
            "E3/2": [[["-e*", "0"], ["0", "-e"]], [["0", "-i"], ["-i", "0"]]],
        },
    ),
    "D3h": tbl(
        table="T 32.7, p. 252",
        eps=3,
        elements=["C3+", "C2'1", "sigma_h"],
        irreps={
            "E'": [
                [["e*", "0"], ["0", "e"]],
                [["0", "-1"], ["-1", "0"]],
                [["1", "0"], ["0", "1"]],
            ],
            "E''": [
                [["e*", "0"], ["0", "e"]],
                [["0", "-1"], ["-1", "0"]],
                [["-1", "0"], ["0", "-1"]],
            ],
            "E1/2": [
                [["-e", "0"], ["0", "-e*"]],
                [["0", "-i"], ["-i", "0"]],
                [["-i", "0"], ["0", "i"]],
            ],
            "E3/2": [
                [["-1", "0"], ["0", "-1"]],
                [["0", "i"], ["i", "0"]],
                [["i", "0"], ["0", "-i"]],
            ],
            "E5/2": [
                [["-e", "0"], ["0", "-e*"]],
                [["0", "-i"], ["-i", "0"]],
                [["i", "0"], ["0", "-i"]],
            ],
        },
    ),
    "D2d": tbl(
        table="T 41.7, p. 367",
        eps=8,
        elements=["S4-", "C2'1"],
        irreps={
            "E": [[["-i", "0"], ["0", "i"]], [["0", "-1"], ["-1", "0"]]],
            "E1/2": [[["e*", "0"], ["0", "e"]], [["0", "-i"], ["-i", "0"]]],
            "E3/2": [[["-e*", "0"], ["0", "-e"]], [["0", "-i"], ["-i", "0"]]],
        },
    ),
    "D6h": tbl(
        table="T 35.7, pp. 279-280",
        eps=3,
        elements=["C6+", "C2'1", "i"],
        irreps={
            "E1g": [
                [["-e", "0"], ["0", "-e*"]],
                [["0", "-1"], ["-1", "0"]],
                [["1", "0"], ["0", "1"]],
            ],
            "E2g": [
                [["e", "0"], ["0", "e*"]],
                [["0", "-1"], ["-1", "0"]],
                [["1", "0"], ["0", "1"]],
            ],
            "E1u": [
                [["-e", "0"], ["0", "-e*"]],
                [["0", "-1"], ["-1", "0"]],
                [["-1", "0"], ["0", "-1"]],
            ],
            "E2u": [
                [["e", "0"], ["0", "e*"]],
                [["0", "-1"], ["-1", "0"]],
                [["-1", "0"], ["0", "-1"]],
            ],
            "E1/2,g": [
                [["ie*", "0"], ["0", "-ie"]],
                [["0", "-i"], ["-i", "0"]],
                [["1", "0"], ["0", "1"]],
            ],
            "E3/2,g": [
                [["-i", "0"], ["0", "i"]],
                [["0", "i"], ["i", "0"]],
                [["1", "0"], ["0", "1"]],
            ],
            "E5/2,g": [
                [["-ie*", "0"], ["0", "ie"]],
                [["0", "-i"], ["-i", "0"]],
                [["1", "0"], ["0", "1"]],
            ],
            "E1/2,u": [
                [["ie*", "0"], ["0", "-ie"]],
                [["0", "-i"], ["-i", "0"]],
                [["-1", "0"], ["0", "-1"]],
            ],
            "E3/2,u": [
                [["-i", "0"], ["0", "i"]],
                [["0", "i"], ["i", "0"]],
                [["-1", "0"], ["0", "-1"]],
            ],
            "E5/2,u": [
                [["-ie*", "0"], ["0", "ie"]],
                [["0", "-i"], ["-i", "0"]],
                [["-1", "0"], ["0", "-1"]],
            ],
        },
    ),
    "C6v": tbl(
        table="T 54.7, p. 499",
        eps=3,
        elements=["C6+", "sigma_v1"],
        irreps={
            "E1": [[["-e", "0"], ["0", "-e*"]], [["0", "1"], ["1", "0"]]],
            "E2": [[["e", "0"], ["0", "e*"]], [["0", "-1"], ["-1", "0"]]],
            "E1/2": [[["ie*", "0"], ["0", "-ie"]], [["0", "-1"], ["1", "0"]]],
            "E3/2": [[["-i", "0"], ["0", "i"]], [["0", "-1"], ["1", "0"]]],
            "E5/2": [[["-ie*", "0"], ["0", "ie"]], [["0", "1"], ["-1", "0"]]],
        },
    ),
}

# Ticket #36: the scanned pages of the remaining axial groups with n <= 10. `symbols` is
# the footnote of the table, letter -> (numerator, denominator) of exp(2 pi i * num / den).
BOOK.update(
    {
        "C5": tbl(
            classes="E C5+ C5^2+ C5^2- C5-",
            symbols={"d": (1, 5), "e": (2, 5)},
            rows=[
                ("A", "1 1 1 1 1", "a"),
                ("1E1", "1 d* e* e d", "b"),
                ("2E1", "1 d e e* d*", "b"),
                ("1E2", "1 e* d d* e", "b"),
                ("2E2", "1 e d* d e*", "b"),
                ("1E1/2", "1 -e* d d* -e", "b"),
                ("2E1/2", "1 -e d* d -e*", "b"),
                ("1E3/2", "1 -d e e* -d*", "b"),
                ("2E3/2", "1 -d* e* e -d", "b"),
                ("A5/2", "1 -1 1 1 -1", "a"),
            ],
        ),
        "C7": tbl(
            classes="E C7+ C7^2+ C7^3+ C7^3- C7^2- C7-",
            symbols={"d": (1, 7), "e": (2, 7), "h": (3, 7)},
            rows=[
                ("A", "1 1 1 1 1 1 1", "a"),
                ("1E1", "1 d* e* h* h e d", "b"),
                ("2E1", "1 d e h h* e* d*", "b"),
                ("1E2", "1 e* h d d* h* e", "b"),
                ("2E2", "1 e h* d* d h e*", "b"),
                ("1E3", "1 h* d e* e d* h", "b"),
                ("2E3", "1 h d* e e* d h*", "b"),
                ("1E1/2", "1 -h* d -e* -e d* -h", "b"),
                ("2E1/2", "1 -h d* -e -e* d -h*", "b"),
                ("1E3/2", "1 -e h* -d* -d h -e*", "b"),
                ("2E3/2", "1 -e* h -d -d* h* -e", "b"),
                ("1E5/2", "1 -d* e* -h* -h e -d", "b"),
                ("2E5/2", "1 -d e -h -h* e* -d*", "b"),
                ("A7/2", "1 -1 1 -1 -1 1 -1", "a"),
            ],
        ),
        "D5": tbl(
            classes="E 2C5 2C5^2 5C2'",
            rows=[
                ("A1", "1 1 1 1", "a"),
                ("A2", "1 1 1 -1", "a"),
                ("E1", "2 2c5^2 2c5^4 0", "a"),
                ("E2", "2 2c5^4 2c5^2 0", "a"),
                ("E1/2", "2 -2c5^4 2c5^2 0", "c"),
                ("E3/2", "2 -2c5^2 2c5^4 0", "c"),
                ("1E5/2", "1 -1 1 i", "b"),
                ("2E5/2", "1 -1 1 -i", "b"),
            ],
        ),
        "D7": tbl(
            classes="E 2C7 2C7^2 2C7^3 7C2'",
            rows=[
                ("A1", "1 1 1 1 1", "a"),
                ("A2", "1 1 1 1 -1", "a"),
                ("E1", "2 2c7^2 2c7^4 2c7^6 0", "a"),
                ("E2", "2 2c7^4 2c7^6 2c7^2 0", "a"),
                ("E3", "2 2c7^6 2c7^2 2c7^4 0", "a"),
                ("E1/2", "2 -2c7^6 2c7^2 -2c7^4 0", "c"),
                ("E3/2", "2 -2c7^4 2c7^6 -2c7^2 0", "c"),
                ("E5/2", "2 -2c7^2 2c7^4 -2c7^6 0", "c"),
                ("1E7/2", "1 -1 1 -1 i", "b"),
                ("2E7/2", "1 -1 1 -1 -i", "b"),
            ],
        ),
        "D8": tbl(
            classes="E 2C8 2C4 2C8^3 C2 4C2' 4C2''",
            rows=[
                ("A1", "1 1 1 1 1 1 1", "a"),
                ("A2", "1 1 1 1 1 -1 -1", "a"),
                ("B1", "1 -1 1 -1 1 1 -1", "a"),
                ("B2", "1 -1 1 -1 1 -1 1", "a"),
                ("E1", "2 sqrt2 0 -sqrt2 -2 0 0", "a"),
                ("E2", "2 0 -2 0 2 0 0", "a"),
                ("E3", "2 -sqrt2 0 sqrt2 -2 0 0", "a"),
                ("E1/2", "2 2c8 sqrt2 2c8^3 0 0 0", "c"),
                ("E3/2", "2 2c8^3 -sqrt2 -2c8 0 0 0", "c"),
                ("E5/2", "2 -2c8^3 -sqrt2 2c8 0 0 0", "c"),
                ("E7/2", "2 -2c8 sqrt2 -2c8^3 0 0 0", "c"),
            ],
        ),
        "C7v": tbl(
            classes="E 2C7 2C7^2 2C7^3 7sigma_v",
            rows=[
                ("A1", "1 1 1 1 1", "a"),
                ("A2", "1 1 1 1 -1", "a"),
                ("E1", "2 2c7^2 2c7^4 2c7^6 0", "a"),
                ("E2", "2 2c7^4 2c7^6 2c7^2 0", "a"),
                ("E3", "2 2c7^6 2c7^2 2c7^4 0", "a"),
                ("E1/2", "2 -2c7^6 2c7^2 -2c7^4 0", "c"),
                ("E3/2", "2 -2c7^4 2c7^6 -2c7^2 0", "c"),
                ("E5/2", "2 -2c7^2 2c7^4 -2c7^6 0", "c"),
                ("1E7/2", "1 -1 1 -1 i", "b"),
                ("2E7/2", "1 -1 1 -1 -i", "b"),
            ],
        ),
        "C8v": tbl(
            classes="E 2C8 2C4 2C8^3 C2 4sigma_v 4sigma_d",
            rows=[
                ("A1", "1 1 1 1 1 1 1", "a"),
                ("A2", "1 1 1 1 1 -1 -1", "a"),
                ("B1", "1 -1 1 -1 1 1 -1", "a"),
                ("B2", "1 -1 1 -1 1 -1 1", "a"),
                ("E1", "2 sqrt2 0 -sqrt2 -2 0 0", "a"),
                ("E2", "2 0 -2 0 2 0 0", "a"),
                ("E3", "2 -sqrt2 0 sqrt2 -2 0 0", "a"),
                ("E1/2", "2 2c8 sqrt2 2c8^3 0 0 0", "c"),
                ("E3/2", "2 2c8^3 -sqrt2 -2c8 0 0 0", "c"),
                ("E5/2", "2 -2c8^3 -sqrt2 2c8 0 0 0", "c"),
                ("E7/2", "2 -2c8 sqrt2 -2c8^3 0 0 0", "c"),
            ],
        ),
    }
)

# --------------------------------------------------------------------------
# Values: book tokens -> exact field elements
# --------------------------------------------------------------------------

TOKEN_RE = re.compile(
    r"^(-)?(\d+)?(i)?(e\*|e|d\*|d|h\*|h|t\*|t|z\*|z|sqrt2|sqrt3|c(\d+)(?:\^(\d+))?)?$"
)


def _symbols(eps, delta, symbols):
    """letter -> (numerator, denominator): the book's root of unity ``exp(2 pi i * num / den)``.

    ``e`` and ``d`` default to ``exp(2 pi i / eps)`` and ``exp(2 pi i / delta)``;
    a table with more symbols (``h`` = eta, ``t`` = theta) or with numerators
    (``e = exp(4 pi i / 5)``) passes its own map.
    """
    out = {"e": (1, eps), "d": (1, delta)}
    if symbols:
        out.update(symbols)
    return out


SYMBOL_NAMES = {"e": "epsilon", "d": "delta", "h": "eta", "t": "theta", "z": "zeta"}


def book_symbols_of(book):
    """The letters of a BOOK table as ``symbols`` (``eps`` / ``delta`` ints folded in)."""
    out = {}
    if book.get("eps"):
        out["e"] = (1, book["eps"])
    if book.get("delta"):
        out["d"] = (1, book["delta"])
    out.update(book.get("symbols") or {})
    return out


def symbol_text(num, den):
    return f"exp(2 pi i / {den})" if num == 1 else f"exp(2 pi i * {num} / {den})"


def parse_symbol_text(text):
    m = re.fullmatch(r"exp\(2 pi i(?: \* (\d+))? / (\d+)\)", text)
    assert m, text
    return (int(m.group(1) or 1), int(m.group(2)))


def _int_after(pattern, text):
    m = re.search(pattern, text)
    assert m, (pattern, text)
    return int(m.group(1))


def token_conductor(tok, eps, delta, symbols=None):
    m = TOKEN_RE.match(tok)
    assert m, tok
    need = 1
    if m.group(3):
        need = math.lcm(need, 4)
    base = m.group(4)
    syms = _symbols(eps, delta, symbols)
    if base and base.rstrip("*") in syms:
        need = math.lcm(need, syms[base.rstrip("*")][1])
    elif base == "sqrt2":
        need = math.lcm(need, 8)
    elif base == "sqrt3":
        need = math.lcm(need, 12)
    elif base and base.startswith("c"):
        need = math.lcm(need, 2 * int(m.group(5)))
    return need


def token_value(F, tok, eps, delta, symbols=None):
    m = TOKEN_RE.match(tok)
    assert m, tok
    sign, coeff, has_i, base = m.group(1), m.group(2), m.group(3), m.group(4)
    x = F.rational(int(coeff) if coeff else 1)
    syms = _symbols(eps, delta, symbols)
    if base and base.rstrip("*") in syms:
        num, n = syms[base.rstrip("*")]
        assert F.conductor % n == 0, (tok, F.conductor)
        z = F.zeta(F.conductor // n * num)
        x = F.mul(x, F.conj(z) if base.endswith("*") else z)
    elif base == "sqrt2":
        x = F.mul(x, F.sqrt2())
    elif base == "sqrt3":
        x = F.mul(x, cb.const(F, "sqrt3"))
    elif base and base.startswith("c"):
        n, mm = int(m.group(5)), int(m.group(6) or 1)
        z = F.zeta(
            F.conductor // (2 * n) * mm
        )  # cos(m pi / n) = (zeta_2n^m + zeta_2n^-m) / 2
        x = F.mul(x, F.scale(F.add(z, F.conj(z)), Fr(1, 2)))
    if has_i:
        x = F.mul(x, F.i())
    return F.neg(x) if sign else x


# --------------------------------------------------------------------------
# Class labels -> invariants
# --------------------------------------------------------------------------

CLASS_RE = re.compile(
    r"^(\d+)?(E|i|sigma_[hvd]|sigma[xyz]?|C(\d+)(?:\^(\d+))?|S(\d+)(?:\^(\d+))?)([+-]|'|''|[xyz])?(\d)?$"
)


def class_invariants(label):
    """(size, element order, parity) of a book class column."""
    m = CLASS_RE.match(label)
    assert m, label
    size = int(m.group(1) or 1)
    body = m.group(2)
    if body == "E":
        return size, 1, 1
    if body == "i":
        return size, 2, -1
    if body.startswith("sigma"):
        return size, 2, -1
    if body.startswith("C"):
        n, k = int(m.group(3)), int(m.group(4) or 1)
        return size, n // math.gcd(n, k), 1
    n, k = int(m.group(5)), int(m.group(6) or 1)
    return size, math.lcm(n // math.gcd(n, k), 2), -1


# --------------------------------------------------------------------------
# Parsing the typeset OCR text
# --------------------------------------------------------------------------

GREEK = {
    "ϵ": "e",
    "δ": "d",
    "η": "h",
    "θ": "t",
    "ζ": "z",
    "σ": "sigma",
    "′": "'",
    "″": "''",
    "−": "-",
    "∗": "*",
}


def normalise(text):
    for a, b in GREEK.items():
        text = text.replace(a, b)
    return text


def parse_typeset(group, text):
    """Parse the table T n.4 of a typeset page into the BOOK format."""
    text = normalise(text)
    # the OCR glues a starred value to a following negative one ("θ∗−iϵ" on the S16 and
    # C8h pages); a star never precedes a minus inside one value
    text = text.replace("*-", "* -")
    seg = text[text.index("Character table") :]
    seg = seg.split("Cartesian tensors")[0]
    toks = seg.replace("\n", " ").split()
    # constants: the footnote "d = exp(2πi/9), e = exp(4πi/9), ..." -> letter -> (num, den)
    symbols = {}
    for k, t in enumerate(toks):
        if t in SYMBOL_NAMES and toks[k + 1] == "=" and toks[k + 2].startswith("exp("):
            m = re.match(r"exp\((\d+)πi/(\d+)\)", toks[k + 2])
            assert m, toks[k + 2]
            num, den = int(m.group(1)), int(m.group(2))
            assert num % 2 == 0, toks[k + 2]
            symbols[t] = (num // 2, den)
    # header: from the group symbol to "tau"
    start = toks.index(group) + 1 if group in toks else toks.index("E") - 0
    tau = toks.index("τ")
    header = toks[start:tau]
    classes = merge_class_tokens(header)
    body = toks[tau + 1 :]
    rows = []
    k = 0
    n = len(classes)
    while k + n + 1 < len(body):
        t = body[k]
        if t in SYMBOL_NAMES and k + 1 < len(body) and body[k + 1] == "=":
            break
        if t.startswith("cm") or t.startswith("c") and t.endswith("m"):
            break
        if t == "T" and re.match(r"^\d+\.\d$", body[k + 1]):
            break
        # "A' 1 ..." is the label A1' followed by its values in D5h but the label A'
        # followed by the value 1 in C5h: try the merged label first, then the bare
        # one, and keep the reading whose row ends on a reality mark
        merged = merge_label_tokens(body, k)
        candidates = [merged] + ([(body[k], k + 1)] if merged[1] == k + 2 else [])
        for label, kk in candidates:
            vals = []
            try:
                while len(vals) < n:
                    v, kk = merge_value_tokens(body, kk)
                    vals.append(v)
                tau_v = body[kk]
            except IndexError:
                continue
            if tau_v in ("a", "b", "c"):
                break
        else:
            raise AssertionError((group, body[k : k + n + 3]))
        k = kk + 1
        rows.append((label, " ".join(vals), tau_v))
    return tbl(classes=" ".join(classes), symbols=symbols or None, rows=rows)


def merge_continued(first, second):
    """One table from the two pages of a wide table: the same rows, more columns."""
    assert [r[0] for r in first["rows"]] == [r[0] for r in second["rows"]], (
        [r[0] for r in first["rows"]],
        [r[0] for r in second["rows"]],
    )
    assert first.get("symbols") == second.get("symbols")
    rows = [
        (lab, f"{vals} {more}", tau)
        for (lab, vals, tau), (_, more, _) in zip(
            first["rows"], second["rows"], strict=True
        )
    ]
    return tbl(
        classes=f"{first['classes']} {second['classes']}",
        symbols=first.get("symbols"),
        rows=rows,
    )


def merge_class_tokens(header):
    """ "2C'" "2" -> "2C2'";  "C+" "6" -> "C6+";  "2C2" "5" -> "2C5^2";  "S3-" "8" -> "S8^3-";  "3sigma" -> "3sigma_h"."""
    out = []
    k = 0
    while k < len(header):
        t = header[k]
        nxt = header[k + 1] if k + 1 < len(header) else ""
        m = re.match(r"^(\d*)([CS])(\d*)([+-]?)('*)$", t)
        if m and nxt.isdigit():
            size, letter, sup, sign, primes = m.groups()
            n = nxt
            label = f"{size}{letter}{n}" + (f"^{sup}" if sup else "") + sign + primes
            out.append(label)
            k += 2
            continue
        t = re.sub(
            r"^(\d*)sigma([hvd]?)$",
            lambda mm: f"{mm.group(1)}sigma_{mm.group(2) or 'h'}",
            t,
        )
        out.append(t)
        k += 1
    return out


def merge_label_tokens(body, k):
    """ "A'" "1" -> "A1'";  "E'" "2" -> "E2'";  "1E1/2,g" stays."""
    t = body[k]
    nxt = body[k + 1] if k + 1 < len(body) else ""
    m = re.match(r"^([12]?[A-Z])('+|'')$", t)
    if m and nxt.isdigit():
        return f"{m.group(1)}{nxt}{m.group(2)}", k + 2
    return t, k + 1


def merge_value_tokens(body, k):
    """ "-" "sqrt" "2" / "√" "2" -> "sqrt2";  "2c2" "5" -> "2c5^2";  "2c10" -> "2c10"."""
    t = body[k]
    if t in ("-", "√", "-√", "sqrt"):
        parts = [t]
        while parts[-1] in ("-", "√", "-√", "sqrt"):
            k += 1
            parts.append(body[k])
        joined = "".join(parts).replace("√", "sqrt")
        return joined, k + 1
    m = re.match(r"^(-?)(\d*)c(\d+)$", t)
    if m:
        # "2c3" "10" is 2c10^3 (the superscript comes first in the OCR text); a bare
        # character never exceeds 4 here, so a following digit >= 5 is a denominator
        nxt = body[k + 1] if k + 1 < len(body) else ""
        if nxt.isdigit() and int(nxt) >= 5 and int(m.group(3)) < int(nxt):
            return f"{m.group(1)}{m.group(2)}c{nxt}^{m.group(3)}", k + 2
        return t, k + 1
    return t.replace("√", "sqrt"), k + 1


# --------------------------------------------------------------------------
# Computed tables
# --------------------------------------------------------------------------


class Computed:
    def __init__(self, symbol):
        self.ctx = osp.Context(symbol)
        ctx = self.ctx
        self.F = ctx.field
        self.classes = []  # (index, size, order, parity, representative element index)
        for cl in ctx.classes:
            rep = cl["representative"]
            g = ctx.G.elements[rep]
            self.classes.append(
                (
                    cl["index"],
                    cl["size"],
                    ag.element_order(ctx.G, g),
                    g[1],
                    rep,
                    cl["regular"],
                )
            )
        self.irreps = ctx.irreps


def lift_to(F, G, x):
    """Move x from field F to field G (conductor of G a multiple, or the fields equal)."""
    if F.conductor == G.conductor:
        return x
    return cb.lift(F, G, x)


def structural_match(group, book, comp):
    """Find class and irrep bijections; return (report lines, ok)."""
    lines = []
    symbols = book_symbols_of(book)
    Fb_cond = 1
    classes = book["classes"].split()
    for _, vals, _ in book["rows"]:
        for tok in vals.split():
            Fb_cond = math.lcm(Fb_cond, token_conductor(tok, None, None, symbols))
    M = math.lcm(Fb_cond, comp.F.conductor)
    F = cb.field(M)
    zero = F.rational(0)
    book_vals = [
        [token_value(F, tok, None, None, symbols) for tok in vals.split()]
        for _, vals, _ in book["rows"]
    ]
    inv_book = [class_invariants(c) for c in classes]
    inv_comp = [(size, order, parity) for _, size, order, parity, _, _ in comp.classes]
    if sorted(inv_book) != sorted(inv_comp):
        lines.append(
            f"  class invariants differ: book {sorted(inv_book)} vs computed {sorted(inv_comp)}"
        )
        return lines, False
    comp_rows = []
    for ir in comp.irreps:
        comp_rows.append(
            [lift_to(comp.F, F, ir.chars[rep]) for _, _, _, _, rep, _ in comp.classes]
        )
    n = len(classes)
    # candidate computed classes per book column
    cands = [[j for j in range(n) if inv_comp[j] == inv_book[i]] for i in range(n)]
    order_ = sorted(range(n), key=lambda i: len(cands[i]))
    proj_book = [("/" in lab) for lab, _, _ in book["rows"]]
    proj_comp = [ir.projective for ir in comp.irreps]
    # the book's rows in the order ordinary first: they fix the class bijection
    # without any sign, and the projective rows then fix one sign per class
    row_order = sorted(range(len(book_vals)), key=lambda r: proj_book[r])

    def match_rows(perm, pos, used, signs):
        """Pair the book's rows (from row_order[pos]) with computed irreps under the
        class bijection perm; a projective row may carry the opposite sign on a
        class, one sign per class (the book's tilde assignment), fixed by the first
        nonzero entry met and honoured afterwards. Returns {book row: irrep}."""
        if pos == len(row_order):
            return {}, signs
        r = row_order[pos]
        vals = book_vals[r]
        for s_, row in enumerate(comp_rows):
            if s_ in used or proj_comp[s_] != proj_book[r]:
                continue
            new_signs = signs
            ok = True
            for i in range(n):
                j = perm[i]
                v, b = row[j], vals[i]
                if not proj_book[r] or b == zero:
                    if v != b:
                        ok = False
                        break
                    continue
                if v == b:
                    want = 1
                elif v == F.neg(b):
                    want = -1
                else:
                    ok = False
                    break
                have = new_signs.get(j)
                if have is None:
                    if new_signs is signs:
                        new_signs = dict(signs)
                    new_signs[j] = want
                elif have != want:
                    ok = False
                    break
            if not ok:
                continue
            rest = match_rows(perm, pos + 1, used | {s_}, new_signs)
            if rest is not None:
                pairing, final_signs = rest
                pairing[r] = s_
                return pairing, final_signs
        return None

    def search(pos, perm, used):
        if pos == n:
            found = match_rows(perm, 0, frozenset(), {})
            if found is None:
                return None
            pairing, signs = found
            return (
                perm[:],
                [signs.get(j, 1) for j in range(n)],
                [pairing[r] for r in range(len(book_vals))],
            )
        i = order_[pos]
        for j in cands[i]:
            if j in used:
                continue
            perm[i] = j
            used.add(j)
            found = search(pos + 1, perm, used)
            if found:
                return found
            used.discard(j)
        return None

    found = search(0, [None] * n, set())
    if not found:
        lines.append(
            "  NO structural match: no class/irrep bijection reproduces the book's table"
        )
        return lines, False
    perm, signs, pairing = found
    lines.append(
        "  classes: "
        + ", ".join(
            f"{classes[i]} = {comp.ctx.classes[perm[i]]['label']}"
            + (" (sign -1)" if signs[perm[i]] < 0 else "")
            for i in range(n)
        )
    )
    lines.append(
        "  irreps:  "
        + ", ".join(
            f"{book['rows'][r][0]} = {comp.irreps[pairing[r]].symbol}"
            for r in range(len(pairing))
        )
    )
    ok = True
    fs_map = {"a": 1, "b": 0, "c": -1}
    for r, (lab, _, tau) in enumerate(book["rows"]):
        fs = comp.irreps[pairing[r]].frobenius_schur
        if fs_map[tau] != fs:
            lines.append(
                f"  tau of {lab} is {tau} but the Frobenius-Schur indicator of {comp.irreps[pairing[r]].symbol} is {fs}"
            )
            ok = False
    return lines, ok


# --------------------------------------------------------------------------
# Files
# --------------------------------------------------------------------------


def toml_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_characters(group, book, typeset):
    table, page, _ = PAGES[group]
    symbol = ag.ALIASES.get(group, group)
    symbols = book_symbols_of(book)
    M = 1
    for _, vals, _ in book["rows"]:
        for tok in vals.split():
            M = math.lcm(M, token_conductor(tok, None, None, symbols))
    F = cb.field(M)
    classes = book["classes"].split()
    lines = [
        f"# Altmann & Herzig, Point-Group Theory Tables (2nd ed.): character table of {group}",
        "# Format: tests/reference/altmann-herzig/characters/README.md",
        f"# {'Parsed from the OCR text of the typeset table' if typeset else 'Read from the page image'} by tools/ah_characters.py; verify against the print.",
        "",
        f"group = {toml_str(symbol)}",
        f'table = "T {table}.4, p. {page - 12}"',
        'status = "done"',
        f"source = {toml_str(SOURCE)}",
        f'transcribed_by = "agent (Claude) from the {"OCR text" if typeset else "page image"} of the PDF; verify against the print"',
        'notes = ""',
        f"conductor = {M}          # field of the `characters` values; zeta = exp(2 pi i / conductor)",
    ]
    for letter, name in SYMBOL_NAMES.items():
        if letter in symbols:
            lines.append(f'book_{name} = "{symbol_text(*symbols[letter])}"')
    for d in DISCREPANCIES.get(group, []):
        lines.append("")
        lines.append("[[known_discrepancy]]")
        for key in ("swap_projective", "swap_labels"):
            if key in d:
                lines.append(f"{key} = [{toml_str(d[key][0])}, {toml_str(d[key][1])}]")
        lines.append(f"note = {toml_str(d['note'])}")
    lines.append("")
    for c in classes:
        size, order, parity = class_invariants(c)
        lines += [
            "[[classes]]",
            f"label = {toml_str(c)}",
            f"size = {size}",
            f"order = {order}",
            f"parity = {parity}",
            "",
        ]
    for lab, vals, tau in book["rows"]:
        toks = vals.split()
        exact = [cb.to_grammar(F, token_value(F, t, None, None, symbols)) for t in toks]
        lines += [
            "[[irreps]]",
            f"label = {toml_str(lab)}",
            f"projective = {'true' if '/' in lab else 'false'}",
            f"book_values = [{', '.join(toml_str(t) for t in toks)}]",
            f"characters = [{', '.join(toml_str(x) for x in exact)}]",
            f"tau = {toml_str(tau)}",
            "",
        ]
    return "\n".join(lines) + "\n"


def element_class(label, classes):
    """The book's class column of a generator: strip the axis index and the sense."""
    root = re.sub(r"[+-]$", "", label)
    # the whole root first ("C10" of C10+), then without an axis index ("C2'" of C2'1)
    for cand in (root, re.sub(r"(sigma_[vd]|C\d+'*|S\d+(?:\^\d+)?)\d$", r"\1", root)):
        for c in classes:
            if re.sub(r"^\d+", "", c) == cand:
                return c
    raise KeyError((label, classes))


def write_matrices(group, mats):
    symbol = ag.ALIASES.get(group, group)
    eps = mats["eps"]
    symbols = mats.get("symbols")
    M = 1
    for irrep in mats["irreps"].values():
        for A in irrep:
            for row in A:
                for tok in row:
                    M = math.lcm(M, token_conductor(tok, eps, 1, symbols))
    F = cb.field(M)
    lines = [
        f"# Altmann & Herzig, Point-Group Theory Tables (2nd ed.): matrix representations of {group} (generators)",
        "# Format: tests/reference/altmann-herzig/matrices/README.md",
        "# Read from the page image by tools/ah_characters.py; verify against the print.",
        "",
        f"group = {toml_str(symbol)}",
        f"table = {toml_str(mats['table'])}",
        'status = "done"',
        f"source = {toml_str(SOURCE)}",
        'transcribed_by = "agent (Claude) from the page image of the PDF; verify against the print"',
        'notes = ""',
        f"conductor = {M}",
    ]
    if symbols:
        lines.append(
            "book_symbols = { "
            + ", ".join(
                f'{k} = "exp({2 * num} pi i / {den})"'
                for k, (num, den) in symbols.items()
            )
            + " }   # the table's footnote: d = delta, e = epsilon, h = eta, t = theta"
        )
    else:
        lines.append(f'book_epsilon = "exp(2 pi i / {eps})"')
    lines += [
        f"elements = [{', '.join(toml_str(e) for e in mats['elements'])}]   # book labels of the generators",
        "",
    ]
    for label, A_list in mats["irreps"].items():
        lines.append(f"[irreps.{toml_str(label)}]")
        lines.append("book_values = [")
        for A in A_list:
            lines.append(
                "  ["
                + ", ".join(
                    "[" + ", ".join(toml_str(t) for t in row) + "]" for row in A
                )
                + "],"
            )
        lines.append("]")
        lines.append("matrices = [")
        for A in A_list:
            rows = [
                "["
                + ", ".join(
                    toml_str(cb.to_grammar(F, token_value(F, t, eps, 1, symbols)))
                    for t in row
                )
                + "]"
                for row in A
            ]
            lines.append("  [" + ", ".join(rows) + "],")
        lines.append("]")
        lines.append("")
    return "\n".join(lines) + "\n"


def book_symbols(data):
    """The `symbols` map of a matrix file (`book_symbols`, or `book_epsilon` alone)."""
    if "book_symbols" in data:
        out = {}
        for k, v in data["book_symbols"].items():
            m = re.fullmatch(r"exp\((\d+) pi i / (\d+)\)", v)
            assert m, v
            out[k] = (int(m.group(1)) // 2, int(m.group(2)))
        return out
    return {"e": (1, _int_after(r"/ (\d+)", data["book_epsilon"]))}


# book element label -> canonical label, for the D_n tables without a character file
DN_ELEMENTS = {
    "C5+": "C5^1",
    "C8+": "C8^1",
    "C9+": "C9^1",
    "C10+": "C10^1",
    "C2'1": "C2'(1)",
}


def check_matrices():
    """Unitarity, D(g)^order = +-1, and trace = the character table's entry, per generator.

    The trace is compared with the transcribed character table when one
    exists, otherwise (the D_n tables of ticket #22) with the computed
    character of the canonical irrep of the same symbol.
    """
    bad = 0
    for group in MATRICES:
        path = MAT_DIR / f"{ag.ALIASES.get(group, group)}.toml"
        data = tomllib.loads(path.read_text())
        symbols = book_symbols(data)
        char_path = CHAR_DIR / path.name
        if char_path.exists():
            chars = tomllib.loads(char_path.read_text())
            F = cb.field(math.lcm(data["conductor"], chars["conductor"]))
            csymbols = {
                letter: parse_symbol_text(chars[f"book_{name}"])
                for letter, name in SYMBOL_NAMES.items()
                if f"book_{name}" in chars
            }
            classes = [c["label"] for c in chars["classes"]]
            rows = {r["label"]: r for r in chars["irreps"]}

            # bound as defaults: the closure is built once per loop iteration
            def expected(label, g, classes=classes, F=F, rows=rows, csymbols=csymbols):
                cl = element_class(g, classes)
                _, order, _ = class_invariants(cl)
                want = token_value(
                    F,
                    rows[label]["book_values"][classes.index(cl)],
                    None,
                    None,
                    csymbols,
                )
                return order, want, cl

        else:
            comp = Computed(ag.ALIASES.get(group, group))
            F = cb.field(math.lcm(data["conductor"], comp.F.conductor))
            by_symbol = {ir.symbol: ir for ir in comp.irreps}

            def expected(label, g, comp=comp, F=F, by_symbol=by_symbol):
                idx = comp.ctx.labels_list.index(DN_ELEMENTS[g])
                order = ag.element_order(comp.ctx.G, comp.ctx.G.elements[idx])
                want = lift_to(comp.F, F, by_symbol[label].chars[idx])
                return order, want, DN_ELEMENTS[g] + " (computed)"

        problems = []
        for label, entry in data["irreps"].items():
            for g, A, Abook in zip(
                data["elements"], entry["matrices"], entry["book_values"]
            ):
                D = [[cb.parse_value(F, x) for x in row] for row in A]
                assert cb.meq(
                    D,
                    [[token_value(F, t, 1, 1, symbols) for t in row] for row in Abook],
                ), (path.name, label, g)
                if not cb.is_unitary(F, D):
                    problems.append(f"{label} {g}: not unitary")
                order, want, cl = expected(label, g)
                P = cb.midentity(F, len(D))
                for _ in range(order):
                    P = cb.mmul(F, P, D)
                sign = (
                    F.rational(-1)
                    if ("/" in label and P[0][0] == F.rational(-1))
                    else F.rational(1)
                )
                if not cb.meq(P, cb.mscale(F, cb.midentity(F, len(D)), sign)):
                    problems.append(f"{label} {g}: D^{order} is not +-1")
                if cb.mtrace(F, D) != want:
                    problems.append(
                        f"{label} {g}: trace differs from the character table entry on {cl}"
                    )
        print(f"{path.name}: {'ok' if not problems else 'PROBLEMS'}")
        for pr in problems:
            print("  " + pr)
        bad += bool(problems)
    return bad


def load_book(ocr_dir):
    """All transcriptions: typeset pages parsed from OCR_DIR, scanned pages from BOOK."""
    out = {}
    for group, (table, page, typeset) in PAGES.items():
        if typeset:
            if not ocr_dir:
                continue
            text = (Path(ocr_dir) / f"p{page:03d}.txt").read_text()
            book = parse_typeset(group, text)
            if group in CONTINUED:
                more = (Path(ocr_dir) / f"p{CONTINUED[group]:03d}.txt").read_text()
                book = merge_continued(book, parse_typeset(group, more))
            out[group] = (book, True)
        elif group in BOOK:
            out[group] = (BOOK[group], False)
    return out


def generate(ocr_dir):
    CHAR_DIR.mkdir(parents=True, exist_ok=True)
    books = load_book(ocr_dir)
    for group, (book, typeset) in books.items():
        symbol = ag.ALIASES.get(group, group)
        (CHAR_DIR / f"{symbol}.toml").write_text(write_characters(group, book, typeset))
    print(f"wrote {len(books)} character files to {CHAR_DIR}")
    MAT_DIR.mkdir(parents=True, exist_ok=True)
    for group, mats in MATRICES.items():
        (MAT_DIR / f"{ag.ALIASES.get(group, group)}.toml").write_text(
            write_matrices(group, mats)
        )
    print(f"wrote {len(MATRICES)} matrix files to {MAT_DIR}")


def read_file(path):
    data = tomllib.loads(path.read_text())
    book = tbl(
        classes=" ".join(c["label"] for c in data["classes"]),
        symbols={
            letter: parse_symbol_text(data[f"book_{name}"])
            for letter, name in SYMBOL_NAMES.items()
            if f"book_{name}" in data
        }
        or None,
        rows=[
            (r["label"], " ".join(r["book_values"]), r["tau"]) for r in data["irreps"]
        ],
    )
    return data["group"], book, data


def apply_discrepancy(book, disc):
    rows = book["rows"]
    if "swap_projective" in disc:
        classes = book["classes"].split()
        a, b = disc["swap_projective"]
        ia, ib = classes.index(a), classes.index(b)
        swapped = []
        for lab, vals, tau in rows:
            toks = vals.split()
            if "/" in lab:
                toks[ia], toks[ib] = toks[ib], toks[ia]
            swapped.append((lab, " ".join(toks), tau))
        rows = swapped
    if "swap_labels" in disc:
        a, b = disc["swap_labels"]
        renamed = {a: b, b: a}
        rows = [(renamed.get(lab, lab), vals, tau) for lab, vals, tau in rows]
    return {**book, "rows": rows}


def check(groups):
    paths = sorted(CHAR_DIR.glob("*.toml"))
    if groups:
        wanted = {ag.ALIASES.get(g, g) for g in groups}
        paths = [p for p in paths if p.stem in wanted]
    bad = 0
    for path in paths:
        symbol, book, data = read_file(path)
        comp = Computed(symbol)
        # the grammar strings must parse back to the token values
        M = data["conductor"]
        F = cb.field(M)
        symbols = book_symbols_of(book)
        for r in data["irreps"]:
            for s, t in zip(r["characters"], r["book_values"]):
                assert cb.parse_value(F, s) == token_value(F, t, None, None, symbols), (
                    symbol,
                    r["label"],
                    s,
                    t,
                )
        lines, ok = structural_match(symbol, book, comp)
        if not ok and "known_discrepancy" in data:
            fixed = book
            for disc in data["known_discrepancy"]:
                fixed = apply_discrepancy(fixed, disc)
            lines2, ok2 = structural_match(symbol, fixed, comp)
            lines = [
                "  documented discrepancy: "
                + data["known_discrepancy"][0]["note"][:120]
                + "...",
                f"  after the documented correction: {'match' if ok2 else 'STILL NO MATCH'}",
            ] + lines2
            ok = ok2
        print(f"{path.name}: {'ok' if ok else 'MISMATCH'}")
        for ln in lines:
            print(ln)
        bad += not ok
    print(f"{len(paths)} file(s), {bad} with mismatches")
    bad += check_matrices()
    return 1 if bad else 0


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    if argv[0] == "generate":
        generate(argv[1] if len(argv) > 1 else None)
        return 0
    if argv[0] == "check":
        return check(argv[1:])
    if argv[0] == "parse":  # debug: show the parsed typeset table of one group
        group, ocr_dir = argv[1], argv[2]
        _, page, _ = PAGES[group]
        print(parse_typeset(group, (Path(ocr_dir) / f"p{page:03d}.txt").read_text()))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
