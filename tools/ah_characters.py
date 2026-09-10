"""Transcription of the Altmann & Herzig character tables and sample matrices
(wayfinder ticket #19).

Two sources feed one set of reference files:

* the **typeset** character tables of the master groups and of the groups that
  own a page of their own (14 groups) are parsed from the OCR text of the
  PDF (pymupdf ``get_text`` of the table T n.4);
* the **scanned** group pages (25 groups) are read from the page images and
  typed into ``BOOK`` below, in the book's own notation transliterated to
  ASCII (``e`` = epsilon, ``d`` = delta, ``2c5^2`` = 2 cos(2 pi / 5), ``sqrt2``,
  ``ie*``); the generator matrices of six axial groups (T n.7) into
  ``MATRICES`` the same way.

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

import itertools
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
}  # fmt: skip

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

# Documented disagreements between the book's table and its own parameter tables.
# swap_projective: two class columns whose projective-row entries the book prints interchanged.
DISCREPANCIES: dict[str, dict[str, Any]] = {
    "D5h": tbl(
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
}

# Generator matrices of T n.7 for the axial samples, read from the page images.
# elements: book labels of the generators; irreps: book label -> one matrix per element,
# entries in the token grammar above (the book's overbar is a minus sign).
MATRICES: dict[str, dict[str, Any]] = {
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

# --------------------------------------------------------------------------
# Values: book tokens -> exact field elements
# --------------------------------------------------------------------------

TOKEN_RE = re.compile(r"^(-)?(\d+)?(i)?(e\*|e|d\*|d|sqrt2|sqrt3|c(\d+)(?:\^(\d+))?)?$")


def _int_after(pattern, text):
    m = re.search(pattern, text)
    assert m, (pattern, text)
    return int(m.group(1))


def token_conductor(tok, eps, delta):
    m = TOKEN_RE.match(tok)
    assert m, tok
    need = 1
    if m.group(3):
        need = math.lcm(need, 4)
    base = m.group(4)
    if base in ("e", "e*"):
        need = math.lcm(need, eps)
    elif base in ("d", "d*"):
        need = math.lcm(need, delta)
    elif base == "sqrt2":
        need = math.lcm(need, 8)
    elif base == "sqrt3":
        need = math.lcm(need, 12)
    elif base and base.startswith("c"):
        need = math.lcm(need, 2 * int(m.group(5)))
    return need


def token_value(F, tok, eps, delta):
    m = TOKEN_RE.match(tok)
    assert m, tok
    sign, coeff, has_i, base = m.group(1), m.group(2), m.group(3), m.group(4)
    x = F.rational(int(coeff) if coeff else 1)
    if base in ("e", "e*", "d", "d*"):
        n = eps if base[0] == "e" else delta
        z = (
            cb.lift(cb.field(n), F, cb.field(n).zeta(1))
            if F.conductor % n == 0
            else None
        )
        assert z is not None, (tok, F.conductor)
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

GREEK = {"ϵ": "e", "δ": "d", "σ": "sigma", "′": "'", "″": "''", "−": "-", "∗": "*"}


def normalise(text):
    for a, b in GREEK.items():
        text = text.replace(a, b)
    return text


def parse_typeset(group, text):
    """Parse the table T n.4 of a typeset page into the BOOK format."""
    text = normalise(text)
    seg = text[text.index("Character table") :]
    seg = seg.split("Cartesian tensors")[0]
    toks = seg.replace("\n", " ").split()
    # constants
    eps = delta = None
    for k, t in enumerate(toks):
        if t == "e" and toks[k + 1] == "=" and toks[k + 2].startswith("exp("):
            eps = _int_after(r"/(\d+)\)", toks[k + 2])
        if t == "d" and toks[k + 1] == "=" and toks[k + 2].startswith("exp("):
            delta = _int_after(r"/(\d+)\)", toks[k + 2])
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
        if t in ("e", "d") and k + 1 < len(body) and body[k + 1] == "=":
            break
        if t.startswith("cm") or t.startswith("c") and t.endswith("m"):
            break
        if t == "T" and re.match(r"^\d+\.\d$", body[k + 1]):
            break
        label, k = merge_label_tokens(body, k)
        vals = []
        while len(vals) < n:
            v, k = merge_value_tokens(body, k)
            vals.append(v)
        tau_v = body[k]
        assert tau_v in ("a", "b", "c"), (group, label, tau_v, vals)
        k += 1
        rows.append((label, " ".join(vals), tau_v))
    return tbl(classes=" ".join(classes), eps=eps, delta=delta, rows=rows)


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
        nxt = body[k + 1] if k + 1 < len(body) else ""
        if (
            nxt.isdigit()
            and int(nxt) in (5, 8, 10, 12, 16, 20)
            and int(m.group(3)) <= 4
        ):
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
    Fb_cond = 1
    eps, delta = book.get("eps") or 1, book.get("delta") or 1
    classes = book["classes"].split()
    for _, vals, _ in book["rows"]:
        for tok in vals.split():
            Fb_cond = math.lcm(Fb_cond, token_conductor(tok, eps, delta))
    M = math.lcm(Fb_cond, comp.F.conductor)
    F = cb.field(M)
    book_vals = [
        [
            lift_to(
                cb.field(Fb_cond), F, token_value(cb.field(Fb_cond), tok, eps, delta)
            )
            if False
            else token_value(F, tok, eps, delta)
            for tok in vals.split()
        ]
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
    regular = [reg for *_, reg in comp.classes]
    # candidate computed classes per book column
    cands = [[j for j in range(n) if inv_comp[j] == inv_book[i]] for i in range(n)]
    order_ = sorted(range(n), key=lambda i: len(cands[i]))
    proj_book = [("/" in lab) for lab, _, _ in book["rows"]]
    proj_comp = [ir.projective for ir in comp.irreps]

    def try_rows(perm, signs):
        # perm[i] = computed class for book column i; signs per computed class
        used = set()
        pairing = []
        for r, vals in enumerate(book_vals):
            target = None
            for s, row in enumerate(comp_rows):
                if s in used or proj_comp[s] != proj_book[r]:
                    continue
                ok = True
                for i in range(n):
                    j = perm[i]
                    v = row[j]
                    if proj_book[r] and signs[j] < 0:
                        v = F.neg(v)
                    if v != vals[i]:
                        ok = False
                        break
                if ok:
                    target = s
                    break
            if target is None:
                return None
            used.add(target)
            pairing.append(target)
        return pairing

    def search(pos, perm, used):
        if pos == n:
            reg_cols = [j for j in range(n) if regular[j]]
            for bits in itertools.product((1, -1), repeat=len(reg_cols)):
                signs = [1] * n
                for j, b in zip(reg_cols, bits):
                    signs[j] = b
                pairing = try_rows(perm, signs)
                if pairing is not None:
                    return perm[:], signs, pairing
            return None
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
    eps, delta = book.get("eps") or 1, book.get("delta") or 1
    M = 1
    for _, vals, _ in book["rows"]:
        for tok in vals.split():
            M = math.lcm(M, token_conductor(tok, eps, delta))
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
    if book.get("eps"):
        lines.append(f'book_epsilon = "exp(2 pi i / {eps})"')
    if book.get("delta"):
        lines.append(f'book_delta = "exp(2 pi i / {delta})"')
    if group in DISCREPANCIES:
        d = DISCREPANCIES[group]
        lines.append("")
        lines.append("[[known_discrepancy]]")
        lines.append(
            f"swap_projective = [{toml_str(d['swap_projective'][0])}, {toml_str(d['swap_projective'][1])}]"
        )
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
        exact = [cb.to_grammar(F, token_value(F, t, eps, delta)) for t in toks]
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
    root = re.sub(r"(sigma_[vd]|C\d+'*|S\d+(?:\^\d+)?)\d$", r"\1", root)
    for c in classes:
        if re.sub(r"^\d+", "", c) == root:
            return c
    raise KeyError((label, classes))


def write_matrices(group, mats):
    symbol = ag.ALIASES.get(group, group)
    eps = mats["eps"]
    M = 1
    for irrep in mats["irreps"].values():
        for A in irrep:
            for row in A:
                for tok in row:
                    M = math.lcm(M, token_conductor(tok, eps, 1))
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
        f'book_epsilon = "exp(2 pi i / {eps})"',
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
                    toml_str(cb.to_grammar(F, token_value(F, t, eps, 1))) for t in row
                )
                + "]"
                for row in A
            ]
            lines.append("  [" + ", ".join(rows) + "],")
        lines.append("]")
        lines.append("")
    return "\n".join(lines) + "\n"


def check_matrices():
    """Unitarity, D(g)^order = +-1, and trace = the character file's entry, per generator."""
    bad = 0
    for group in MATRICES:
        path = MAT_DIR / f"{ag.ALIASES.get(group, group)}.toml"
        data = tomllib.loads(path.read_text())
        chars = tomllib.loads((CHAR_DIR / path.name).read_text())
        F = cb.field(math.lcm(data["conductor"], chars["conductor"]))
        eps = _int_after(r"/ (\d+)", data["book_epsilon"])
        ceps = (
            _int_after(r"/ (\d+)", chars["book_epsilon"])
            if "book_epsilon" in chars
            else 1
        )
        classes = [c["label"] for c in chars["classes"]]
        rows = {r["label"]: r for r in chars["irreps"]}
        problems = []
        for label, entry in data["irreps"].items():
            for g, A, Abook in zip(
                data["elements"], entry["matrices"], entry["book_values"]
            ):
                D = [[cb.parse_value(F, x) for x in row] for row in A]
                assert cb.meq(
                    D, [[token_value(F, t, eps, 1) for t in row] for row in Abook]
                ), (path.name, label, g)
                if not cb.is_unitary(F, D):
                    problems.append(f"{label} {g}: not unitary")
                cl = element_class(g, classes)
                _, order, _ = class_invariants(cl)
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
                want = token_value(
                    F, rows[label]["book_values"][classes.index(cl)], ceps, 1
                )
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
            out[group] = (parse_typeset(group, text), True)
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
        eps=_int_after(r"/ (\d+)", data["book_epsilon"])
        if "book_epsilon" in data
        else None,
        delta=_int_after(r"/ (\d+)", data["book_delta"])
        if "book_delta" in data
        else None,
        rows=[
            (r["label"], " ".join(r["book_values"]), r["tau"]) for r in data["irreps"]
        ],
    )
    return data["group"], book, data


def apply_discrepancy(book, disc):
    classes = book["classes"].split()
    a, b = disc["swap_projective"]
    ia, ib = classes.index(a), classes.index(b)
    rows = []
    for lab, vals, tau in book["rows"]:
        toks = vals.split()
        if "/" in lab:
            toks[ia], toks[ib] = toks[ib], toks[ia]
        rows.append((lab, " ".join(toks), tau))
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
        eps, delta = book["eps"] or 1, book["delta"] or 1
        for r in data["irreps"]:
            for s, t in zip(r["characters"], r["book_values"]):
                assert cb.parse_value(F, s) == token_value(F, t, eps, delta), (
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
