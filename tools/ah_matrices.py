"""Transcription of the irrep matrices of Altmann & Herzig, *Point-Group Theory
Tables* (2nd ed.), for the polyhedral groups (wayfinder ticket #17).

The book prints the matrices of every operation of ``O`` (T 69.7, pp. 582-584)
and ``T`` (T 70.7, p. 592) in every irrep of dimension >= 2, ordinary and
spinor, and for ``I`` the matrices of three generators (T 74.7b, p. 655)
together with a word in those generators for every operation (T 74.7a,
p. 654).  ``Td`` reads T 69.7 with its own operation names (T 73.7); ``Th``,
``Oh`` and ``Ih`` are direct products with ``Ci`` and reuse the tables of
``T``, ``O``, ``I`` with ``D(i g) = +-D(g)`` (book section 16-7, p. 78).

Every printed entry is transcribed below in the book's own symbols (``d``
for delta, ``e`` for epsilon, ``n`` for eta, ``w`` for omega, a trailing
``s`` for the complex conjugate, a leading ``-`` for the book's overline),
and evaluated exactly in ``Q(zeta_M)`` with the value grammar of
``docs/spec/basis-conventions.md``.  The irrational constants of T 74.7b,
printed as ``r exp(i arctan x)``, are identified as exact surds
(``IDENTIFIED`` below) and re-checked against the printed formula in floating
point.

    uv run python tools/ah_matrices.py check       # verify the transcription
    uv run python tools/ah_matrices.py generate    # write tests/reference/altmann-herzig/matrices/

Checks (all exact in the field): every matrix is unitary; the matrices of the
book's untilded operations (for ``I``: the products along the book's words)
multiply as ``D(a) D(b) = sign * D(c)`` with the sign of the quaternion
product, i.e. ``-1`` exactly where the book prints a tilde for a spinor
irrep; the trace of every matrix equals the character of exactly one
canonical irrep of the same dimension (identifying the book's symbol); the
words of T 74.7a reproduce the parameters of T 75.1.
"""

from __future__ import annotations

import cmath
import math
import sys
import tomllib
from fractions import Fraction as Fr
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
sys.path.insert(0, str(ROOT / "tools"))
import ah_book as ab
import canonical_basis_prototype as cb

SOURCE = ab.SOURCE
OUT_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "matrices"
LABEL_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "labels"

# --------------------------------------------------------------------------
# The transcription
# --------------------------------------------------------------------------
# A matrix is "prefactor | row; row; ..." (prefactor optional), entries separated
# by blanks (or by commas where an entry contains a binary + or -), in the value
# grammar extended by the group's symbols.

# T 69.7 (O, pp. 582-584): delta = exp(2 pi i/24), epsilon = exp(2 pi i/8),
# eta = exp(2 pi i/3).  Conductor 24.
O_SYMBOLS = {
    "d": "zeta",
    "e": "zeta^3",
    "n": "zeta^8",
}
O_IRREPS = ["E", "T1", "T2", "E1/2", "E5/2", "F3/2"]
O_TABLE = {
    "E": [
        "1 0; 0 1",
        "1 0 0; 0 1 0; 0 0 1",
        "1 0 0; 0 1 0; 0 0 1",
        "1 0; 0 1",
        "1 0; 0 1",
        "1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1",
    ],
    "C2x": [
        "1 0; 0 1",
        "-1 0 0; 0 -1 0; 0 0 1",
        "-1 0 0; 0 -1 0; 0 0 1",
        "0 -i; -i 0",
        "0 -i; -i 0",
        "0 -i 0 0; -i 0 0 0; 0 0 0 -i; 0 0 -i 0",
    ],
    "C2y": [
        "1 0; 0 1",
        "1 0 0; 0 -1 0; 0 0 -1",
        "1 0 0; 0 -1 0; 0 0 -1",
        "0 -1; 1 0",
        "0 -1; 1 0",
        "0 -1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 1 0",
    ],
    "C2z": [
        "1 0; 0 1",
        "-1 0 0; 0 1 0; 0 0 -1",
        "-1 0 0; 0 1 0; 0 0 -1",
        "-i 0; 0 i",
        "-i 0; 0 i",
        "-i 0 0 0; 0 i 0 0; 0 0 -i 0; 0 0 0 i",
    ],
    "C3,1+": [
        "n 0; 0 ns",
        "0 0 -i; -i 0 0; 0 -1 0",
        "0 0 -i; -i 0 0; 0 -1 0",
        "sqrt2/2 | es -e; es e",
        "sqrt2/2 | es -e; es e",
        "sqrt2/2 | i*ds ds 0 0; i*ds -ds 0 0; 0 0 -d i*d; 0 0 -d -i*d",
    ],
    "C3,2+": [
        "n 0; 0 ns",
        "0 0 -i; i 0 0; 0 1 0",
        "0 0 -i; i 0 0; 0 1 0",
        "sqrt2/2 | es e; -es e",
        "sqrt2/2 | es e; -es e",
        "sqrt2/2 | i*ds -ds 0 0; -i*ds -ds 0 0; 0 0 -d -i*d; 0 0 d -i*d",
    ],
    "C3,3+": [
        "n 0; 0 ns",
        "0 0 i; -i 0 0; 0 1 0",
        "0 0 i; -i 0 0; 0 1 0",
        "sqrt2/2 | e es; -e es",
        "sqrt2/2 | e es; -e es",
        "sqrt2/2 | -ds i*ds 0 0; ds i*ds 0 0; 0 0 -i*d -d; 0 0 i*d -d",
    ],
    "C3,4+": [
        "n 0; 0 ns",
        "0 0 i; i 0 0; 0 -1 0",
        "0 0 i; i 0 0; 0 -1 0",
        "sqrt2/2 | e -es; e es",
        "sqrt2/2 | e -es; e es",
        "sqrt2/2 | -ds -i*ds 0 0; -ds i*ds 0 0; 0 0 -i*d d; 0 0 -i*d -d",
    ],
    "C3,1-": [
        "ns 0; 0 n",
        "0 i 0; 0 0 -1; i 0 0",
        "0 i 0; 0 0 -1; i 0 0",
        "sqrt2/2 | e e; -es es",
        "sqrt2/2 | e e; -es es",
        "sqrt2/2 | -i*d -i*d 0 0; d -d 0 0; 0 0 -ds -ds; 0 0 -i*ds i*ds",
    ],
    "C3,2-": [
        "ns 0; 0 n",
        "0 -i 0; 0 0 1; i 0 0",
        "0 -i 0; 0 0 1; i 0 0",
        "sqrt2/2 | e -e; es es",
        "sqrt2/2 | e -e; es es",
        "sqrt2/2 | -i*d i*d 0 0; -d -d 0 0; 0 0 -ds ds; 0 0 i*ds i*ds",
    ],
    "C3,3-": [
        "ns 0; 0 n",
        "0 i 0; 0 0 1; -i 0 0",
        "0 i 0; 0 0 1; -i 0 0",
        "sqrt2/2 | es -es; e e",
        "sqrt2/2 | es -es; e e",
        "sqrt2/2 | -d d 0 0; -i*d -i*d 0 0; 0 0 i*ds -i*ds; 0 0 -ds -ds",
    ],
    "C3,4-": [
        "ns 0; 0 n",
        "0 -i 0; 0 0 -1; -i 0 0",
        "0 -i 0; 0 0 -1; -i 0 0",
        "sqrt2/2 | es es; -e e",
        "sqrt2/2 | es es; -e e",
        "sqrt2/2 | -d -d 0 0; i*d -i*d 0 0; 0 0 i*ds i*ds; 0 0 ds -ds",
    ],
    "C4x+": [
        "0 ns; n 0",
        "0 -i 0; -i 0 0; 0 0 1",
        "0 i 0; i 0 0; 0 0 -1",
        "sqrt2/2 | 1 -i; -i 1",
        "sqrt2/2 | -1 i; i -1",
        "sqrt2/2 | 0 0 ns -i*ns; 0 0 -i*ns ns; n -i*n 0 0; -i*n n 0 0",
    ],
    "C4y+": [
        "0 n; ns 0",
        "1 0 0; 0 0 1; 0 -1 0",
        "-1 0 0; 0 0 -1; 0 1 0",
        "sqrt2/2 | 1 -1; 1 1",
        "sqrt2/2 | -1 1; -1 -1",
        "sqrt2/2 | 0 0 n -n; 0 0 n n; ns -ns 0 0; ns ns 0 0",
    ],
    "C4z+": [
        "0 1; 1 0",
        "0 0 -i; 0 1 0; -i 0 0",
        "0 0 i; 0 -1 0; i 0 0",
        "es 0; 0 e",
        "-es 0; 0 -e",
        "0 0 es 0; 0 0 0 e; es 0 0 0; 0 e 0 0",
    ],
    "C4x-": [
        "0 ns; n 0",
        "0 i 0; i 0 0; 0 0 1",
        "0 -i 0; -i 0 0; 0 0 -1",
        "sqrt2/2 | 1 i; i 1",
        "sqrt2/2 | -1 -i; -i -1",
        "sqrt2/2 | 0 0 ns i*ns; 0 0 i*ns ns; n i*n 0 0; i*n n 0 0",
    ],
    "C4y-": [
        "0 n; ns 0",
        "1 0 0; 0 0 -1; 0 1 0",
        "-1 0 0; 0 0 1; 0 -1 0",
        "sqrt2/2 | 1 1; -1 1",
        "sqrt2/2 | -1 -1; 1 -1",
        "sqrt2/2 | 0 0 n n; 0 0 -n n; ns ns 0 0; -ns ns 0 0",
    ],
    "C4z-": [
        "0 1; 1 0",
        "0 0 i; 0 1 0; i 0 0",
        "0 0 -i; 0 -1 0; -i 0 0",
        "e 0; 0 es",
        "-e 0; 0 -es",
        "0 0 e 0; 0 0 0 es; e 0 0 0; 0 es 0 0",
    ],
    "C2'a": [
        "0 1; 1 0",
        "0 0 -i; 0 -1 0; i 0 0",
        "0 0 i; 0 1 0; -i 0 0",
        "0 -e; es 0",
        "0 e; -es 0",
        "0 0 0 -e; 0 0 es 0; 0 -e 0 0; es 0 0 0",
    ],
    "C2'b": [
        "0 1; 1 0",
        "0 0 i; 0 -1 0; -i 0 0",
        "0 0 -i; 0 1 0; i 0 0",
        "0 -es; e 0",
        "0 es; -e 0",
        "0 0 0 -es; 0 0 e 0; 0 -es 0 0; e 0 0 0",
    ],
    "C2'c": [
        "0 n; ns 0",
        "-1 0 0; 0 0 -1; 0 -1 0",
        "1 0 0; 0 0 1; 0 1 0",
        "sqrt2/2 | -i -i; -i i",
        "sqrt2/2 | i i; i -i",
        "sqrt2/2 | 0 0 -i*n -i*n; 0 0 -i*n i*n; -i*ns -i*ns 0 0; -i*ns i*ns 0 0",
    ],
    "C2'd": [
        "0 ns; n 0",
        "0 i 0; -i 0 0; 0 0 -1",
        "0 -i 0; i 0 0; 0 0 1",
        "sqrt2/2 | i 1; -1 -i",
        "sqrt2/2 | -i -1; 1 i",
        "sqrt2/2 | 0 0 i*ns ns; 0 0 -ns -i*ns; i*n n 0 0; -n -i*n 0 0",
    ],
    "C2'e": [
        "0 n; ns 0",
        "-1 0 0; 0 0 1; 0 1 0",
        "1 0 0; 0 0 -1; 0 -1 0",
        "sqrt2/2 | i -i; -i -i",
        "sqrt2/2 | -i i; i i",
        "sqrt2/2 | 0 0 i*n -i*n; 0 0 -i*n -i*n; i*ns -i*ns 0 0; -i*ns -i*ns 0 0",
    ],
    "C2'f": [
        "0 ns; n 0",
        "0 -i 0; i 0 0; 0 0 -1",
        "0 i 0; -i 0 0; 0 0 1",
        "sqrt2/2 | -i 1; -1 i",
        "sqrt2/2 | i -1; 1 -i",
        "sqrt2/2 | 0 0 -i*ns ns; 0 0 -ns i*ns; -i*n n 0 0; -n i*n 0 0",
    ],
}

# T 70.7 (T, p. 592): delta = 2^(-1/2) exp(2 pi i/24), epsilon = 2^(-1/2) exp(2 pi i/8)
# (the 1/sqrt2 is inside the symbols here).  Conductor 24.
T_SYMBOLS = {
    "d": "sqrt2/2*zeta",
    "e": "sqrt2/2*zeta^3",
}
T_IRREPS = ["T", "E1/2", "1F3/2", "2F3/2"]
T_TABLE = {
    "E": ["1 0 0; 0 1 0; 0 0 1", "1 0; 0 1", "1 0; 0 1", "1 0; 0 1"],
    "C2x": ["-1 0 0; 0 -1 0; 0 0 1", "0 -i; -i 0", "0 -i; -i 0", "0 -i; -i 0"],
    "C2y": ["1 0 0; 0 -1 0; 0 0 -1", "0 -1; 1 0", "0 -1; 1 0", "0 -1; 1 0"],
    "C2z": ["-1 0 0; 0 1 0; 0 0 -1", "-i 0; 0 i", "-i 0; 0 i", "-i 0; 0 i"],
    "C3,1+": [
        "0 0 -i; -i 0 0; 0 -1 0",
        "es -e; es e",
        "i*ds ds; i*ds -ds",
        "-d i*d; -d -i*d",
    ],
    "C3,2+": [
        "0 0 -i; i 0 0; 0 1 0",
        "es e; -es e",
        "i*ds -ds; -i*ds -ds",
        "-d -i*d; d -i*d",
    ],
    "C3,3+": [
        "0 0 i; -i 0 0; 0 1 0",
        "e es; -e es",
        "-ds i*ds; ds i*ds",
        "-i*d -d; i*d -d",
    ],
    "C3,4+": [
        "0 0 i; i 0 0; 0 -1 0",
        "e -es; e es",
        "-ds -i*ds; -ds i*ds",
        "-i*d d; -i*d -d",
    ],
    "C3,1-": [
        "0 i 0; 0 0 -1; i 0 0",
        "e e; -es es",
        "-i*d -i*d; d -d",
        "-ds -ds; -i*ds i*ds",
    ],
    "C3,2-": [
        "0 -i 0; 0 0 1; i 0 0",
        "e -e; es es",
        "-i*d i*d; -d -d",
        "-ds ds; i*ds i*ds",
    ],
    "C3,3-": [
        "0 i 0; 0 0 1; -i 0 0",
        "es -es; e e",
        "-d d; -i*d -i*d",
        "i*ds -i*ds; -ds -ds",
    ],
    "C3,4-": [
        "0 -i 0; 0 0 -1; -i 0 0",
        "es es; -e e",
        "-d -d; i*d -i*d",
        "i*ds i*ds; ds -ds",
    ],
}

# T 74.7b (I, p. 655).  Conductor 120.  The book's definitions, verbatim:
#   delta = exp(2 pi i/24), eta = exp(2 pi i/8), omega = exp(2 pi i/3),
#   g+- = (sqrt5 +- 1)/2, s+- = (sqrt3 +- 1)/sqrt2, t = sqrt5,
#   epsilon = (8 + 3 sqrt3)^(1/2) exp[i arctan{(-6416 + 13216 sqrt3 + 2352 sqrt5 - 5075 sqrt15)/13189}],
#   zeta    = (8 - 3 sqrt3)^(1/2) exp[i arctan{(6416 + 13216 sqrt3 - 2352 sqrt5 - 5075 sqrt15)/13189}],
#   theta   = sqrt5 exp[i arctan{(4 - sqrt5)/sqrt3}],
#   lambda  = exp(i arctan sqrt(5/3)),
#   mu      = (4 + sqrt3)^(1/2) exp[i arctan{(84 - 16 sqrt3 + 72 sqrt5 + 5 sqrt15)/131}],
#   nu      = (4 - sqrt3)^(1/2) exp[i arctan{(84 + 16 sqrt3 + 72 sqrt5 - 5 sqrt15)/131}],
#   xi      = sqrt5 exp(i arctan 2),
#   rho     = sqrt7 exp[i arctan{(80 sqrt3 + 49 sqrt15)/57}],
#   sigma   = sqrt5 exp{i arctan(28 + 16 sqrt3 + 12 sqrt5 + 7 sqrt15)}.
# The arctan forms are identified (PSLQ over 1, sqrt2, sqrt3, sqrt5 and their
# products, 80 digits) as the surds below; `check` re-evaluates the printed
# formulas in floating point against them.
I_PRINTED = {  # name -> (modulus, argument) as printed, evaluated in floating point
    "eps": (
        math.sqrt(8 + 3 * math.sqrt(3)),
        math.atan(
            (-6416 + 13216 * math.sqrt(3) + 2352 * math.sqrt(5) - 5075 * math.sqrt(15))
            / 13189
        ),
    ),
    "zet": (
        math.sqrt(8 - 3 * math.sqrt(3)),
        math.atan(
            (6416 + 13216 * math.sqrt(3) - 2352 * math.sqrt(5) - 5075 * math.sqrt(15))
            / 13189
        ),
    ),
    "th": (math.sqrt(5), math.atan((4 - math.sqrt(5)) / math.sqrt(3))),
    "l": (1.0, math.atan(math.sqrt(5 / 3))),
    "mu": (
        math.sqrt(4 + math.sqrt(3)),
        math.atan(
            (84 - 16 * math.sqrt(3) + 72 * math.sqrt(5) + 5 * math.sqrt(15)) / 131
        ),
    ),
    "nu": (
        math.sqrt(4 - math.sqrt(3)),
        math.atan(
            (84 + 16 * math.sqrt(3) + 72 * math.sqrt(5) - 5 * math.sqrt(15)) / 131
        ),
    ),
    "xi": (math.sqrt(5), math.atan(2)),
    "rho": (
        math.sqrt(7),
        math.atan((80 * math.sqrt(3) + 49 * math.sqrt(15)) / 57),
    ),
    "sig": (
        math.sqrt(5),
        math.atan(28 + 16 * math.sqrt(3) + 12 * math.sqrt(5) + 7 * math.sqrt(15)),
    ),
}
I_SYMBOLS = {
    "d": "zeta^5",
    "eta": "zeta^15",
    "w": "zeta^40",
    "gp": "(sqrt5 + 1)/2",
    "gm": "(sqrt5 - 1)/2",
    "sp": "(sqrt3 + 1)*sqrt2/2",
    "sm": "(sqrt3 - 1)*sqrt2/2",
    "t": "sqrt5",
    # IDENTIFIED surds of the printed arctan forms
    "eps": "(28 + 7*sqrt3 + 6*sqrt5 + sqrt3*sqrt5)/16 + i*(7 + 14*sqrt3 - 3*sqrt5 - 4*sqrt3*sqrt5)/16",
    "zet": "(28 - 7*sqrt3 + 6*sqrt5 - sqrt3*sqrt5)/16 + i*(-7 + 14*sqrt3 + 3*sqrt5 - 4*sqrt3*sqrt5)/16",
    "th": "(5*sqrt3 + sqrt3*sqrt5)/8 + i*(15 - sqrt5)/8",
    "l": "sqrt2*sqrt3/4 + i*sqrt2*sqrt5/4",
    "mu": "(3*sqrt2 + 3*sqrt2*sqrt3 + sqrt2*sqrt5 - sqrt2*sqrt3*sqrt5)/8 + i*(3*sqrt2 - sqrt2*sqrt3 + 3*sqrt2*sqrt5 + sqrt2*sqrt3*sqrt5)/8",
    "nu": "(3*sqrt2 - 3*sqrt2*sqrt3 + sqrt2*sqrt5 + sqrt2*sqrt3*sqrt5)/8 + i*(3*sqrt2 + sqrt2*sqrt3 + 3*sqrt2*sqrt5 - sqrt2*sqrt3*sqrt5)/8",
    "xi": "1 + 2*i",
    "rho": "(-7*sqrt3 + 5*sqrt3*sqrt5)/16 + i*(35 + 3*sqrt5)/16",
    "sig": "(15*sqrt2 - 5*sqrt2*sqrt3 - sqrt2*sqrt5 - sqrt2*sqrt3*sqrt5)/16 + i*(15*sqrt2 + 5*sqrt2*sqrt3 - sqrt2*sqrt5 + sqrt2*sqrt3*sqrt5)/16",
}
I_IRREPS = ["T1", "T2", "F", "H", "E1/2", "E7/2", "F3/2", "I5/2"]
I_GENERATORS = ["C5,1+", "C3,1+", "C2a"]
I_TABLE = {
    "C5,1+": [
        "1/2 | gm -i -i*gp; -i gp -gm; -i*gp -gm 1",
        "1/2 | -gp -i i*gm; -i -gm gp; i*gm gp 1",
        "1/4 | -1 -t -i*t -i*t; -t -1 3*i -i; -i*t 3*i 1 1; -i*t -i 1 -3",
        "1/2 | 0 l^2*ws -l -i*l*ws -i*l*w; ls^2*w 0 -ls -i*ls*w -i*ls*ws; -ls -l 1 0 i; -i*ls*w -i*l*ws 0 -1 -1; -i*ls*ws -i*l*w i -1 0",
        "1/2 | gp - i, -i*gm; -i*gm, gp + i",
        "1/2 | -gm - i, i*gp; i*gp, -gm + i",
        "sqrt2/4 | eta -i*sm mus -ls; i*sp etas -ls nu; nus l eta i*sp; l mu -i*sm etas",
        "1/6 | -xis -i -i*sigs -sp*th -sig -sm*ths; -i -xi sm*th -sigs sp*ths i*sig; -sig sp*ths -sqrt2*eta i*sqrt2*sm -zets -rhos; -sm*ths i*sig -i*sqrt2*sp -sqrt2*etas -rhos -epss; -i*sigs sm*th -eps rho -sqrt2*eta -i*sqrt2*sp; -sp*th -sigs rho -zet i*sqrt2*sm -sqrt2*etas",
    ],
    "C3,1+": [
        "0 0 -i; -i 0 0; 0 -1 0",
        "0 0 -i; -i 0 0; 0 -1 0",
        "1 0 0 0; 0 0 0 -i; 0 -i 0 0; 0 0 -1 0",
        "w 0 0 0 0; 0 ws 0 0 0; 0 0 0 0 -i; 0 0 -i 0 0; 0 0 0 -1 0",
        "sqrt2/2 | etas -eta; etas eta",
        "sqrt2/2 | etas -eta; etas eta",
        "sqrt2/2 | d -i*d 0 0; d i*d 0 0; 0 0 -i*ds -ds; 0 0 -i*ds ds",
        "sqrt2/2 | etas -eta 0 0 0 0; etas eta 0 0 0 0; 0 0 d -i*d 0 0; 0 0 d i*d 0 0; 0 0 0 0 -i*ds -ds; 0 0 0 0 -i*ds ds",
    ],
    "C2a": [
        "-1 0 0; 0 -1 0; 0 0 1",
        "-1 0 0; 0 -1 0; 0 0 1",
        "1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1",
        "1 0 0 0 0; 0 1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0; 0 0 0 0 1",
        "0 -i; -i 0",
        "0 -i; -i 0",
        "0 -i 0 0; -i 0 0 0; 0 0 0 -i; 0 0 -i 0",
        "0 -i 0 0 0 0; -i 0 0 0 0 0; 0 0 0 -i 0 0; 0 0 -i 0 0 0; 0 0 0 0 0 -i; 0 0 0 0 -i 0",
    ],
}
# T 74.7a (I, p. 654): every operation as a word in the generators (left to right as printed).
I_WORDS = {
    "E": "C2a C2a",
    "C5,1+": "C5,1+",
    "C5,2+": "C5,1+ C5,1+ C2a C3,1+",
    "C5,3+": "C5,1+ C2a C3,1+ C5,1+",
    "C5,4+": "C3,1+ C5,1+ C2a",
    "C5,5+": "C5,1+ C2a C5,1+ C2a",
    "C5,6+": "C2a C5,1+ C2a C5,1+",
    "C5,1-": "C3,1+ C5,1+ C3,1+",
    "C5,2-": "C2a C5,1+ C2a",
    "C5,3-": "C5,1+ C3,1+ C3,1+",
    "C5,4-": "C2a C3,1+ C5,1+",
    "C5,5-": "C3,1+ C3,1+ C5,1+",
    "C5,6-": "C5,1+ C2a C3,1+",
    "C5,1^2+": "C5,1+ C5,1+",
    "C5,2^2+": "C3,1+ C2a C5,1+",
    "C5,3^2+": "C3,1+ C2a C3,1+ C5,1+",
    "C5,4^2+": "C5,1+ C3,1+ C2a",
    "C5,5^2+": "C2a C3,1+ C5,1+ C3,1+",
    "C5,6^2+": "C3,1+ C5,1+ C3,1+ C2a",
    "C5,1^2-": "C5,1+ C5,1+ C5,1+",
    "C5,2^2-": "C2a C5,1+ C5,1+ C2a",
    "C5,3^2-": "C3,1+ C3,1+ C2a C5,1+ C2a",
    "C5,4^2-": "C2a C5,1+ C3,1+",
    "C5,5^2-": "C5,1+ C2a",
    "C5,6^2-": "C2a C5,1+",
    "C3,1+": "C3,1+",
    "C3,2+": "C3,1+ C2a",
    "C3,3+": "C2a C3,1+ C2a",
    "C3,4+": "C2a C3,1+",
    "C3,5+": "C5,1+ C2a C5,1+",
    "C3,6+": "C3,1+ C2a C5,1+ C3,1+ C2a",
    "C3,7+": "C2a C3,1+ C2a C5,1+",
    "C3,8+": "C2a C5,1+ C2a C3,1+",
    "C3,9+": "C3,1+ C5,1+ C2a C3,1+",
    "C3,10+": "C3,1+ C2a C5,1+ C2a",
    "C3,1-": "C3,1+ C3,1+",
    "C3,2-": "C2a C3,1+ C3,1+",
    "C3,3-": "C3,1+ C2a C3,1+",
    "C3,4-": "C3,1+ C3,1+ C2a",
    "C3,5-": "C5,1+ C2a C3,1+ C2a",
    "C3,6-": "C3,1+ C5,1+ C2a C5,1+",
    "C3,7-": "C2a C5,1+ C5,1+",
    "C3,8-": "C3,1+ C5,1+ C5,1+",
    "C3,9-": "C5,1+ C5,1+ C3,1+",
    "C3,10-": "C5,1+ C5,1+ C2a",
    "C2a": "C2a",
    "C2b": "C3,1+ C2a C3,1+ C3,1+",
    "C2c": "C3,1+ C3,1+ C2a C3,1+",
    "C2d": "C3,1+ C5,1+",
    "C2e": "C3,1+ C2a C5,1+ C3,1+",
    "C2f": "C5,1+ C3,1+ C2a C5,1+",
    "C2g": "C5,1+ C3,1+",
    "C2h": "C5,1+ C5,1+ C3,1+ C2a",
    "C2i": "C2a C5,1+ C3,1+ C3,1+",
    "C2j": "C2a C5,1+ C3,1+ C2a",
    "C2k": "C5,1+ C3,1+ C3,1+ C2a",
    "C2l": "C2a C5,1+ C5,1+ C3,1+",
    "C2m": "C2a C3,1+ C5,1+ C2a",
    "C2n": "C3,1+ C3,1+ C2a C5,1+",
    "C2o": "C3,1+ C5,1+ C2a C5,1+ C2a",
}

# Errata of T 74.7b, found by the checks below and applied before any matrix is
# used or written: with delta = exp(2 pi i/24) as printed, the C3,1+ matrices of
# F3/2 and I5/2 cube to +1 on their delta blocks (a spinor C3 must cube to -1,
# and T 70.7 prints the same blocks with the opposite sign), and the C5,1+
# matrices of F3/2 and I5/2 are printed with rows and columns interchanged.
# Reading delta as -delta in those C3,1+ matrices and transposing those C5,1+
# matrices makes every product along the words of T 74.7a agree with the
# book's factor system and every trace with the character; no other single
# sign, conjugation or transposition does (tools/ah_matrices.py history).
ERRATA_NOTE = {
    "negate d": (
        "With delta = exp(2 pi i/24) as printed the matrix cubes to +1 on its delta blocks; "
        "a spinor C3 must cube to -1, and T 70.7 prints the same blocks with the opposite "
        "sign. Read as -delta; with the transposed C5,1+ every product along the words of "
        "T 74.7a then agrees with the factor system and every trace with the character."
    ),
    "transpose": (
        "Printed with rows and columns interchanged: as printed no single sign, conjugation "
        "or transposition of any other entry makes the words of T 74.7a close; the "
        "transpose, with the corrected C3,1+, reproduces the factor system and the character."
    ),
}
ERRATA = {
    "I": [
        ("C3,1+", "F3/2", "negate d"),
        ("C3,1+", "I5/2", "negate d"),
        ("C5,1+", "F3/2", "transpose"),
        ("C5,1+", "I5/2", "transpose"),
    ],
}

GROUPS: dict[str, dict[str, Any]] = {
    "O": {
        "kind": "O",
        "conductor": 24,
        "table": "T 69.7, pp. 582-584",
        "symbols": O_SYMBOLS,
        "irreps": O_IRREPS,
        "matrices": O_TABLE,
        "generators": ["C4z+", "C2'd"],
        "words": None,
        "printed": {
            "delta": "exp(2 pi i/24)",
            "epsilon": "exp(2 pi i/8)",
            "eta": "exp(2 pi i/3)",
        },
    },
    "T": {
        "kind": "T",
        "conductor": 24,
        "table": "T 70.7, p. 592",
        "symbols": T_SYMBOLS,
        "irreps": T_IRREPS,
        "matrices": T_TABLE,
        "generators": ["C3,1+", "C2x"],
        "words": None,
        "printed": {
            "delta": "2^(-1/2) exp(2 pi i/24)",
            "epsilon": "2^(-1/2) exp(2 pi i/8)",
        },
    },
    "I": {
        "kind": "I",
        "conductor": 120,
        "table": "T 74.7a, p. 654 (words) and T 74.7b, p. 655 (generator matrices)",
        "symbols": I_SYMBOLS,
        "irreps": I_IRREPS,
        "matrices": I_TABLE,
        "generators": I_GENERATORS,
        "words": I_WORDS,
        "printed": {
            "delta": "exp(2 pi i/24)",
            "eta": "exp(2 pi i/8)",
            "omega": "exp(2 pi i/3)",
            "g+-": "(sqrt5 +- 1)/2",
            "s+-": "(sqrt3 +- 1)/sqrt2",
            "t": "sqrt5",
            "epsilon": "(8 + 3 sqrt3)^(1/2) exp[i arctan{(-6416 + 13216 sqrt3 + 2352 sqrt5 - 5075 sqrt15)/13189}]",
            "zeta": "(8 - 3 sqrt3)^(1/2) exp[i arctan{(6416 + 13216 sqrt3 - 2352 sqrt5 - 5075 sqrt15)/13189}]",
            "theta": "sqrt5 exp[i arctan{(4 - sqrt5)/sqrt3}]",
            "lambda": "exp(i arctan sqrt(5/3))",
            "mu": "(4 + sqrt3)^(1/2) exp[i arctan{(84 - 16 sqrt3 + 72 sqrt5 + 5 sqrt15)/131}]",
            "nu": "(4 - sqrt3)^(1/2) exp[i arctan{(84 + 16 sqrt3 + 72 sqrt5 - 5 sqrt15)/131}]",
            "xi": "sqrt5 exp(i arctan 2)",
            "rho": "sqrt7 exp[i arctan{(80 sqrt3 + 49 sqrt15)/57}]",
            "sigma": "sqrt5 exp{i arctan(28 + 16 sqrt3 + 12 sqrt5 + 7 sqrt15)}",
        },
    },
}

# --------------------------------------------------------------------------
# Evaluation in Q(zeta_M)
# --------------------------------------------------------------------------


class BookParser(cb.Parser):
    """The value grammar plus the book's symbols (and their conjugates, suffix s)."""

    def __init__(self, F, text, env):
        super().__init__(F, text)
        self.env = env

    def atom(self):
        kind, val = self.peek()
        if kind == "name" and val in self.env:
            self.take()
            return self.env[val]
        return super().atom()


def parse(F, text, env):
    p = BookParser(F, text, env)
    x = p.expr()
    assert p.pos == len(p.toks), text
    return x


def environment(F, symbols):
    env = {}
    for name, defn in symbols.items():
        x = cb.parse_value(F, defn)
        env[name] = x
        env[name + "s"] = F.conj(x)
    return env


def parse_matrix(F, text, env):
    pre = F.rational(1)
    if "|" in text:
        head, text = text.split("|")
        pre = cb.parse_value(F, head.strip())
    rows = [
        [a.strip() for a in r.split(",")] if "," in r else r.split()
        for r in text.split(";")
    ]
    d = len(rows)
    assert all(len(r) == d for r in rows), text
    return [[F.mul(pre, parse(F, a, env)) for a in r] for r in rows]


def galois(F, x, k):
    out = F.zero()
    for j, a in enumerate(x):
        if a:
            out = F.add(out, F.scale(F.zeta(j * k % F.conductor), a))
    return out


def in_subfield(F, x, d):
    M = F.conductor
    return all(
        galois(F, x, k) == x
        for k in range(1, M)
        if math.gcd(k, M) == 1 and (k - 1) % d == 0
    )


def minimal_conductor(F, values):
    M = F.conductor
    for d in sorted(k for k in range(1, M + 1) if M % k == 0):
        if all(in_subfield(F, x, d) for x in values):
            return d
    return M


# --------------------------------------------------------------------------
# Output in the value grammar
# --------------------------------------------------------------------------

SURDS = [
    ("1", ()),
    ("sqrt2", ("sqrt2",)),
    ("sqrt3", ("sqrt3",)),
    ("sqrt5", ("sqrt5",)),
    ("sqrt2*sqrt3", ("sqrt2", "sqrt3")),
    ("sqrt2*sqrt5", ("sqrt2", "sqrt5")),
    ("sqrt3*sqrt5", ("sqrt3", "sqrt5")),
    ("sqrt2*sqrt3*sqrt5", ("sqrt2", "sqrt3", "sqrt5")),
]
REQUIRES = {"sqrt2": 8, "sqrt3": 12, "sqrt5": 5}


def surd_basis(F):
    M = F.conductor
    out = []
    for name, parts in SURDS:
        if all(M % REQUIRES[p] == 0 for p in parts):
            x = F.rational(1)
            for p in parts:
                x = F.mul(x, cb.const(F, p))
            out.append((name, x))
    return out


def fmt_q(c):
    return cb.fmt_q(c)


def _combination(cs, names):
    """Rational combination as '(a + b*sqrt3)/den' with a common denominator."""
    den = 1
    for c in cs:
        den = den * c.denominator // math.gcd(den, c.denominator)
    terms = []
    for c, name in zip(cs, names):
        if c == 0:
            continue
        k = c * den
        assert k.denominator == 1
        k = k.numerator
        sign = "-" if k < 0 else "+"
        mag = abs(k)
        body = (
            fmt_q(Fr(mag)) if name == "1" else (name if mag == 1 else f"{mag}*{name}")
        )
        terms.append((sign, body))
    if not terms:
        return "0"
    s = ("-" if terms[0][0] == "-" else "") + terms[0][1]
    for sign, body in terms[1:]:
        s += f" {sign} {body}"
    if len(terms) > 1 and den != 1:
        s = f"({s})"
    if den != 1:
        s += f"/{den}"
    return s


def to_grammar(F, x):
    """Rational | r*zeta^k | r*surd*zeta^k | re + i*im over the surds | power basis."""
    M = F.conductor
    if all(a == 0 for a in x[1:]):
        return fmt_q(x[0])
    # r * zeta^k and r * surd * zeta^k with r rational
    for sname, sval in [("", F.rational(1))] + [
        (n, v) for n, v in surd_basis(F) if "*" not in n and n != "1"
    ]:
        for k in range(M):
            y = F.mul(x, F.conj(F.mul(sval, F.zeta(k))))  # x / (surd zeta^k) * |surd|^2
            if all(a == 0 for a in y[1:]):
                r = y[0] / (F.mul(sval, F.conj(sval))[0])
                if k == 0:
                    root = ""
                elif 4 * k == M:
                    root = "i"
                elif 2 * k == M:
                    root, r = "", -r
                elif 4 * k == 3 * M:
                    root, r = "i", -r
                else:
                    root = f"zeta^{k}" if k != 1 else "zeta"
                parts = [p for p in (sname, root) if p]
                num, den = abs(r.numerator), r.denominator
                head = "" if num == 1 and parts else str(num)
                body = "*".join(([head] if head else []) + parts)
                if den != 1:
                    body += f"/{den}"
                return ("-" if r < 0 else "") + body
    # re + i*im over the surd basis
    basis = surd_basis(F)
    names = [n for n, _ in basis]
    vecs = [v for _, v in basis]
    two = Fr(2)
    re_ = F.scale(F.add(x, F.conj(x)), Fr(1, 2))
    im = F.scale(F.mul(F.sub(x, F.conj(x)), F.conj(F.i())), Fr(1, 2))
    del two
    cre = cb.coords(F, re_, vecs)
    cim = cb.coords(F, im, vecs)
    if cre is not None and cim is not None:
        parts = []
        if any(cre):
            parts.append(_combination(cre, names))
        if any(cim):
            s = _combination(cim, names)
            if s.startswith("-") and not s.startswith("("):
                parts.append("- i*" + s[1:] if parts else "-i*" + s[1:])
            elif s.startswith("-("):
                parts.append(("- i*" if parts else "-i*") + s[1:])
            else:
                parts.append(("+ i*" if parts else "i*") + s)
        return " ".join(parts)
    return cb._terms(
        (a, f"zeta^{j}" if j > 1 else "zeta" if j == 1 else "1")
        for j, a in enumerate(x)
    )


# --------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------


class Transcription:
    def __init__(self, name):
        spec = GROUPS[name]
        self.name = name
        self.kind = spec["kind"]
        self.spec = spec
        self.M = spec["conductor"]
        self.F = cb.field(self.M)
        self.env = environment(self.F, spec["symbols"])
        self.irreps = spec["irreps"]
        self.book, self.ops = ab.book_group(name)
        self.labels = self.book.labels
        # the printed matrices
        self.errata = ERRATA.get(name, [])
        self.printed = {}  # label -> {irrep: matrix}, errata applied
        for label, mats in spec["matrices"].items():
            assert label in self.book.index, label
            assert len(mats) == len(self.irreps), label
            self.printed[label] = {
                ir: self.corrected(label, ir, text)
                for ir, text in zip(self.irreps, mats)
            }
        # the canonical group and its cover
        self.P = cb.Poly(self.kind, cb.BASE_CONDUCTOR[self.kind])
        self.q_of = {}  # book label -> exact quaternion of the book's untilded operation
        numeric = {}
        for q in self.P.Q:
            numeric[self._key(self._float(q))] = q
        for r, p, label in self.book.elements:
            assert p == +1
            self.q_of[label] = numeric[self._key(self.ops[r].quaternion)]

    def corrected(self, label, ir, text):
        env = self.env
        transpose = False
        for lab, irrep, fix in self.errata:
            if (lab, irrep) != (label, ir):
                continue
            if fix.startswith("negate "):
                sym = fix.split()[1]
                env = dict(env)
                env[sym] = self.F.neg(env[sym])
                env[sym + "s"] = self.F.neg(env[sym + "s"])
            elif fix == "transpose":
                transpose = True
            else:
                raise ValueError(fix)
        A = parse_matrix(self.F, text, env)
        if transpose:
            A = [list(col) for col in zip(*A)]
        return A

    def _float(self, q):
        alpha = {"T": 1.0, "O": math.sqrt(2), "I": (1 + math.sqrt(5)) / 2}[self.kind]
        return tuple(float(a) + float(b) * alpha for a, b in q)

    @staticmethod
    def _key(v):
        return tuple(round(c, 6) + 0.0 for c in v)

    def qmul(self, a, b):
        return self.P.cov.mul(a, b)

    @staticmethod
    def qneg(q):
        return tuple((-a, -b) for a, b in q)

    def lift(self, x):
        return cb.lift(self.P.F, self.F, x)

    # ---- all elements: printed or from the words ----

    def all_matrices(self, report):
        """label -> (quaternion, {irrep: matrix}); for I the products along the words."""
        words = self.spec["words"]
        if words is None:
            assert set(self.printed) == set(self.labels)
            return {lab: (self.q_of[lab], self.printed[lab]) for lab in self.labels}
        assert set(words) == set(self.labels)
        out = {}
        tilded = []
        for lab in self.labels:
            gens = words[lab].split()
            q = None
            mats = None
            for g in gens:
                qg = self.q_of[g]
                q = qg if q is None else self.qmul(q, qg)
                mg = self.printed[g]
                mats = (
                    dict(mg)
                    if mats is None
                    else {ir: cb.mmul(self.F, mats[ir], mg[ir]) for ir in self.irreps}
                )
            if q == self.q_of[lab]:
                pass
            elif q == self.qneg(self.q_of[lab]):
                tilded.append(lab)
            else:
                report.fail(f"word for {lab} ({words[lab]}) is not {lab} or its tilde")
            out[lab] = (q, mats)
        report.note(
            f"words of T 74.7a: {len(self.labels) - len(tilded)} give the untilded operation, "
            f"{len(tilded)} the tilded one"
            + (f" ({', '.join(tilded)})" if tilded else "")
        )
        return out

    # ---- checks ----

    def check(self, report):
        F = self.F
        elems = self.all_matrices(report)
        index = {q: lab for lab, (q, _) in elems.items()}
        neg_index = {self.qneg(q): lab for lab, (q, _) in elems.items()}
        # unitarity of every printed matrix
        for lab, mats in self.printed.items():
            for ir, A in mats.items():
                if not cb.is_unitary(F, A):
                    report.fail(f"{lab} {ir}: not unitary")
        # multiplication with the sign of the quaternion product
        n_pairs = 0
        right = self.labels if self.spec["words"] is None else self.spec["generators"]
        for a, (qa, ma) in elems.items():
            for b in right:
                qb, mb = elems[b]
                q = self.qmul(qa, qb)
                if q in index:
                    c, sign = index[q], 1
                else:
                    c, sign = neg_index[q], -1
                _, mc = elems[c]
                for ir in self.irreps:
                    spinor = "/" in ir
                    prod = cb.mmul(F, ma[ir], mb[ir])
                    expect = mc[ir]
                    if spinor and sign < 0:
                        expect = cb.mscale(F, expect, F.rational(-1))
                    if not cb.meq(prod, expect):
                        report.fail(
                            f"{ir}: D({a}) D({b}) != {'-' if spinor and sign < 0 else ''}D({c})"
                        )
                n_pairs += 1
        report.note(
            f"{n_pairs} products checked"
            + (
                ""
                if self.spec["words"] is None
                else " (every operation times every generator)"
            )
        )
        # the tilde marks of the labels file agree with the quaternion signs
        rows = self.label_rows()
        if rows is not None:
            for a, (qa, _) in elems.items():
                for b, (qb, _) in elems.items():
                    q = self.qmul(qa, qb)
                    printed = rows[self.labels.index(a)][self.labels.index(b)]
                    expect = index[q] if q in index else "~" + neg_index[q]
                    if self.spec["words"] is None and printed != expect:
                        report.fail(
                            f"tilde mark {a} * {b}: labels file {printed}, quaternions {expect}"
                        )
        # traces identify one canonical irrep per book symbol
        self.identified = {}
        for ir in self.irreps:
            d = len(self.printed[self.labels[1]][ir])
            hits = []
            for cir in self.P.ctx.irreps:
                if cir.dimension != d:
                    continue
                chi, _ = self.P.character(cir.symbol)
                if all(
                    cb.mtrace(F, mats[ir]) == self.lift(chi(q))
                    for lab, (q, mats) in elems.items()
                ):
                    hits.append(cir.symbol)
            if len(hits) != 1:
                report.fail(f"{ir}: trace matches canonical irreps {hits}")
            else:
                self.identified[ir] = hits[0]
        report.note(
            "book symbol -> canonical: "
            + ", ".join(f"{k} -> {v}" for k, v in self.identified.items())
        )
        # minimal conductor per irrep
        conds = {}
        for ir in self.irreps:
            vals = [
                a for lab, mats in self.printed.items() for row in mats[ir] for a in row
            ]
            conds[ir] = minimal_conductor(F, vals)
        self.conductors = conds
        report.note(
            "minimal conductor per irrep: "
            + ", ".join(f"{k}: {v}" for k, v in conds.items())
        )
        # the identified surds reproduce the printed formulas
        if self.name == "I":
            for name, (r, phase) in I_PRINTED.items():
                z = r * cmath.exp(1j * phase)
                w = F.complex(self.env[name])
                if abs(z - w) > 1e-9:
                    report.fail(f"constant {name}: identified {w} != printed {z}")
            report.note(
                f"{len(I_PRINTED)} printed constants agree with their identified surds"
            )
        # every generator matrix round-trips through the grammar
        for lab, mats in self.printed.items():
            for ir, A in mats.items():
                for row in A:
                    for a in row:
                        s = to_grammar(F, a)
                        if cb.parse_value(F, s) != a:
                            report.fail(
                                f"{lab} {ir}: grammar round trip failed for {s}"
                            )
        return elems

    def label_rows(self):
        path = LABEL_DIR / f"{self.name}.toml"
        if not path.exists():
            return None
        data = tomllib.loads(path.read_text())
        assert data["book"]["elements"] == self.labels
        return data["multiplication"]["rows"]


class Report:
    def __init__(self):
        self.failures = []
        self.notes = []

    def fail(self, msg):
        self.failures.append(msg)
        print("   FAIL", msg)

    def note(self, msg):
        self.notes.append(msg)
        print("  ", msg)


def check_all():
    ok = True
    results = {}
    for name in GROUPS:
        print(f"{name}:")
        rep = Report()
        tr = Transcription(name)
        tr.check(rep)
        results[name] = (tr, rep)
        if rep.failures:
            ok = False
            print(f"   {len(rep.failures)} failure(s)")
        else:
            print("   ok")
    return ok, results


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


def toml_str(s):
    return ab.toml_str(s)


def matrix_lines(F, A, indent):
    out = []
    for row in A:
        out.append(
            indent + "[" + ", ".join(toml_str(to_grammar(F, a)) for a in row) + "],"
        )
    return out


def generate(results, out_dir=OUT_DIR):
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, (tr, rep) in results.items():
        spec = tr.spec
        F = tr.F
        page = spec["table"]
        lines = [
            f"# Altmann & Herzig, Point-Group Theory Tables (2nd ed.): irrep matrices of {name}",
            "# Format: tests/reference/altmann-herzig/matrices/README.md",
            "# Written by tools/ah_matrices.py from the transcription it holds; every matrix is checked",
            "# exactly (unitarity, products against the book's factor system, traces against the",
            "# canonical characters) before this file is written.",
            "",
            f"group = {toml_str(name)}",
            f"table = {toml_str(page)}",
            'status = "done"',
            f"source = {toml_str(SOURCE)}",
            'transcribed_by = "agent (Claude) from page images of the PDF; verify against the print"',
            f"conductor = {tr.M}",
            'notes = """',
        ]
        if spec["words"] is None:
            lines.append(
                "The book prints one matrix per operation; all of them are transcribed, in the book's"
            )
            lines.append(
                "element order, for every irrep of dimension >= 2 (one-dimensional irreps are their characters)."
            )
        else:
            lines.append(
                "The book prints matrices for three generators only (T 74.7b) and a word in those"
            )
            lines.append(
                "generators for every operation (T 74.7a, [words] below, left to right as printed);"
            )
            lines.append(
                "the matrix of an operation is the product of the generator matrices along its word."
            )
        for n in rep.notes:
            lines.append(n + ".")
        if tr.errata:
            lines.append(
                "Errata: with delta = exp(2 pi i/24) the printed C3,1+ matrices of F3/2 and I5/2 cube to +1"
            )
            lines.append(
                "on their delta blocks (a spinor C3 must cube to -1; T 70.7 prints the same blocks with the"
            )
            lines.append(
                "opposite sign), and the printed C5,1+ matrices of F3/2 and I5/2 have rows and columns"
            )
            lines.append(
                "interchanged. The [[known_discrepancy]] entries record the printed forms; the matrices below carry"
            )
            lines.append(
                "the corrections, which are the only ones that make the products along every word agree"
            )
            lines.append("with the factor system and the traces with the characters.")
        lines.append(
            "Entries are exact in Q(zeta_M), M = conductor, in the value grammar of"
        )
        lines.append(
            "docs/spec/basis-conventions.md; the book's symbols are listed in [book.constants]."
        )
        lines.append('"""')
        lines.append("")
        lines.append("[book]")
        lines.append(f"irreps = [{', '.join(toml_str(s) for s in tr.irreps)}]")
        lines.append(
            f"elements = [{', '.join(toml_str(s) for s in spec['matrices'])}]"
            "  # book labels of the transcribed matrices"
        )
        lines.append(
            f"generators = [{', '.join(toml_str(s) for s in spec['generators'])}]"
            "  # a generating set among them (book labels; canonical labels via labels/<Group>.toml)"
        )
        lines.append("")
        lines.append("[book.constants]  # as printed in the table footnote")
        for k, v in spec["printed"].items():
            lines.append(f"{toml_str(k)} = {toml_str(v)}")
        lines.append("")
        lines.append(
            "[book.exact]  # the same constants in the value grammar (conjugates carry a trailing s)"
        )
        for k, v in spec["symbols"].items():
            lines.append(f"{k} = {toml_str(v)}")
        lines.append("")
        if tr.errata:
            for lab, ir, fix in tr.errata:
                lines.append(
                    "[[known_discrepancy]]  # a printed matrix corrected before use (docs/spec/verification.md)"
                )
                lines.append(f"element = {toml_str(lab)}")
                lines.append(f"irrep = {toml_str(ir)}")
                lines.append(f"fix = {toml_str(fix)}")
                lines.append(
                    f"printed = {toml_str(spec['matrices'][lab][tr.irreps.index(ir)])}"
                    "  # as printed, in the symbols of [book.exact]"
                )
                lines.append(f"note = {toml_str(ERRATA_NOTE[fix])}")
                lines.append("")
        lines.append(
            "[canonical]  # book irrep symbol -> canonical Mulliken symbol, by the traces"
        )
        for k, v in tr.identified.items():
            lines.append(f"{toml_str(k)} = {toml_str(v)}")
        lines.append("")
        lines.append(
            "[minimal_conductor]  # smallest Q(zeta_d) holding every printed entry of the irrep"
        )
        for k, v in tr.conductors.items():
            lines.append(f"{toml_str(k)} = {v}")
        if spec["words"] is not None:
            lines.append("")
            lines.append(
                "[words]  # T 74.7a: every operation as a product of generators, left to right"
            )
            for lab in tr.labels:
                lines.append(
                    f"{toml_str(lab)} = [{', '.join(toml_str(g) for g in spec['words'][lab].split())}]"
                )
        for ir in tr.irreps:
            lines.append("")
            lines.append(f"[irreps.{toml_str(ir)}]")
            lines.append(
                "matrices = [  # one per element of [book].elements, rows then columns"
            )
            for lab in spec["matrices"]:
                A = tr.printed[lab][ir]
                lines.append(f"  [  # {lab}")
                lines.extend(matrix_lines(F, A, "    "))
                lines.append("  ],")
            lines.append("]")
        (out_dir / f"{name}.toml").write_text("\n".join(lines) + "\n")
        print(f"wrote {out_dir / f'{name}.toml'}")


def main(argv):
    if len(argv) < 2 or argv[1] not in ("check", "generate"):
        print(__doc__)
        return 2
    ok, results = check_all()
    if not ok:
        print("transcription has failures; nothing written")
        return 1
    if argv[1] == "generate":
        generate(results)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
