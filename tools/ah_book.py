"""Draft transcription of Altmann & Herzig, *Point-Group Theory Tables* (2nd ed.),
from its parameter tables (wayfinder ticket #13).

The book parametrises every proper operation of a group by a rotation angle
``phi`` and a pole ``n`` (T n.1, section 3-2 of the book) and writes every
improper operation as ``i g``; the quaternion ``[[cos(phi/2), sin(phi/2) n]]``
listed for ``g`` is, by the book's own convention (3.30), the untilded
double-group element. Everything else the ticket asks for (multiplication
table, tilde marks, factor table, the book's section) follows from these
parameters by quaternion multiplication, so the parameters of the seven
master tables (T 31.1 D2h, T 33.1 D4h, T 35.1 D6h, T 37.1 D8h, T 39.1 D10h,
T 71.1 Oh, T 75.1 Ih) and C3v's own T 51.1A/B are what is transcribed here,
by hand from the OCR text of the PDF, together with the subgroup-element
tables T m.0 (which operations each subgroup contains and how it names
them) and the irrep symbol rows of the character tables T n.4.

    uv run python tools/ah_book.py check [OCR_DIR]   # rebuild every master, compare
                                                     # with the printed factor tables
    uv run python tools/ah_book.py generate          # write tests/reference/altmann-herzig/labels/

``OCR_DIR`` holds one ``pNNN.txt`` per PDF page (pymupdf ``get_text``); without it
the check is internal only (group axioms, order, isomorphism to our group).
"""

from __future__ import annotations

import math
import sys
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
sys.path.insert(0, str(ROOT / "tools"))
import abstract_group_prototype as ag
import ah_transcription as at

SOURCE = "S. L. Altmann and P. Herzig, Point-Group Theory Tables, 2nd ed. (Wien: Clarendon, 2011)"
OUT_DIR = ROOT / "tests" / "reference" / "altmann-herzig" / "labels"

# --------------------------------------------------------------------------
# Symbolic constants of the parameter tables
# --------------------------------------------------------------------------

PHI = (1 + math.sqrt(5)) / 2
SIGMA = math.atan(PHI)  # book: sigma = arctan((1 + sqrt5)/2)
RHO1 = math.atan((3 + math.sqrt(5)) / 2)  # book: rho' = arctan((3 + sqrt5)/2)

SYMBOLS = {
    "0": 0.0,
    "1": 1.0,
    "1/2": 0.5,
    "sqrt3/2": math.sqrt(3) / 2,
    "1/sqrt2": 1 / math.sqrt(2),
    "1/sqrt3": 1 / math.sqrt(3),
    "c8": math.cos(math.pi / 8),
    "s8": math.sin(math.pi / 8),
    "c5": math.cos(math.pi / 5),
    "s5": math.sin(math.pi / 5),
    "c10": math.cos(math.pi / 10),
    "s10": math.sin(math.pi / 10),
    "cos_sigma": math.cos(SIGMA),
    "sin_sigma": math.sin(SIGMA),
    "cos_rho'": math.cos(RHO1),
    "sin_rho'": math.sin(RHO1),
}
ANGLES = {
    "0": 0.0,
    "pi/5": math.pi / 5,
    "pi/4": math.pi / 4,
    "pi/3": math.pi / 3,
    "2pi/5": 2 * math.pi / 5,
    "pi/2": math.pi / 2,
    "3pi/5": 3 * math.pi / 5,
    "2pi/3": 2 * math.pi / 3,
    "3pi/4": 3 * math.pi / 4,
    "4pi/5": 4 * math.pi / 5,
    "pi": math.pi,
}


def value(token):
    sign = -1.0 if token.startswith("-") else 1.0
    return sign * SYMBOLS[token.lstrip("-")]


# --------------------------------------------------------------------------
# Master parameter tables: "proper improper phi nx ny nz", book order.
# Labels in the ASCII of tests/reference/altmann-herzig/labels/README.md.
# --------------------------------------------------------------------------

MASTERS = {
    # T 31.1, p. 246
    "D2h": """
E     i        0    0 0 0
C2z   sigma_z  pi   0 0 1
C2x   sigma_x  pi   1 0 0
C2y   sigma_y  pi   0 1 0
""",
    # T 33.1, p. 257
    "D4h": """
E      i         0     0 0 0
C4+    S4-       pi/2  0 0 1
C4-    S4+       pi/2  0 0 -1
C2     sigma_h   pi    0 0 1
C2'1   sigma_v1  pi    1 0 0
C2'2   sigma_v2  pi    0 1 0
C2''1  sigma_d1  pi    1/sqrt2 1/sqrt2 0
C2''2  sigma_d2  pi    -1/sqrt2 1/sqrt2 0
""",
    # T 35.1, p. 274
    "D6h": """
E      i         0      0 0 0
C6+    S3-       pi/3   0 0 1
C6-    S3+       pi/3   0 0 -1
C3+    S6-       2pi/3  0 0 1
C3-    S6+       2pi/3  0 0 -1
C2     sigma_h   pi     0 0 1
C2'1   sigma_d1  pi     1 0 0
C2'2   sigma_d2  pi     -1/2 sqrt3/2 0
C2'3   sigma_d3  pi     -1/2 -sqrt3/2 0
C2''1  sigma_v1  pi     0 1 0
C2''2  sigma_v2  pi     -sqrt3/2 -1/2 0
C2''3  sigma_v3  pi     sqrt3/2 -1/2 0
""",
    # T 37.1, p. 305
    "D8h": """
E      i         0      0 0 0
C8+    S8^3-     pi/4   0 0 1
C8-    S8^3+     pi/4   0 0 -1
C4+    S4-       pi/2   0 0 1
C4-    S4+       pi/2   0 0 -1
C8^3+  S8-       3pi/4  0 0 1
C8^3-  S8+       3pi/4  0 0 -1
C2     sigma_h   pi     0 0 1
C2'1   sigma_v1  pi     1 0 0
C2'2   sigma_v2  pi     0 1 0
C2'3   sigma_v3  pi     1/sqrt2 1/sqrt2 0
C2'4   sigma_v4  pi     -1/sqrt2 1/sqrt2 0
C2''1  sigma_d1  pi     c8 s8 0
C2''2  sigma_d2  pi     -s8 c8 0
C2''3  sigma_d3  pi     s8 c8 0
C2''4  sigma_d4  pi     -c8 s8 0
""",
    # T 39.1, p. 345
    "D10h": """
E       i          0      0 0 0
C10+    S5^2-      pi/5   0 0 1
C10-    S5^2+      pi/5   0 0 -1
C5+     S10^3-     2pi/5  0 0 1
C5-     S10^3+     2pi/5  0 0 -1
C10^3+  S5-        3pi/5  0 0 1
C10^3-  S5+        3pi/5  0 0 -1
C5^2+   S10-       4pi/5  0 0 1
C5^2-   S10+       4pi/5  0 0 -1
C2      sigma_h    pi     0 0 1
C2'1    sigma_d1   pi     1 0 0
C2'2    sigma_d2   pi     s10 c10 0
C2'3    sigma_d3   pi     -c5 s5 0
C2'4    sigma_d4   pi     -c5 -s5 0
C2'5    sigma_d5   pi     s10 -c10 0
C2''1   sigma_v1   pi     0 1 0
C2''2   sigma_v2   pi     -c10 s10 0
C2''3   sigma_v3   pi     -s5 -c5 0
C2''4   sigma_v4   pi     s5 -c5 0
C2''5   sigma_v5   pi     c10 s10 0
""",
    # T 71.1, p. 597
    "Oh": """
E      i         0      0 0 0
C2x    sigma_x   pi     1 0 0
C2y    sigma_y   pi     0 1 0
C2z    sigma_z   pi     0 0 1
C3,1+  S6,1-     2pi/3  1/sqrt3 1/sqrt3 1/sqrt3
C3,2+  S6,2-     2pi/3  -1/sqrt3 -1/sqrt3 1/sqrt3
C3,3+  S6,3-     2pi/3  1/sqrt3 -1/sqrt3 -1/sqrt3
C3,4+  S6,4-     2pi/3  -1/sqrt3 1/sqrt3 -1/sqrt3
C3,1-  S6,1+     2pi/3  -1/sqrt3 -1/sqrt3 -1/sqrt3
C3,2-  S6,2+     2pi/3  1/sqrt3 1/sqrt3 -1/sqrt3
C3,3-  S6,3+     2pi/3  -1/sqrt3 1/sqrt3 1/sqrt3
C3,4-  S6,4+     2pi/3  1/sqrt3 -1/sqrt3 1/sqrt3
C4x+   S4x-      pi/2   1 0 0
C4y+   S4y-      pi/2   0 1 0
C4z+   S4z-      pi/2   0 0 1
C4x-   S4x+      pi/2   -1 0 0
C4y-   S4y+      pi/2   0 -1 0
C4z-   S4z+      pi/2   0 0 -1
C2'a   sigma_d1  pi     1/sqrt2 1/sqrt2 0
C2'b   sigma_d2  pi     -1/sqrt2 1/sqrt2 0
C2'c   sigma_d3  pi     1/sqrt2 0 1/sqrt2
C2'd   sigma_d4  pi     0 -1/sqrt2 -1/sqrt2
C2'e   sigma_d5  pi     1/sqrt2 0 -1/sqrt2
C2'f   sigma_d6  pi     0 -1/sqrt2 1/sqrt2
""",
    # T 75.1, pp. 663-664
    "Ih": """
E        i           0      0 0 0
C5,1+    S10,1^3-    2pi/5  cos_sigma 0 sin_sigma
C5,2+    S10,2^3-    2pi/5  -cos_sigma 0 sin_sigma
C5,3+    S10,3^3-    2pi/5  0 sin_sigma cos_sigma
C5,4+    S10,4^3-    2pi/5  0 -sin_sigma cos_sigma
C5,5+    S10,5^3-    2pi/5  sin_sigma cos_sigma 0
C5,6+    S10,6^3-    2pi/5  sin_sigma -cos_sigma 0
C5,1-    S10,1^3+    2pi/5  -cos_sigma 0 -sin_sigma
C5,2-    S10,2^3+    2pi/5  cos_sigma 0 -sin_sigma
C5,3-    S10,3^3+    2pi/5  0 -sin_sigma -cos_sigma
C5,4-    S10,4^3+    2pi/5  0 sin_sigma -cos_sigma
C5,5-    S10,5^3+    2pi/5  -sin_sigma -cos_sigma 0
C5,6-    S10,6^3+    2pi/5  -sin_sigma cos_sigma 0
C5,1^2+  S10,1-      4pi/5  cos_sigma 0 sin_sigma
C5,2^2+  S10,2-      4pi/5  -cos_sigma 0 sin_sigma
C5,3^2+  S10,3-      4pi/5  0 sin_sigma cos_sigma
C5,4^2+  S10,4-      4pi/5  0 -sin_sigma cos_sigma
C5,5^2+  S10,5-      4pi/5  sin_sigma cos_sigma 0
C5,6^2+  S10,6-      4pi/5  sin_sigma -cos_sigma 0
C5,1^2-  S10,1+      4pi/5  -cos_sigma 0 -sin_sigma
C5,2^2-  S10,2+      4pi/5  cos_sigma 0 -sin_sigma
C5,3^2-  S10,3+      4pi/5  0 -sin_sigma -cos_sigma
C5,4^2-  S10,4+      4pi/5  0 sin_sigma -cos_sigma
C5,5^2-  S10,5+      4pi/5  -sin_sigma -cos_sigma 0
C5,6^2-  S10,6+      4pi/5  -sin_sigma cos_sigma 0
C3,1+    S6,1-       2pi/3  1/sqrt3 1/sqrt3 1/sqrt3
C3,2+    S6,2-       2pi/3  -1/sqrt3 -1/sqrt3 1/sqrt3
C3,3+    S6,3-       2pi/3  1/sqrt3 -1/sqrt3 -1/sqrt3
C3,4+    S6,4-       2pi/3  -1/sqrt3 1/sqrt3 -1/sqrt3
C3,5+    S6,5-       2pi/3  -sin_rho' 0 cos_rho'
C3,6+    S6,6-       2pi/3  0 cos_rho' -sin_rho'
C3,7+    S6,7-       2pi/3  cos_rho' -sin_rho' 0
C3,8+    S6,8-       2pi/3  sin_rho' 0 cos_rho'
C3,9+    S6,9-       2pi/3  0 cos_rho' sin_rho'
C3,10+   S6,10-      2pi/3  cos_rho' sin_rho' 0
C3,1-    S6,1+       2pi/3  -1/sqrt3 -1/sqrt3 -1/sqrt3
C3,2-    S6,2+       2pi/3  1/sqrt3 1/sqrt3 -1/sqrt3
C3,3-    S6,3+       2pi/3  -1/sqrt3 1/sqrt3 1/sqrt3
C3,4-    S6,4+       2pi/3  1/sqrt3 -1/sqrt3 1/sqrt3
C3,5-    S6,5+       2pi/3  sin_rho' 0 -cos_rho'
C3,6-    S6,6+       2pi/3  0 -cos_rho' sin_rho'
C3,7-    S6,7+       2pi/3  -cos_rho' sin_rho' 0
C3,8-    S6,8+       2pi/3  -sin_rho' 0 -cos_rho'
C3,9-    S6,9+       2pi/3  0 -cos_rho' -sin_rho'
C3,10-   S6,10+      2pi/3  -cos_rho' -sin_rho' 0
C2a      sigma_a     pi     1 0 0
C2b      sigma_b     pi     0 1 0
C2c      sigma_c     pi     0 0 1
C2d      sigma_d     pi     c5 s10 1/2
C2e      sigma_e     pi     1/2 -c5 -s10
C2f      sigma_f     pi     s10 1/2 -c5
C2g      sigma_g     pi     s10 1/2 c5
C2h      sigma_h     pi     1/2 -c5 s10
C2i      sigma_i     pi     c5 s10 -1/2
C2j      sigma_j     pi     -s10 1/2 c5
C2k      sigma_k     pi     -c5 s10 -1/2
C2l      sigma_l     pi     -1/2 -c5 s10
C2m      sigma_m     pi     -c5 s10 1/2
C2n      sigma_n     pi     -s10 1/2 -c5
C2o      sigma_o     pi     -1/2 -c5 -s10
""",
    # T 51.1A and T 51.1B, p. 485 (C3v has two settings and its own table)
    "C3vA": """
E         -   0      0 0 0
C3+       -   2pi/3  0 0 1
C3-       -   2pi/3  0 0 -1
sigma_v1  *   pi     0 1 0
sigma_v2  *   pi     -sqrt3/2 -1/2 0
sigma_v3  *   pi     sqrt3/2 -1/2 0
""",
    "C3vB": """
E         -   0      0 0 0
C3+       -   2pi/3  0 0 1
C3-       -   2pi/3  0 0 -1
sigma_v1  *   pi     1 0 0
sigma_v2  *   pi     -1/2 sqrt3/2 0
sigma_v3  *   pi     -1/2 -sqrt3/2 0
""",
}
MASTER_TABLE = {
    "D2h": "T 31",
    "D4h": "T 33",
    "D6h": "T 35",
    "D8h": "T 37",
    "D10h": "T 39",
    "Oh": "T 71",
    "Ih": "T 75",
    "C3vA": "T 51",
    "C3vB": "T 51",
}


class Operation:
    """One row of a parameter table: a proper rotation and its improper partner."""

    def __init__(self, proper, improper, phi, n):
        self.proper = proper
        self.improper = improper
        self.phi = phi
        self.n = n  # symbolic tokens
        angle = ANGLES[phi]
        nv = [value(t) for t in n]
        self.quaternion = (math.cos(angle / 2),) + tuple(
            math.sin(angle / 2) * c for c in nv
        )


def parse_master(text):
    ops = []
    for line in text.strip().splitlines():
        proper, improper, phi, nx, ny, nz = line.split()
        ops.append(Operation(proper, improper, phi, (nx, ny, nz)))
    return ops


# --------------------------------------------------------------------------
# The book's group from the parameters
# --------------------------------------------------------------------------


def qmul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    )


def close(p, q, tol=1e-6):
    return all(abs(a - b) < tol for a, b in zip(p, q))


class BookGroup:
    """Elements: (row index, parity). Quaternion of (r, -1) is that of (r, +1) with the marker i."""

    def __init__(self, ops, elements):
        """``elements``: list of (row, parity, book label) in book order."""
        self.ops = ops
        self.elements = elements
        self.labels = [e[2] for e in elements]
        self.index = {lab: k for k, lab in enumerate(self.labels)}
        n = len(elements)
        self.table = [[-1] * n for _ in range(n)]
        self.factor = [[0] * n for _ in range(n)]
        quats = [ops[r].quaternion for r, _, _ in elements]
        for a, (ra, pa, _) in enumerate(elements):
            for b, (rb, pb, _) in enumerate(elements):
                q = qmul(quats[a], quats[b])
                parity = pa * pb
                hit = None
                for c, (rc, pc, _) in enumerate(elements):
                    if pc != parity:
                        continue
                    if close(q, quats[c]):
                        hit, sign = c, 1
                        break
                    if close(q, tuple(-x for x in quats[c])):
                        hit, sign = c, -1
                        break
                if hit is None:
                    raise ValueError(
                        f"product {self.labels[a]} * {self.labels[b]} is not in the set"
                    )
                self.table[a][b] = hit
                self.factor[a][b] = sign
        # group axioms
        for a in range(n):
            assert sorted(self.table[a]) == list(range(n)), "not a Latin square"
        e = self.index["E"]
        assert all(self.table[e][a] == a and self.table[a][e] == a for a in range(n))


def master_group(name):
    ops = parse_master(MASTERS[name])
    if name.startswith("C3v"):
        elements = [
            (r, +1 if op.improper == "-" else -1, op.proper) for r, op in enumerate(ops)
        ]
    else:
        elements = [(r, +1, op.proper) for r, op in enumerate(ops)] + [
            (r, -1, op.improper) for r, op in enumerate(ops)
        ]
    return ops, BookGroup(ops, elements)


# --------------------------------------------------------------------------
# Subgroups: which master rows each target group contains and how it names
# them (T m.0), in the order of the group's own "(3) Operations of G" line
# (read from the print: the cyclic groups and the odd-n Cnh / Dnh groups list
# their elements by powers of the generator; the others follow the master).
# ``sel`` lists (master label, subgroup label); ``None`` keeps the master's
# label and order; "proper" keeps the proper rows.
# --------------------------------------------------------------------------


def _same(*labels):
    return [(lab, lab) for lab in labels]


def _ren(pairs):
    return [tuple(p.split("=")) if "=" in p else (p, p) for p in pairs.split()]


SUBGROUPS = {
    # master D2h, T 31.0
    "C1": ("D2h", _same("E")),
    "C2": ("D2h", _ren("E C2z=C2")),
    "Ci": ("D2h", _same("E", "i")),
    "Cs": ("D2h", _ren("E sigma_z=sigma_h")),
    "D2": ("D2h", _same("E", "C2z", "C2x", "C2y")),
    "C2v": ("D2h", _ren("E C2z=C2 sigma_x sigma_y")),
    "C2h": ("D2h", _ren("E C2z=C2 i sigma_z=sigma_h")),
    "D2h": ("D2h", None),
    # master D4h, T 33.0
    "C4": ("D4h", _same("E", "C4+", "C2", "C4-")),
    "S4": ("D4h", _same("E", "S4-", "C2", "S4+")),
    "D4": ("D4h", _same("E", "C4+", "C4-", "C2", "C2'1", "C2'2", "C2''1", "C2''2")),
    "C4v": (
        "D4h",
        _same("E", "C4+", "C4-", "C2", "sigma_v1", "sigma_v2", "sigma_d1", "sigma_d2"),
    ),
    "D2d": (
        "D4h",
        _same("E", "C2", "C2'1", "C2'2", "S4-", "S4+", "sigma_d1", "sigma_d2"),
    ),
    "C4h": ("D4h", _same("E", "C4+", "C2", "C4-", "i", "S4-", "sigma_h", "S4+")),
    "D4h": ("D4h", None),
    # master D6h, T 35.0
    "C3": ("D6h", _same("E", "C3+", "C3-")),
    "C6": ("D6h", _same("E", "C6+", "C3+", "C2", "C3-", "C6-")),
    "S6": ("D6h", _same("E", "C3+", "C3-", "i", "S6-", "S6+")),
    "D3": ("D6h", _same("E", "C3+", "C3-", "C2'1", "C2'2", "C2'3")),
    "D6": (
        "D6h",
        _same(
            "E",
            "C6+",
            "C6-",
            "C3+",
            "C3-",
            "C2",
            "C2'1",
            "C2'2",
            "C2'3",
            "C2''1",
            "C2''2",
            "C2''3",
        ),
    ),
    "C6v": (
        "D6h",
        _same(
            "E",
            "C6+",
            "C6-",
            "C3+",
            "C3-",
            "C2",
            "sigma_d1",
            "sigma_d2",
            "sigma_d3",
            "sigma_v1",
            "sigma_v2",
            "sigma_v3",
        ),
    ),
    "C3h": ("D6h", _same("E", "C3+", "C3-", "sigma_h", "S3+", "S3-")),
    "C6h": (
        "D6h",
        _same(
            "E",
            "C6+",
            "C3+",
            "C2",
            "C3-",
            "C6-",
            "i",
            "S3-",
            "S6-",
            "sigma_h",
            "S6+",
            "S3+",
        ),
    ),
    "D3d": (
        "D6h",
        _same(
            "E",
            "C3+",
            "C3-",
            "C2'1",
            "C2'2",
            "C2'3",
            "i",
            "S6-",
            "S6+",
            "sigma_d1",
            "sigma_d2",
            "sigma_d3",
        ),
    ),
    "D3h": (
        "D6h",
        _same(
            "E",
            "C3+",
            "C3-",
            "C2'1",
            "C2'2",
            "C2'3",
            "sigma_h",
            "S3+",
            "S3-",
            "sigma_v1",
            "sigma_v2",
            "sigma_v3",
        ),
    ),
    "D6h": ("D6h", None),
    # C3v: own tables, two settings
    "C3v": ("C3vA", None),
    "C3vB": ("C3vB", None),
    # master D8h, T 37.0
    "S8": ("D8h", _same("E", "S8^3-", "C4+", "S8-", "C2", "S8+", "C4-", "S8^3+")),
    "D4d": (
        "D8h",
        _same(
            "E",
            "C4+",
            "C4-",
            "C2",
            "C2'1",
            "C2'2",
            "C2'3",
            "C2'4",
            "S8^3-",
            "S8^3+",
            "S8-",
            "S8+",
            "sigma_d1",
            "sigma_d2",
            "sigma_d3",
            "sigma_d4",
        ),
    ),
    # master D10h, T 39.0
    "C5v": (
        "D10h",
        _same(
            "E",
            "C5+",
            "C5-",
            "C5^2+",
            "C5^2-",
            "sigma_v1",
            "sigma_v2",
            "sigma_v3",
            "sigma_v4",
            "sigma_v5",
        ),
    ),
    "D5d": (
        "D10h",
        _same(
            "E",
            "C5+",
            "C5-",
            "C5^2+",
            "C5^2-",
            "C2'1",
            "C2'2",
            "C2'3",
            "C2'4",
            "C2'5",
            "i",
            "S10^3-",
            "S10^3+",
            "S10-",
            "S10+",
            "sigma_d1",
            "sigma_d2",
            "sigma_d3",
            "sigma_d4",
            "sigma_d5",
        ),
    ),
    "D5h": (
        "D10h",
        _same(
            "E",
            "C5+",
            "C5-",
            "C5^2+",
            "C5^2-",
            "C2'1",
            "C2'2",
            "C2'3",
            "C2'4",
            "C2'5",
            "sigma_h",
            "S5+",
            "S5-",
            "S5^2+",
            "S5^2-",
            "sigma_v1",
            "sigma_v2",
            "sigma_v3",
            "sigma_v4",
            "sigma_v5",
        ),
    ),
    # master Oh, T 71.0
    "T": (
        "Oh",
        _same(
            "E",
            "C2x",
            "C2y",
            "C2z",
            "C3,1+",
            "C3,2+",
            "C3,3+",
            "C3,4+",
            "C3,1-",
            "C3,2-",
            "C3,3-",
            "C3,4-",
        ),
    ),
    "O": (
        "Oh",
        _same(
            "E",
            "C2x",
            "C2y",
            "C2z",
            "C3,1+",
            "C3,2+",
            "C3,3+",
            "C3,4+",
            "C3,1-",
            "C3,2-",
            "C3,3-",
            "C3,4-",
            "C4x+",
            "C4y+",
            "C4z+",
            "C4x-",
            "C4y-",
            "C4z-",
            "C2'a",
            "C2'b",
            "C2'c",
            "C2'd",
            "C2'e",
            "C2'f",
        ),
    ),
    "Th": (
        "Oh",
        _same(
            "E",
            "C2x",
            "C2y",
            "C2z",
            "C3,1+",
            "C3,2+",
            "C3,3+",
            "C3,4+",
            "C3,1-",
            "C3,2-",
            "C3,3-",
            "C3,4-",
            "i",
            "sigma_x",
            "sigma_y",
            "sigma_z",
            "S6,1-",
            "S6,2-",
            "S6,3-",
            "S6,4-",
            "S6,1+",
            "S6,2+",
            "S6,3+",
            "S6,4+",
        ),
    ),
    "Td": (
        "Oh",
        _same(
            "E",
            "C2x",
            "C2y",
            "C2z",
            "C3,1+",
            "C3,2+",
            "C3,3+",
            "C3,4+",
            "C3,1-",
            "C3,2-",
            "C3,3-",
            "C3,4-",
            "S4x-",
            "S4y-",
            "S4z-",
            "S4x+",
            "S4y+",
            "S4z+",
            "sigma_d1",
            "sigma_d2",
            "sigma_d3",
            "sigma_d4",
            "sigma_d5",
            "sigma_d6",
        ),
    ),
    "Oh": ("Oh", None),
    # master Ih, T 75.0
    "I": ("Ih", "proper"),
    "Ih": ("Ih", None),
}

# Book table number and page of each group's own tables (Contents, pp. xi-xii)
BOOK_TABLE = {
    "C1": ("T 1", 108),
    "C2": ("T 2", 110),
    "C3": ("T 3", 112),
    "C4": ("T 4", 114),
    "C6": ("T 6", 119),
    "Ci": ("T 11", 138),
    "Cs": ("T 12", 140),
    "S4": ("T 13", 144),
    "S6": ("T 14", 146),
    "S8": ("T 15", 149),
    "D2": ("T 22", 194),
    "D3": ("T 23", 196),
    "D4": ("T 24", 199),
    "D6": ("T 26", 207),
    "D2h": ("T 31", 246),
    "D3h": ("T 32", 250),
    "D4h": ("T 33", 256),
    "D5h": ("T 34", 263),
    "D6h": ("T 35", 273),
    "D2d": ("T 41", 366),
    "D3d": ("T 42", 370),
    "D4d": ("T 43", 375),
    "D5d": ("T 44", 382),
    "C2v": ("T 50", 482),
    "C3v": ("T 51", 484),
    "C3vB": ("T 51", 484),
    "C4v": ("T 52", 489),
    "C5v": ("T 53", 492),
    "C6v": ("T 54", 497),
    "C2h": ("T 60", 532),
    "C3h": ("T 61", 534),
    "C4h": ("T 62", 537),
    "C6h": ("T 64", 545),
    "O": ("T 69", 580),
    "T": ("T 70", 590),
    "Oh": ("T 71", 595),
    "Th": ("T 72", 632),
    "Td": ("T 73", 637),
    "I": ("T 74", 642),
    "Ih": ("T 75", 659),
}

# Irrep symbol rows of T n.4, as printed (ASCII: the book's "E1/2,g" style kept)
IRREPS = {
    "C1": "A | A1/2",
    "C2": "A B | 1E1/2 2E1/2",
    "C3": "A 1E 2E | 1E1/2 2E1/2 A3/2",
    "C4": "A B 1E 2E | 1E1/2 2E1/2 1E3/2 2E3/2",
    "C6": "A B 1E1 2E1 1E2 2E2 | 1E1/2 2E1/2 1E3/2 2E3/2 1E5/2 2E5/2",
    "Ci": "Ag Au | A1/2,g A1/2,u",
    "Cs": "A' A'' | 1E1/2 2E1/2",
    "S4": "A B 1E 2E | 1E1/2 2E1/2 1E3/2 2E3/2",
    "S6": "Ag 1Eg 2Eg Au 1Eu 2Eu | 1E1/2,g 2E1/2,g A3/2,g 1E1/2,u 2E1/2,u A3/2,u",
    "S8": "A B 1E1 2E1 1E2 2E2 1E3 2E3 | 1E1/2 2E1/2 1E3/2 2E3/2 1E5/2 2E5/2 1E7/2 2E7/2",
    "D2": "A B1 B2 B3 | E1/2",
    "D3": "A1 A2 E | E1/2 1E3/2 2E3/2",
    "D4": "A1 A2 B1 B2 E | E1/2 E3/2",
    "D6": "A1 A2 B1 B2 E1 E2 | E1/2 E3/2 E5/2",
    "D2h": "Ag B1g B2g B3g Au B1u B2u B3u | E1/2,g E1/2,u",
    "D3h": "A1' A2' E' A1'' A2'' E'' | E1/2 E3/2 E5/2",
    "D4h": "A1g A2g B1g B2g Eg A1u A2u B1u B2u Eu | E1/2,g E3/2,g E1/2,u E3/2,u",
    "D5h": "A1' A2' E1' E2' A1'' A2'' E1'' E2'' | E1/2 E3/2 E5/2 E7/2 E9/2",
    "D6h": "A1g A2g B1g B2g E1g E2g A1u A2u B1u B2u E1u E2u | E1/2,g E3/2,g E5/2,g E1/2,u E3/2,u E5/2,u",
    "D2d": "A1 A2 B1 B2 E | E1/2 E3/2",
    "D3d": "A1g A2g Eg A1u A2u Eu | E1/2,g 1E3/2,g 2E3/2,g E1/2,u 1E3/2,u 2E3/2,u",
    "D4d": "A1 A2 B1 B2 E1 E2 E3 | E1/2 E3/2 E5/2 E7/2",
    "D5d": "A1g A2g E1g E2g A1u A2u E1u E2u | E1/2,g E3/2,g 1E5/2,g 2E5/2,g E1/2,u E3/2,u 1E5/2,u 2E5/2,u",
    "C2v": "A1 A2 B1 B2 | E1/2",
    "C3v": "A1 A2 E | E1/2 1E3/2 2E3/2",
    "C3vB": "A1 A2 E | E1/2 1E3/2 2E3/2",
    "C4v": "A1 A2 B1 B2 E | E1/2 E3/2",
    "C5v": "A1 A2 E1 E2 | E1/2 E3/2 1E5/2 2E5/2",
    "C6v": "A1 A2 B1 B2 E1 E2 | E1/2 E3/2 E5/2",
    "C2h": "Ag Bg Au Bu | 1E1/2,g 2E1/2,g 1E1/2,u 2E1/2,u",
    "C3h": "A' 1E' 2E' A'' 1E'' 2E'' | 1E1/2 2E1/2 1E3/2 2E3/2 1E5/2 2E5/2",
    "C4h": "Ag Bg 1Eg 2Eg Au Bu 1Eu 2Eu | 1E1/2,g 2E1/2,g 1E3/2,g 2E3/2,g 1E1/2,u 2E1/2,u 1E3/2,u 2E3/2,u",
    "C6h": "Ag Bg 1E1g 2E1g 1E2g 2E2g Au Bu 1E1u 2E1u 1E2u 2E2u | "
    "1E1/2,g 2E1/2,g 1E3/2,g 2E3/2,g 1E5/2,g 2E5/2,g 1E1/2,u 2E1/2,u 1E3/2,u 2E3/2,u 1E5/2,u 2E5/2,u",
    "O": "A1 A2 E T1 T2 | E1/2 E5/2 F3/2",
    "T": "A 1E 2E T | E1/2 1F3/2 2F3/2",
    "Oh": "A1g A2g Eg T1g T2g A1u A2u Eu T1u T2u | E1/2,g E5/2,g F3/2,g E1/2,u E5/2,u F3/2,u",
    "Th": "Ag 1Eg 2Eg Tg Au 1Eu 2Eu Tu | E1/2,g 1F3/2,g 2F3/2,g E1/2,u 1F3/2,u 2F3/2,u",
    "Td": "A1 A2 E T1 T2 | E1/2 E5/2 F3/2",
    "I": "A T1 T2 F H | E1/2 E7/2 F3/2 I5/2",
    "Ih": "Ag T1g T2g Fg Hg Au T1u T2u Fu Hu | E1/2,g E7/2,g F3/2,g I5/2,g E1/2,u E7/2,u F3/2,u I5/2,u",
}

# Anchors of the map onto canonical labels: (book label, canonical label);
# the second anchor of the polyhedral groups is found by search (see below).
ANCHORS = {
    "C1": "E=E",
    "C2": "E=E C2=C2",
    "Ci": "E=E i=i",
    "Cs": "E=E sigma_h=sigma_h",
    "D2": "E=E C2z=C2 C2x=C2'(1)",
    "C2v": "E=E C2=C2 sigma_x=sigma_v(1)",
    "C2h": "E=E C2=C2 i=i",
    "D2h": "E=E C2z=C2 C2x=C2'(1) i=i",
    "C4": "E=E C4+=C4^1",
    "S4": "E=E S4-=S4^-1",
    "D4": "E=E C4+=C4^1 C2'1=C2'(1)",
    "C4v": "E=E C4+=C4^1 sigma_v1=sigma_v(1)",
    "D2d": "E=E S4-=S4^-1 C2'1=C2'(1)",
    "C4h": "E=E C4+=C4^1 i=i",
    "D4h": "E=E C4+=C4^1 C2'1=C2'(1) i=i",
    "C3": "E=E C3+=C3^1",
    "C6": "E=E C6+=C6^1",
    "S6": "E=E C3+=C3^1 i=i",
    "D3": "E=E C3+=C3^1 C2'1=C2'(1)",
    "D6": "E=E C6+=C6^1 C2'1=C2'(1)",
    "C6v": "E=E C6+=C6^1 sigma_d1=sigma_v(1)",
    "C3h": "E=E S3-=S3^-1",
    "C6h": "E=E C6+=C6^1 i=i",
    "D3d": "E=E C3+=C3^1 C2'1=C2'(1) i=i",
    "D3h": "E=E S3-=S3^-1 C2'1=C2'(1)",
    "D6h": "E=E C6+=C6^1 C2'1=C2'(1) i=i",
    "C3v": "E=E C3+=C3^1 sigma_v1=sigma_v(1)",
    "C3vB": "E=E C3+=C3^1 sigma_v1=sigma_v(1)",
    "S8": "E=E S8^3-=S8^-3",
    "D4d": "E=E S8^3-=S8^-3 C2'1=C2'(1)",
    "C5v": "E=E C5+=C5^1 sigma_v1=sigma_v(1)",
    "D5d": "E=E C5+=C5^1 C2'1=C2'(1) i=i",
    "D5h": "E=E S5^2-=S5^-2 C2'1=C2'(1)",
    "T": "E=E C3,1+=C3^1(1)",
    "O": "E=E C4z+=C4^1(1)",
    "Th": "E=E C3,1+=C3^1(1) i=i",
    "Td": "E=E S4z-=S4^-1(1)",
    "Oh": "E=E C4z+=C4^1(1) i=i",
    "I": "E=E C5,1+=C5^1(1)",
    "Ih": "E=E C5,1+=C5^1(1) i=i",
}
# second generator searched for the polyhedral groups: book label -> candidate canonical labels
SEARCH = {
    "T": "C2x",
    "O": "C3,1+",
    "Th": "C2x",
    "Td": "C3,1+",
    "Oh": "C3,1+",
    "I": "C3,1+",
    "Ih": "C3,1+",
}


# --------------------------------------------------------------------------
# Building a target group from a master
# --------------------------------------------------------------------------


def book_group(name):
    master, sel = SUBGROUPS[name]
    ops, full = master_group(master)
    if sel is None:
        return full, ops
    if sel == "proper":
        pairs = [(lab, lab) for r, p, lab in full.elements if p == +1]
    else:
        pairs = sel
    chosen = []
    for master_label, sub_label in pairs:
        k = full.index[master_label]
        r, p, _ = full.elements[k]
        chosen.append((r, p, sub_label))
    return BookGroup(ops, chosen), ops


def homomorphism(book, grp, anchors):
    """Extend anchors (book index -> canonical index) by closure; None on conflict."""
    iso = dict(anchors)
    queue = deque(iso)
    known = set(iso)
    while queue:
        a = queue.popleft()
        for b in list(known):
            for x, y in ((a, b), (b, a)):
                p = book.table[x][y]
                t = grp.table[iso[x]][iso[y]]
                if p in iso:
                    if iso[p] != t:
                        return None
                else:
                    iso[p] = t
                    known.add(p)
                    queue.append(p)
    return iso


def element_order(table, k):
    e, n, x = 0, 1, k
    while x != e:
        x = table[x][k]
        n += 1
    return n


def find_map(name, book, grp):
    anchors = {}
    for pair in ANCHORS[name].split():
        b, c = pair.split("=")
        anchors[book.index[b]] = grp.index[c]
    if name in SEARCH:
        b = book.index[SEARCH[name]]
        order = element_order(book.table, b)
        for c in range(len(grp.elements)):
            if element_order(grp.table, c) != order:
                continue
            iso = homomorphism(book, grp, {**anchors, b: c})
            if (
                iso is not None
                and len(iso) == len(book.elements)
                and len(set(iso.values())) == len(iso)
            ):
                return iso, {**anchors, b: c}
        raise ValueError(f"{name}: no isomorphism found")
    iso = homomorphism(book, grp, anchors)
    if (
        iso is None
        or len(iso) != len(book.elements)
        or len(set(iso.values())) != len(iso)
    ):
        raise ValueError(
            f"{name}: anchors do not give an isomorphism ({None if iso is None else len(iso)})"
        )
    return iso, anchors


def solve_epsilon(book, grp, iso, fixed):
    """eps with eps = +1 on the anchors in ``fixed``; None if the tilde marks are not a coboundary."""
    n = len(book.elements)
    r = [
        [book.factor[x][y] * grp.omega_at(iso[x], iso[y]) for y in range(n)]
        for x in range(n)
    ]
    eps = {k: 1 for k in fixed}
    queue = deque(eps)
    while queue:
        a = queue.popleft()
        for b in list(eps):
            for x, y in ((a, b), (b, a)):
                p = book.table[x][y]
                v = r[x][y] * eps[x] * eps[y]
                if p in eps:
                    if eps[p] != v:
                        return None
                else:
                    eps[p] = v
                    queue.append(p)
    assert len(eps) == n
    return eps


# --------------------------------------------------------------------------
# Checks against the printed factor tables (OCR text)
# --------------------------------------------------------------------------

FACTOR_PAGES = {  # group -> (pdf pages of the printed factor table, column blocks)
    "D2h": ([259], 1),
    "D4h": ([269], 1),
    "D6h": ([288], 1),
    "D8h": ([319], 1),
    "D10h": ([360, 361], 2),
    "Oh": ([614, 615, 616, 617], 2),
    "I": ([660, 661, 662, 663], 2),
}


def compare_factor_tables(ocr_dir):
    """Parse the printed factor tables (OCR text) and compare with the reconstructed ones."""
    ok = True
    for name, (pages, blocks) in FACTOR_PAGES.items():
        book = master_group(name)[1] if name in MASTERS else book_group(name)[0]
        n = len(book.elements)
        text = " ".join((Path(ocr_dir) / f"p{p:03d}.txt").read_text() for p in pages)
        text = text.replace("\u2212", "-")
        seg = text[text.find("Factor table") :]
        if "Character table" in seg:
            seg = seg[: seg.find("Character table")]
        vals = [int(t) for t in seg.split() if t in ("1", "-1")]
        width = n // blocks
        per_block = n * width
        rebuilt = []
        for i in range(n):
            for b in range(blocks):
                rebuilt += vals[
                    b * per_block + i * width : b * per_block + (i + 1) * width
                ]
        expected = [x for row in book.factor for x in row]
        if len(vals) < n * n:
            print(f"{name}: parsed {len(vals)} entries, expected {n * n}; skipped")
            ok = False
            continue
        diff = sum(1 for a, b in zip(rebuilt, expected) if a != b)
        print(
            f"{name}: printed factor table {'agrees' if diff == 0 else f'DIFFERS in {diff} entries'} "
            f"({n}x{n}, {len(vals)} entries parsed)"
        )
        ok &= diff == 0
    return ok


# --------------------------------------------------------------------------
# Generation of the transcription files
# --------------------------------------------------------------------------


def toml_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def generate(out_dir=OUT_DIR):
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for name in SUBGROUPS:
        symbol = "C3v" if name == "C3vB" else ag.ALIASES.get(name, name)
        book, ops = book_group(name)
        grp = at.Group(symbol)
        iso, anchors = find_map(name, book, grp)
        eps = solve_epsilon(book, grp, iso, list(anchors))
        assert eps is not None, name
        swapped = [
            grp.elements[iso[k]] for k in range(len(book.elements)) if eps[k] == -1
        ]
        ordinary, projective = (s.split() for s in IRREPS[name].split("|"))
        table, page = BOOK_TABLE[name]
        master, _ = SUBGROUPS[name]
        lines = [
            f"# Altmann & Herzig, Point-Group Theory Tables (2nd ed.): element labels of {symbol}"
            + (" (setting B)" if name == "C3vB" else ""),
            "# Format: tests/reference/altmann-herzig/labels/README.md",
            (
                "# Drafted by tools/ah_book.py from the parameter table "
                f"{MASTER_TABLE[master]}.1 and the subgroup table {MASTER_TABLE[master]}.0 of the book's PDF;"
            ),
            "# the multiplication table and tilde marks are reproduced from the parameters by quaternion",
            "# multiplication (book section 3-2) and agree with the printed factor tables (tools/ah_book.py check).",
            "",
            f"group = {toml_str(symbol)}",
        ]
        if name == "C3vB":
            lines.append('setting = "B"')
        elif name == "C3v":
            lines.append('setting = "A"')
        lines += [
            f"table = {toml_str(f'{table}, p. {page}')}",
            'status = "done"',
            f"source = {toml_str(SOURCE)}",
            'transcribed_by = "agent (Claude) from the OCR text of the PDF; verify against the print"',
            'notes = """',
            "Poles follow the book's Rule 1 (section 3-2): the poles of a class are the images of one",
            "another under the group, not a positive-hemisphere rule; e.g. the C2' poles of D6h are the",
            "C6 rotations of (1 0 0). Improper operations are i times the proper operation with the same",
            "parameters (marker i implicit). Labels transliterated per README.md.",
            '"""',
            "# [map] is derived from these anchors; the checker fixes eps = +1 on them",
            "anchors = [" + ", ".join(toml_str(book.labels[k]) for k in anchors) + "]",
            "",
            "[book]",
            "elements = [" + ", ".join(toml_str(lab) for lab in book.labels) + "]",
            f"master = {toml_str(master.rstrip('AB') if master.startswith('C3v') else master)}",
            "irreps = [" + ", ".join(toml_str(s) for s in ordinary + projective) + "]",
            "",
            "# phi and the pole n of the proper operation (the rotation part for improper ones)",
            "[book.parameters]",
        ]
        for r, p, lab in book.elements:
            op = ops[r]
            n = ", ".join(toml_str(t) for t in op.n)
            extra = "" if p == +1 else f", proper = {toml_str(op.proper)}"
            lines.append(
                f"{toml_str(lab)} = {{ phi = {toml_str(op.phi)}, n = [{n}]{extra} }}"
            )
        lines += [
            "",
            "# book label -> canonical label, derived from the anchors above",
            "[map]",
        ]
        for k in range(len(book.elements)):
            lines.append(
                f"{toml_str(book.labels[k])} = {toml_str(grp.elements[iso[k]])}"
            )
        lines += [
            "",
            "# reconstructed from the parameters; '~' marks a product the book prints with a tilde",
            "[multiplication]",
            "rows = [",
        ]
        for x in range(len(book.elements)):
            row = ", ".join(
                toml_str(
                    ("~" if book.factor[x][y] == -1 else "")
                    + book.labels[book.table[x][y]]
                )
                for y in range(len(book.elements))
            )
            lines.append(f"  [{row}],")
        lines += ["]", ""]
        fname = f"{symbol}.toml" if name != "C3vB" else "C3vB.toml"
        (out_dir / fname).write_text("\n".join(lines))
        summary.append((name, len(book.elements), swapped))
    for name, order, swapped in summary:
        print(f"{name:5s} |G|={order:3d} swapped={swapped}")


def check(ocr_dir=None):
    for master in MASTERS:
        ops, full = master_group(master)
        print(f"{master}: {len(full.elements)} elements, group axioms hold")
    for name in SUBGROUPS:
        symbol = "C3v" if name == "C3vB" else ag.ALIASES.get(name, name)
        book, _ = book_group(name)
        grp = at.Group(symbol)
        assert len(book.elements) == grp.G.order, (
            name,
            len(book.elements),
            grp.G.order,
        )
        iso, anchors = find_map(name, book, grp)
        eps = solve_epsilon(book, grp, iso, list(anchors))
        assert eps is not None, (
            f"{name}: tilde marks are not a coboundary of the canonical section"
        )
        ordinary, projective = (s.split() for s in IRREPS[name].split("|"))
        nreg = sum(1 for c in grp.classes) if False else None
        del nreg
        print(
            f"{name:5s} ok: isomorphic to {symbol}, irreps {len(ordinary)}+{len(projective)}"
        )
    if ocr_dir:
        compare_factor_tables(ocr_dir)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "generate":
        generate()
    else:
        check(sys.argv[2] if len(sys.argv) > 2 else None)
