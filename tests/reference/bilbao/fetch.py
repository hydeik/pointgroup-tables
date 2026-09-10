"""Fetch the Bilbao character tables of the 32 crystallographic point groups
(wayfinder ticket #20).

Two programs of the Bilbao Crystallographic Server feed one pair of reference
files per group:

* **POINT** (``/cgi-bin/rep/programs/sam/point.py?sg=<sg>&num=<n>``): the
  character table of the point group with the server's class headers
  (Seitz-like symbols ``3+``, ``2_100``, ``m_110``), the Mulliken and the
  Bradley & Cracknell ``Gamma`` label of every irrep;
* **REPRESENTATIONS DPG** (``/cgi-bin/cryst/programs/representations_out.pl
  ?tipogrupo=dbg&pointspace=point&num=<sg>&super=<n>``): the character table
  of the double point group with its classes ``C1 .. Ck`` listed element by
  element (``d`` marks the operations that differ by a 2 pi spin rotation),
  three labels per irrep (Bradley & Cracknell ``Gamma``, Mulliken, CDML
  ``Gamma``), the reality of every irrep, and the 3 x 3 and SU(2) matrices of
  every operation.

Both are written as ``characters/<group>.toml`` (``double = false``) and
``characters/<group>-double.toml`` (``double = true``) in the format of
``docs/spec/verification.md`` ("Reference data"). The script is run once; the
tests read the TOML only, and the script stays as the record of provenance.

    uv run python tests/reference/bilbao/fetch.py [--raw-dir DIR] [--turnstile-cookie VALUE] [GROUP ...]

``--raw-dir`` caches the fetched HTML (a page present there is not fetched
again). The double-group program sits behind the server's Cloudflare
Turnstile check ("please verify you are human"); the check sets a cookie
``turnstile_passed=<unix time>`` for an hour, and the program refuses (HTTP
500) without it. Pass the challenge once in a browser and give the cookie's
value with ``--turnstile-cookie``; the POINT pages need no cookie.

Values are converted into the value grammar of
``docs/spec/basis-conventions.md`` by the exact field arithmetic of the
basis-convention prototype; the printed expression is kept next to it as
``server_values``. Element orders of the double group are computed from the
server's own matrices (order of the 3 x 3 part and of the SU(2) part) and
checked against the Seitz symbol. Nothing here is production code.
"""

from __future__ import annotations

import argparse
import cmath
import datetime as dt
import html
import math
import os
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "prototypes"))
import canonical_basis_prototype as cb

HOST = "https://www.cryst.ehu.eus"
POINT_URL = HOST + "/cgi-bin/rep/programs/sam/point.py?sg={sg}&num={num}"
DPG_URL = (
    HOST
    + "/cgi-bin/cryst/programs/representations_out.pl?tipogrupo=dbg&pointspace=point&num={sg}&super={num}"
)
OUT_DIR = Path(__file__).resolve().parent / "characters"
USER_AGENT = "pointgroup-tables reference fetch (github.com/hydeik/pointgroup-tables)"
DELAY = 1.5  # seconds between requests

SOURCE_POINT = (
    "Bilbao Crystallographic Server, Point Group Tables (POINT): "
    "M. I. Aroyo, A. Kirov, C. Capillas, J. M. Perez-Mato and H. Wondratschek, "
    "Acta Cryst. A62 (2006) 115-128"
)
SOURCE_DPG = (
    "Bilbao Crystallographic Server, REPRESENTATIONS DPG: "
    "L. Elcoro, B. Bradlyn, Z. Wang, M. G. Vergniory, J. Cano, C. Felser, B. A. Bernevig, "
    "D. Orobengoa, G. de la Flor and M. I. Aroyo, J. Appl. Cryst. 50 (2017) 1457-1477"
)

# group symbol -> (Bilbao point-group number, space group of the POINT link,
#                  space group of the REPRESENTATIONS DPG link, Hermann-Mauguin symbol)
# The two index pages link C3v to different space groups (P3m1 and P31m); each
# file cites the page its own index links to.
GROUPS = {
    "C1": (1, 1, 1, "1"),
    "S2": (2, 2, 2, "-1"),
    "C2": (3, 3, 3, "2"),
    "C1h": (4, 6, 6, "m"),
    "C2h": (5, 10, 10, "2/m"),
    "D2": (6, 16, 16, "222"),
    "C2v": (7, 25, 25, "mm2"),
    "D2h": (8, 47, 47, "mmm"),
    "C4": (9, 75, 75, "4"),
    "S4": (10, 81, 81, "-4"),
    "C4h": (11, 83, 83, "4/m"),
    "D4": (12, 89, 89, "422"),
    "C4v": (13, 99, 99, "4mm"),
    "D2d": (14, 111, 111, "-42m"),
    "D4h": (15, 123, 123, "4/mmm"),
    "C3": (16, 143, 143, "3"),
    "S6": (17, 147, 147, "-3"),
    "D3": (18, 149, 149, "32"),
    "C3v": (19, 156, 157, "3m"),
    "D3d": (20, 162, 162, "-3m"),
    "C6": (21, 168, 168, "6"),
    "C3h": (22, 174, 174, "-6"),
    "C6h": (23, 175, 175, "6/m"),
    "D6": (24, 177, 177, "622"),
    "C6v": (25, 183, 183, "6mm"),
    "D3h": (26, 187, 187, "-6m2"),
    "D6h": (27, 191, 191, "6/mmm"),
    "T": (28, 195, 195, "23"),
    "Th": (29, 200, 200, "m-3"),
    "O": (30, 207, 207, "432"),
    "Td": (31, 215, 215, "-43m"),
    "Oh": (32, 221, 221, "m-3m"),
}
ALIASES = {"Ci": "S2", "Cs": "C1h", "C3i": "S6"}


# --------------------------------------------------------------------------
# Fetching
# --------------------------------------------------------------------------


def fetch(url, raw_path, cookie):
    if raw_path is not None and raw_path.exists():
        return raw_path.read_bytes()
    headers = {"User-Agent": USER_AGENT}
    if cookie:
        headers["Cookie"] = f"turnstile_passed={cookie}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=120) as resp:
        data = resp.read()
    if raw_path is not None:
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_bytes(data)
    time.sleep(DELAY)
    return data


# --------------------------------------------------------------------------
# HTML helpers
# --------------------------------------------------------------------------

OVERLINE_FONT = re.compile(
    r'<font style=\s*"text-decoration:\s*overline;">(.*?)</font>', re.DOTALL
)
RADIC = re.compile(
    r'(?:&radic;|√)\s*<span style="text-decoration:\s*overline;">(\d+)</span>',
    re.DOTALL,
)


def marked(s):
    """HTML fragment -> text with ^{..}, _{..}, bar{..} and sqrtN kept."""
    s = RADIC.sub(r"sqrt\1", s)
    s = OVERLINE_FONT.sub(r"bar{\1}", s)
    s = re.sub(r"<sup>(.*?)</sup>", r"^{\1}", s, flags=re.DOTALL)
    s = re.sub(r"<sub>(.*?)</sub>", r"_{\1}", s, flags=re.DOTALL)
    s = re.sub(r"<[^>]+>", "", s)
    s = html.unescape(s).replace("\xa0", " ")
    return re.sub(r"\s+", " ", s).strip()


def seitz_marked(s):
    """Seitz symbol fragment -> marked text; the matrices table typesets the
    minus signs of rotoinversions and of direction indices as overlines."""
    return marked(OVERLINE_FONT.sub(r"-\1", s))


def cells(row_html):
    return re.findall(r"<td[^>]*>(.*?)</td>", row_html, flags=re.DOTALL)


def rows(table_html):
    return re.findall(r"<tr>(.*?)</tr>", table_html, flags=re.DOTALL)


# --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------

LABEL_RE = re.compile(
    r"^(?:\^\{(?P<pre>\d)\})?(?P<bar>bar\{)?(?P<base>[A-Z]+|Γ)\}?"
    r"(?:_\{(?P<sub>[^}]*)\})?(?:\^\{(?P<sup>[^}]*)\})?(?:_\{(?P<sub2>[^}]*)\})?$"
)


def label_ascii(text):
    """Server irrep label (marked text) -> (ASCII label, barred).

    The overline of a double-valued irrep becomes ``~`` before the barred
    token, as the ``~`` of a barred element: ``^{1}bar{E}_{1g}`` -> ("1~E1g",
    True); ``A^{''}_{1}`` -> ("A1''", False); ``GM_{5}^{+}`` and ``Γ_{5}^{+}``
    -> ("Gamma5+", False); ``bar{GM}_{7}`` -> ("~Gamma7", True).
    """
    m = LABEL_RE.match(text)
    assert m, text
    base = m.group("base")
    if base in ("GM", "Γ"):
        base = "Gamma"
    if m.group("bar"):
        base = "~" + base
    sub = (m.group("sub") or "") + (m.group("sub2") or "")
    sup = m.group("sup") or ""
    primes = "".join(ch for ch in sub if ch == "'") + sup
    digits = "".join(ch for ch in sub if ch != "'")
    return (m.group("pre") or "") + base + digits + primes, bool(m.group("bar"))


# --------------------------------------------------------------------------
# Values: the server's expressions -> exact field elements / complex numbers
# --------------------------------------------------------------------------

TOKEN = re.compile(r"\s*(?:(1j)|(\d+)|(sqrt\d+)|(e\^\{[^}]*\})|(i)|(w)|([()+\-*/^]))")


def tokenize(text):
    text = re.sub(
        r"\^\{(\d+)\}", r"^\1", text
    )  # w^{2} -> w^2 (e^{..} keeps its braces)
    pos, out = 0, []
    while pos < len(text):
        m = TOKEN.match(text, pos)
        assert m and m.end() > pos, (text, pos)
        pos = m.end()
        tok = next(g for g in m.groups() if g is not None)
        out.append(tok)
    return out


class Exact:
    """Arithmetic in Q(zeta_M) through the prototype field."""

    def __init__(self, F):
        self.F = F

    def integer(self, n):
        return self.F.rational(n)

    def i(self):
        return self.F.i()

    def sqrt(self, n):
        if n == 2:
            return self.F.sqrt2()
        assert n == 3, n
        return cb.const(self.F, "sqrt3")

    def w(self):  # the POINT pages' w = exp(2 pi i / 3)
        return self.F.zeta(self.F.conductor // 3)

    def exp(self, k, n):  # exp(i k pi / n)
        M = self.F.conductor
        assert (2 * n) and M % (2 * n) == 0, (k, n, M)
        return self.F.zeta((M // (2 * n)) * k)

    def add(self, x, y):
        return self.F.add(x, y)

    def sub(self, x, y):
        return self.F.sub(x, y)

    def mul(self, x, y):
        return self.F.mul(x, y)

    def div(self, x, y):
        return self.F.mul(x, cb.finv(self.F, y))

    def neg(self, x):
        return self.F.neg(x)

    def pow(self, x, n):
        out = self.F.rational(1)
        for _ in range(n):
            out = self.F.mul(out, x)
        return out


class Numeric:
    def integer(self, n):
        return complex(n)

    def i(self):
        return 1j

    def sqrt(self, n):
        return complex(math.sqrt(n))

    def w(self):
        return cmath.exp(2j * math.pi / 3)

    def exp(self, k, n):
        return cmath.exp(1j * math.pi * k / n)

    add = staticmethod(lambda x, y: x + y)
    sub = staticmethod(lambda x, y: x - y)
    mul = staticmethod(lambda x, y: x * y)
    div = staticmethod(lambda x, y: x / y)
    neg = staticmethod(lambda x: -x)
    pow = staticmethod(lambda x, n: x**n)


EXP_RE = re.compile(r"^e\^\{(-)?i(\d*)π(?:/(\d+))?\}$")


class Parser:
    """expr := term (('+'|'-') term)*; term := factor (('*'|'/')? factor)*;
    factor := '-'? atom ('^' int)?; atom := int | i | 1j | w | sqrtN | e^{..} | '(' expr ')'.
    Juxtaposition multiplies (``i sqrt3``, ``2 i``)."""

    def __init__(self, arith, text):
        self.A = arith
        self.toks = tokenize(text)
        self.pos = 0

    def peek(self):
        return self.toks[self.pos] if self.pos < len(self.toks) else None

    def take(self):
        tok = self.toks[self.pos]
        self.pos += 1
        return tok

    def parse(self):
        x = self.expr()
        assert self.peek() is None, self.toks
        return x

    def expr(self):
        x = self.term()
        while self.peek() in ("+", "-"):
            if self.take() == "+":
                x = self.A.add(x, self.term())
            else:
                x = self.A.sub(x, self.term())
        return x

    def term(self):
        x = self.factor()
        while True:
            tok = self.peek()
            if tok == "*":
                self.take()
                x = self.A.mul(x, self.factor())
            elif tok == "/":
                self.take()
                x = self.A.div(x, self.factor())
            elif tok is not None and tok not in ("+", "-", ")", "^"):
                x = self.A.mul(x, self.factor())
            else:
                return x

    def factor(self):
        if self.peek() == "-":
            self.take()
            return self.A.neg(self.factor())
        x = self.atom()
        if self.peek() == "^":
            self.take()
            x = self.A.pow(x, int(self.take()))
        return x

    def atom(self):
        tok = self.take()
        if tok == "(":
            x = self.expr()
            assert self.take() == ")"
            return x
        if tok.isdigit():
            return self.A.integer(int(tok))
        if tok in ("i", "1j"):
            return self.A.i()
        if tok == "w":
            return self.A.w()
        if tok.startswith("sqrt"):
            return self.A.sqrt(int(tok[4:]))
        m = EXP_RE.match(tok)
        assert m, tok
        k = int(m.group(2) or 1) * (-1 if m.group(1) else 1)
        return self.A.exp(k, int(m.group(3) or 1))


def conductor_need(text):
    """Smallest conductor whose field holds every symbol of the expression."""
    need = 1
    for tok in tokenize(text):
        if tok in ("i", "1j"):
            need = math.lcm(need, 4)
        elif tok == "w":
            need = math.lcm(need, 3)
        elif tok == "sqrt2":
            need = math.lcm(need, 8)
        elif tok == "sqrt3":
            need = math.lcm(need, 12)
        elif tok.startswith("e^"):
            m = EXP_RE.match(tok)
            assert m, tok
            need = math.lcm(need, 2 * int(m.group(3) or 1))
    return need


def exact_value(F, text):
    return Parser(Exact(F), text).parse()


def numeric_value(text):
    return Parser(Numeric(), text).parse()


def server_ascii(text):
    """The printed expression in the value grammar's spelling (``i*sqrt3``)."""
    toks = tokenize(text)
    out = []
    for k, tok in enumerate(toks):
        if tok == "1j":
            tok = "i"
        if k and out[-1] not in "(+-*/^" and tok not in ")+-*/^":
            out.append("*")
        out.append(tok)
    return "".join(out)


# --------------------------------------------------------------------------
# Seitz symbols -> invariants
# --------------------------------------------------------------------------

SEITZ_RE = re.compile(
    r"^(?:\^\{d\})?(?P<d>d)?(?P<imp>-)?(?P<body>[12346]|m)(?:\^\{(?P<sense>[+-])\})?(?:_\{(?P<dir>[^}]*)\})?$"
)


def seitz_ascii(text):
    """``^{d}3^{+}_{001}`` -> ``d3+_001``; ``m_{1-10}`` -> ``m_1-10``; ``^{d}1_{}`` -> ``d1``."""
    m = SEITZ_RE.match(text.replace("^{d}", "d"))
    assert m, text
    out = (
        ("d" if m.group("d") else "")
        + ("-" if m.group("imp") else "")
        + m.group("body")
    )
    out += m.group("sense") or ""
    if m.group("dir"):
        out += "_" + m.group("dir")
    return out


def seitz_invariants(text):
    """(order in the double group, parity, barred) from the Seitz symbol alone;
    used as a cross-check of the matrices."""
    m = SEITZ_RE.match(text.replace("^{d}", "d"))
    assert m, text
    barred = bool(m.group("d"))
    improper = bool(m.group("imp")) or m.group("body") == "m"
    n = 2 if m.group("body") == "m" else int(m.group("body"))
    if n == 1:
        su2 = 2 if barred else 1
    elif barred:
        su2 = 2 * n // math.gcd(n + 1, 2 * n)
    else:
        su2 = 2 * n
    return (math.lcm(su2, 2) if improper else su2, -1 if improper else 1, barred)


def point_class_invariants(text):
    """POINT class header (``3^{+}``, ``-4^{-}``, ``m_{100}``) -> (order in G, parity)."""
    m = SEITZ_RE.match(text)
    assert m and not m.group("d"), text
    body = m.group("body")
    if body == "m":
        return 2, -1
    n = int(body)
    if m.group("imp"):
        return {1: 2, 2: 2, 3: 6, 4: 4, 6: 6}[n], -1
    return n, 1


# --------------------------------------------------------------------------
# POINT pages
# --------------------------------------------------------------------------


def parse_point(page):
    h = page.decode("utf-8", errors="replace")
    i = h.find("Character Table of the group")
    j = h.find("</table>", i)
    assert i > 0 and j > i, "character table not found"
    trs = rows(h[i:j])
    header = [marked(c) for c in cells(trs[0])]
    assert header[1] == "#", header
    end = -1 if header[-1] == "functions" else None  # C1 has no functions column
    class_labels = header[2:end]
    body = trs[1:]
    sizes = [1] * len(class_labels)
    if marked(cells(body[0])[0]) == "Mult.":
        sizes = [int(marked(c)) for c in cells(body[0])[2:end]]
        body = body[1:]
    classes = []
    for lab, size in zip(class_labels, sizes):
        order, parity = point_class_invariants(lab)
        classes.append(
            {"label": seitz_ascii(lab), "size": size, "order": order, "parity": parity}
        )
    irreps = []
    for tr in body:
        cs = cells(tr)
        # a complex-conjugate pair is printed as two lines in one row
        parts = [c.split("<br>") for c in cs[:end]]
        depth = len(parts[0])
        assert all(len(p) == depth for p in parts), [marked(c) for c in cs]
        for k in range(depth):
            label, _ = label_ascii(marked(parts[0][k]))
            bethe, _ = label_ascii(marked(parts[1][k]))
            values = [marked(p[k]) for p in parts[2:]]
            irreps.append({"label": label, "bethe": bethe, "values": values})
    tail = h[j : h.find("Subgroups of the group", j)]
    notes = [
        re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", n))).strip()
        for n in re.findall(
            r"\[\s*<b>Note:</b>(.*?)(?:</small>\s*\]|\])", tail, flags=re.DOTALL
        )
    ]
    w = re.search(r"w = exp\((\d+)iπ/(\d+)\)", html.unescape(tail))
    return {
        "classes": classes,
        "irreps": irreps,
        "notes": notes,
        "w": f"exp({w.group(1)} pi i / {w.group(2)})" if w else None,
    }


# --------------------------------------------------------------------------
# REPRESENTATIONS DPG pages
# --------------------------------------------------------------------------


def matrix_cells(inner):
    vals = [
        marked(c) for c in re.findall(r"<td[^>]*>(.*?)</td>", inner, flags=re.DOTALL)
    ]
    n = round(math.sqrt(len(vals)))
    assert n * n == len(vals), vals
    return [[numeric_value(vals[r * n + c]) for c in range(n)] for r in range(n)]


def mat_mul(A, B):
    n = len(A)
    return [
        [sum(A[r][k] * B[k][c] for k in range(n)) for c in range(n)] for r in range(n)
    ]


def is_identity(A, tol=1e-9):
    return all(
        abs(A[r][c] - (1 if r == c else 0)) < tol
        for r in range(len(A))
        for c in range(len(A))
    )


def matrix_order(A, limit=48):
    P = A
    for k in range(1, limit + 1):
        if is_identity(P):
            return k
        P = mat_mul(P, A)
    raise AssertionError("matrix of order > 48")


def det3(A):
    return (
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
    )


def parse_dpg(page):
    h = page.decode("latin-1")
    # character table
    i = h.find("Table of characters")
    j = h.find("</table>", i)
    trs = rows(h[i:j])
    header = [marked(c) for c in cells(trs[0])]
    assert header[:3] == ["(1)", "(2)", "(3)"], header
    class_names = header[3:]
    irreps = []
    for tr in trs[1:]:
        cs = [marked(c) for c in cells(tr)]
        bethe, bar1 = label_ascii(cs[0])
        label, bar2 = label_ascii(cs[1])
        cdml, bar3 = label_ascii(cs[2])
        assert not bar1 and bar2 == bar3, cs[:3]
        irreps.append(
            {
                "bethe": bethe,
                "label": label,
                "cdml": cdml,
                "barred": bar2,
                "values": cs[3:],
            }
        )
    # classes
    k = h.find("Lists of symmetry operations in the conjugacy classes")
    l = h.find("</table>", k)
    classes = []
    for td in re.findall(r"<td>(.*?)</td>", h[k:l], flags=re.DOTALL):
        name, _, elements = marked(td).partition(":")
        classes.append(
            {
                "label": name.replace("_{", "").replace("}", ""),
                "elements_raw": [e.strip() for e in elements.split(",")],
            }
        )
    assert [c["label"] for c in classes] == [
        n.replace("_{", "").replace("}", "") for n in class_names
    ], (classes, class_names)
    # matrices: reality per irrep from the header, order and parity per element
    m = h.find("Matrices of the representations of the group")
    section = h[m:]
    marker = re.compile(
        r'<tr><td style="border-right: 0px"><center>(\d+)</center></td>'
    )
    marks = list(marker.finditer(section))
    head = section[: marks[0].start()]
    reality = {}
    for lab, r in re.findall(
        r"<td><center>(.*?)\((-?[01])\)</center></td>", head, flags=re.DOTALL
    ):
        name, bar = label_ascii(marked(lab))
        reality[(name, bar)] = int(r)
    elements = {}
    for a, b in zip(marks, marks[1:] + [None]):
        seg = section[a.end() : b.start() if b else None]
        mats = re.findall(r"<pre><table>(.*?)</table></pre>", seg, flags=re.DOTALL)
        assert len(mats) >= 2, seg[:300]  # the 3 x 3 and the SU(2) matrix come first
        W, U = matrix_cells(mats[0]), matrix_cells(mats[1])
        sm = re.search(
            r'<td style=white-space:nowrap;"><center>(.*?)</center></td>',
            seg,
            flags=re.DOTALL,
        )
        assert sm, seg[:300]
        seitz = seitz_marked(sm.group(1))
        det = det3(W).real
        assert abs(abs(det) - 1) < 1e-9, det
        # (order in the double group, parity) from the matrices; the same from the symbol
        elements[seitz_ascii(seitz)] = (
            (math.lcm(matrix_order(W), matrix_order(U)), 1 if det > 0 else -1),
            seitz_invariants(seitz),
        )
    # classes: invariants from the matrices, checked against the Seitz symbols
    seen = []
    for c in classes:
        names = [seitz_ascii(e) for e in c["elements_raw"]]
        infos = [elements[n] for n in names]
        orders = {computed[0] for computed, _ in infos}
        parities = {computed[1] for computed, _ in infos}
        assert len(orders) == 1 and len(parities) == 1, (c["label"], names)
        for n, (computed, from_symbol) in zip(names, infos):
            assert from_symbol[:2] == computed, (n, computed, from_symbol)
        barred = {from_symbol[2] for _, from_symbol in infos}
        c["elements"] = names
        c["size"] = len(names)
        c["order"] = orders.pop()
        c["parity"] = parities.pop()
        c["barred"] = (
            "mixed" if len(barred) == 2 else ("all" if barred.pop() else "none")
        )
        seen += names
    assert sorted(seen) == sorted(elements), (
        "class lists do not partition the operations"
    )
    for ir in irreps:
        ir["reality"] = reality[(ir["cdml"], ir["barred"])]
    return {"classes": classes, "irreps": irreps}


# --------------------------------------------------------------------------
# TOML output
# --------------------------------------------------------------------------


def toml_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def toml_list(items):
    return "[" + ", ".join(toml_str(x) for x in items) + "]"


def convert_values(irreps, extra_need=1):
    M = extra_need
    for ir in irreps:
        for v in ir["values"]:
            M = math.lcm(M, conductor_need(v))
    F = cb.field(M)
    for ir in irreps:
        exact = [exact_value(F, v) for v in ir["values"]]
        for x, v in zip(exact, ir["values"]):
            approx = numeric_value(v)
            z = sum(
                float(a) * cmath.exp(2j * math.pi * k / M) for k, a in enumerate(x) if a
            )
            assert abs(z - approx) < 1e-9, (v, z, approx)
        ir["characters"] = [cb.to_grammar(F, x) for x in exact]
        ir["server_values"] = [server_ascii(v) for v in ir["values"]]
        assert all(cb.parse_value(F, s) == x for s, x in zip(ir["characters"], exact))
    return M


def write_point(group, url, data, fetched):
    hm = GROUPS[group][3]
    M = convert_values(data["irreps"])
    lines = [
        f"# Bilbao Crystallographic Server, Point Group Tables (POINT): character table of {group} ({hm})",
        "# Format: tests/reference/bilbao/characters/README.md",
        "# Written by tests/reference/bilbao/fetch.py from the server's HTML; the tests read this file only.",
        "",
        f"group = {toml_str(group)}",
        f"table = {toml_str(url)}",
        'status = "done"',
        f"source = {toml_str(SOURCE_POINT)}",
        f"fetched = {toml_str(fetched)}",
        "double = false",
        f"hermann_mauguin = {toml_str(hm)}",
        'notes = ""',
        f"conductor = {M}          # field of the `characters` values; zeta = exp(2 pi i / conductor)",
    ]
    if data["w"]:
        lines.append(
            f"server_w = {toml_str(data['w'])}    # the page's definition of w in `server_values`"
        )
    if data["notes"]:
        lines.append(f"server_notes = {toml_list(data['notes'])}")
    lines.append("")
    for c in data["classes"]:
        lines += [
            "[[classes]]",
            f"label = {toml_str(c['label'])}",
            f"size = {c['size']}",
            f"order = {c['order']}",
            f"parity = {c['parity']}",
            "",
        ]
    for ir in data["irreps"]:
        lines += [
            "[[irreps]]",
            f"label = {toml_str(ir['label'])}",
            f"bethe = {toml_str(ir['bethe'])}",
            "projective = false",
            f"server_values = {toml_list(ir['server_values'])}",
            f"characters = {toml_list(ir['characters'])}",
            "",
        ]
    return "\n".join(lines) + "\n"


def write_dpg(group, url, data, fetched):
    hm = GROUPS[group][3]
    M = convert_values(data["irreps"])
    dim = [int(ir["characters"][0]) for ir in data["irreps"]]
    d1 = next(k for k, c in enumerate(data["classes"]) if c["elements"] == ["d1"])
    lines = [
        f"# Bilbao Crystallographic Server, REPRESENTATIONS DPG: character table of the double group of {group} ({hm})",
        "# Format: tests/reference/bilbao/characters/README.md",
        "# Written by tests/reference/bilbao/fetch.py from the server's HTML; the tests read this file only.",
        "",
        f"group = {toml_str(group)}",
        f"table = {toml_str(url)}",
        'status = "done"',
        f"source = {toml_str(SOURCE_DPG)}",
        f"fetched = {toml_str(fetched)}",
        "double = true",
        f"hermann_mauguin = {toml_str(hm)}",
        'notes = ""',
        f"conductor = {M}          # field of the `characters` values; zeta = exp(2 pi i / conductor)",
        "",
    ]
    for c in data["classes"]:
        lines += [
            "[[classes]]",
            f"label = {toml_str(c['label'])}",
            f"elements = {toml_list(c['elements'])}",
            f"size = {c['size']}",
            f"order = {c['order']}",
            f"parity = {c['parity']}",
            f"barred = {toml_str(c['barred'])}",
            "",
        ]
    for ir, d in zip(data["irreps"], dim):
        F = cb.field(M)
        projective = cb.parse_value(F, ir["characters"][d1]) == F.rational(-d)
        assert projective == ir["barred"], (group, ir["label"])
        lines += [
            "[[irreps]]",
            f"label = {toml_str(ir['label'])}",
            f"bethe = {toml_str(ir['bethe'])}",
            f"cdml = {toml_str(ir['cdml'])}",
            f"projective = {'true' if projective else 'false'}",
            f"reality = {ir['reality']}",
            f"server_values = {toml_list(ir['server_values'])}",
            f"characters = {toml_list(ir['characters'])}",
            "",
        ]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main(argv):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("groups", nargs="*", help="Schoenflies symbols (default: all 32)")
    ap.add_argument("--raw-dir", type=Path, help="cache directory for the fetched HTML")
    ap.add_argument(
        "--turnstile-cookie",
        default=os.environ.get("BILBAO_TURNSTILE_COOKIE"),
        help="value of the server's turnstile_passed cookie (needed for the double-group pages)",
    )
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    args = ap.parse_args(argv)
    groups = [ALIASES.get(g, g) for g in args.groups] or list(GROUPS)
    fetched = dt.datetime.now(tz=dt.UTC).date().isoformat()
    args.out.mkdir(parents=True, exist_ok=True)
    for group in groups:
        num, sg_point, sg_dpg, _ = GROUPS[group]
        url = POINT_URL.format(sg=sg_point, num=num)
        raw = args.raw_dir / f"point_{sg_point}_{num}.html" if args.raw_dir else None
        data = parse_point(fetch(url, raw, None))
        (args.out / f"{group}.toml").write_text(write_point(group, url, data, fetched))
        url = DPG_URL.format(sg=sg_dpg, num=num)
        raw = args.raw_dir / f"dpg_{sg_dpg}_{num}.html" if args.raw_dir else None
        data = parse_dpg(fetch(url, raw, args.turnstile_cookie))
        (args.out / f"{group}-double.toml").write_text(
            write_dpg(group, url, data, fetched)
        )
        print(
            f"{group}: {len(data['classes'])} double-group classes, {len(data['irreps'])} irreps"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
