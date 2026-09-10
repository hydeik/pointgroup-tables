"""Prototype: canonical polyhedral bases (wayfinder ticket #16).

Throwaway code that constructs, verifies and freezes the canonical basis data
of the polyhedral families (docs/spec/basis-conventions.md, "Canonical basis:
polyhedral families"): for each of 2T, 2O, 2I, the images of the generator
pair q1, q2 in every irrep of R-dagger of dimension >= 2, adapted to the
principal element (D(c) diagonal, exponents descending), with the remaining
phases chosen so that every matrix lies in the smallest cyclotomic field the
irrep allows.

Construction path (all exact, in Q(zeta_M) with Fraction coefficients):

* E1/2: the spinor map psi(q) of the quaternion;
* the rotation irrep (T of T, T1 of O and I): the 3 x 3 rotation matrix of q;
* twists by a one-dimensional character (T2 = A2 x T1 and E5/2 = A2 x E1/2 of O)
  and Galois conjugates sqrt5 -> -sqrt5 (T2 and E7/2 of I);
* every other irrep: the isotypic block of a tensor product, cut out by the
  central projector of its character.

Every irrep is then adapted to c: the eigenvectors of D(s(c)) inside the block
are taken in descending exponent order, which fixes the basis up to one
scalar per vector; the scalars are found by solving |s|^2 = h_1 / h_i in the
real subfield of Q(zeta_M'), for the smallest conductor M' (a multiple of the
group's base conductor) in which a solution exists.

    uv run python prototypes/canonical_basis_prototype.py            # build, verify, report
    uv run python prototypes/canonical_basis_prototype.py --write    # also write the canonical TOML

Characters for the projectors come from the output-schema prototype (O, I) and
from the Koster scaffolding (T). Nothing here is production code.
"""

from __future__ import annotations

import itertools
import math
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
sys.path.insert(0, str(ROOT / "tools"))
import abstract_group_prototype as ag
import koster_bethe  # noqa: F401  (extends the schema prototype to T)
import output_schema_prototype as osp

OUT = (
    ROOT
    / "src"
    / "pointgroup_tables"
    / "data"
    / "conventions"
    / "basis"
    / "canonical.toml"
)

BASE_CONDUCTOR = {"T": 12, "O": 8, "I": 20}
PRINCIPAL = {"T": "C3^1(1)", "O": "C4^1(1)", "I": "C5^1(1)"}
ALPHA = {"T": None, "O": "sqrt2", "I": "phi"}

# --------------------------------------------------------------------------
# Field helpers on top of osp.Field (tuples of Fractions over the power basis)
# --------------------------------------------------------------------------


def field(M):
    return osp.Field(M)


def const(F, name):
    M = F.conductor
    if name == "sqrt2":
        return F.sqrt2()
    if name == "sqrt3":
        assert M % 12 == 0
        return F.add(F.zeta(M // 12), F.zeta(-(M // 12)))
    if name == "sqrt5":
        return F.sqrt5()
    if name == "phi":
        return F.phi()
    if name == "i":
        return F.i()
    raise KeyError(name)


def solve_q(A, b):
    """Solve A x = b over Fractions (A square, invertible)."""
    n = len(A)
    M = [list(A[r]) + [b[r]] for r in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        p = M[col][col]
        M[col] = [v / p for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * c for a, c in zip(M[r], M[col])]
    return [M[r][n] for r in range(n)]


def finv(F, x):
    n = F.degree
    cols = [F.mul(x, F.zeta(j)) for j in range(n)]
    A = [[cols[j][r] for j in range(n)] for r in range(n)]
    e0 = [Fr(1)] + [Fr(0)] * (n - 1)
    return tuple(solve_q(A, e0))


def lift(F, F2, x):
    """Embed x in Q(zeta_M) into Q(zeta_M'), M | M'."""
    step = F2.conductor // F.conductor
    assert F2.conductor % F.conductor == 0
    out = F2.zero()
    for j, a in enumerate(x):
        if a:
            out = F2.add(out, F2.scale(F2.zeta(j * step), a))
    return out


def real_basis(F):
    """A Q-basis of the real subfield: 1 and zeta^k + zeta^-k."""
    out = [F.rational(1)]
    for k in range(1, F.degree // 2):
        out.append(F.add(F.zeta(k), F.zeta(-k)))
    return out


def coords(F, x, basis):
    """Rational coordinates of x over a Q-linearly independent list, or None."""
    n = F.degree
    m = len(basis)
    # least squares free: solve the n x m system exactly by elimination on an augmented matrix
    rows = [[basis[j][r] for j in range(m)] + [x[r]] for r in range(n)]
    pivots = []
    r0 = 0
    for col in range(m):
        piv = next((r for r in range(r0, n) if rows[r][col] != 0), None)
        if piv is None:
            continue
        rows[r0], rows[piv] = rows[piv], rows[r0]
        p = rows[r0][col]
        rows[r0] = [v / p for v in rows[r0]]
        for r in range(n):
            if r != r0 and rows[r][col] != 0:
                f = rows[r][col]
                rows[r] = [a - f * c for a, c in zip(rows[r], rows[r0])]
        pivots.append(col)
        r0 += 1
    for r in range(r0, n):
        if rows[r][m] != 0:
            return None
    sol = [Fr(0)] * m
    for k, col in enumerate(pivots):
        sol[col] = rows[k][m]
    return sol


# --------------------------------------------------------------------------
# Matrices over F
# --------------------------------------------------------------------------


def mmul(F, A, B):
    return osp.mat_mul(F, A, B)


def mdag(F, A):
    return [[F.conj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def meq(A, B):
    return all(a == b for ra, rb in zip(A, B) for a, b in zip(ra, rb))


def midentity(F, d):
    return osp.mat_identity(F, d)


def mkron(F, A, B):
    da, db = len(A), len(B)
    out = []
    for i in range(da):
        for k in range(db):
            row = []
            for j in range(da):
                for l in range(db):
                    row.append(F.mul(A[i][j], B[k][l]))
            out.append(row)
    return out


def mscale(F, A, s):
    return [[F.mul(a, s) for a in row] for row in A]


def madd(F, A, B):
    return [[F.add(a, b) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mtrace(F, A):
    return osp.mat_trace(F, A)


def is_unitary(F, A):
    return meq(mmul(F, mdag(F, A), A), midentity(F, len(A)))


def nullspace(F, A):
    """Basis of {v : A v = 0} over F, by elimination with field inverses."""
    m, n = len(A), len(A[0])
    rows = [list(r) for r in A]
    pivots = []
    r0 = 0
    for col in range(n):
        piv = next((r for r in range(r0, m) if not F.is_zero(rows[r][col])), None)
        if piv is None:
            continue
        rows[r0], rows[piv] = rows[piv], rows[r0]
        inv = finv(F, rows[r0][col])
        rows[r0] = [F.mul(v, inv) for v in rows[r0]]
        for r in range(m):
            if r != r0 and not F.is_zero(rows[r][col]):
                f = rows[r][col]
                rows[r] = [F.sub(a, F.mul(f, c)) for a, c in zip(rows[r], rows[r0])]
        pivots.append(col)
        r0 += 1
    free = [c for c in range(n) if c not in pivots]
    basis = []
    for fcol in free:
        v = [F.zero()] * n
        v[fcol] = F.rational(1)
        for k, pcol in enumerate(pivots):
            v[pcol] = F.neg(rows[k][fcol])
        basis.append(v)
    return basis


def rank(F, vectors):
    if not vectors:
        return 0
    A = [list(v) for v in vectors]
    return len(A) - len(
        nullspace(F, [[A[r][c] for r in range(len(A))] for c in range(len(A[0]))])
    )


# --------------------------------------------------------------------------
# Groups, quaternions and the reference constructions
# --------------------------------------------------------------------------


class Poly:
    """One polyhedral rotation group R with its binary cover R-dagger."""

    def __init__(self, kind, M):
        self.kind = kind
        self.ctx = osp.Context(kind)  # proper group: elements are (perm, +1)
        self.G = self.ctx.G
        self.cov = self.ctx.D.cover
        self.Q = list(self.cov.elements)  # quaternions of R-dagger
        self.index = {q: k for k, q in enumerate(self.Q)}
        self.F = field(M)
        self.alpha = const(self.F, ALPHA[kind]) if ALPHA[kind] else None
        self.labels = self.ctx.labels
        self.r_index = {g[0]: k for k, g in enumerate(self.G.elements)}
        c = next(g for g in self.G.elements if self.labels[g] == PRINCIPAL[kind])
        self.c = self.cov.section(c[0])
        self.N = ag.element_order(self.G, c)
        self.q1, self.q2 = self.generators()

    def generators(self):
        # the spec's pairs, in the prototype's coordinate encoding (a + b alpha)
        spec = {
            "T": (((Fr(1, 2), 0),) * 4, ((0, 0), (0, 0), (0, 0), (Fr(1), 0))),
            "O": (
                ((0, Fr(1, 2)), (0, 0), (0, 0), (0, Fr(1, 2))),
                ((0, 0), (0, 0), (0, Fr(1, 2)), (0, Fr(1, 2))),
            ),
            "I": (
                ((0, Fr(1, 2)), (Fr(-1, 2), Fr(1, 2)), (0, 0), (Fr(1, 2), 0)),
                ((Fr(1, 2), 0), (Fr(-1, 2), Fr(1, 2)), (0, Fr(1, 2)), (0, 0)),
            ),
        }[self.kind]
        out = []
        for q in spec:
            q = tuple((Fr(a), Fr(b)) for a, b in q)
            assert q in self.index, ("generator not in the cover", self.kind, q)
            out.append(q)
        return out

    def coord(self, ab):
        a, b = ab
        x = self.F.rational(a)
        if b and self.alpha is not None:
            x = self.F.add(x, self.F.scale(self.alpha, Fr(b)))
        elif b:
            x = self.F.add(x, self.F.rational(b))  # 2T: alpha = 1
        return x

    def wxyz(self, q):
        return [self.coord(ab) for ab in q]

    def spinor(self, q):
        F = self.F
        w, x, y, z = self.wxyz(q)
        i = F.i()
        return [
            [F.add(w, F.mul(i, z)), F.add(F.neg(y), F.mul(i, x))],
            [F.add(y, F.mul(i, x)), F.sub(w, F.mul(i, z))],
        ]

    def rotation(self, q):
        F = self.F
        w, x, y, z = self.wxyz(q)
        two = Fr(2)

        def sq(a):
            return F.mul(a, a)

        def t(a, b, c, d, sign):  # 2 (a b + sign c d)
            p = F.mul(a, b)
            r = F.mul(c, d)
            return F.scale(F.add(p, r) if sign > 0 else F.sub(p, r), two)

        one = F.rational(1)
        return [
            [
                F.sub(one, F.scale(F.add(sq(y), sq(z)), two)),
                t(x, y, w, z, -1),
                t(x, z, w, y, +1),
            ],
            [
                t(x, y, w, z, +1),
                F.sub(one, F.scale(F.add(sq(x), sq(z)), two)),
                t(y, z, w, x, -1),
            ],
            [
                t(x, z, w, y, -1),
                t(y, z, w, x, +1),
                F.sub(one, F.scale(F.add(sq(x), sq(y)), two)),
            ],
        ]

    def character(self, symbol):
        """chi as a function of the quaternion, from the schema prototype's irreps."""
        ir = next(ir for ir in self.ctx.irreps if ir.symbol == symbol)
        F = self.F

        def chi(q):
            r = self.cov.proj(q)
            x = ir.chars[self.r_index[r]]
            if ir.projective and self.cov.section(r) != q:
                x = F.neg(x)
            return x

        return chi, ir


class Rep:
    """A representation of R-dagger as a dict quaternion -> matrix, plus its symbol."""

    def __init__(self, P, symbol, mats):
        self.P = P
        self.symbol = symbol
        self.mats = mats
        self.dim = len(next(iter(mats.values())))

    @staticmethod
    def from_function(P, symbol, fn):
        return Rep(P, symbol, {q: fn(q) for q in P.Q})

    def tensor(self, other):
        F = self.P.F
        return Rep(
            self.P,
            f"{self.symbol}x{other.symbol}",
            {q: mkron(F, self.mats[q], other.mats[q]) for q in self.P.Q},
        )

    def twist(self, symbol, chi):
        F = self.P.F
        return Rep(
            self.P, symbol, {q: mscale(F, self.mats[q], chi(q)) for q in self.P.Q}
        )

    def galois(self, symbol, k):
        """Apply zeta -> zeta^k entrywise (k coprime to M)."""
        F = self.P.F
        M = F.conductor
        assert math.gcd(k, M) == 1

        def sigma(x):
            out = F.zero()
            for j, a in enumerate(x):
                if a:
                    out = F.add(out, F.scale(F.zeta(j * k % M), a))
            return out

        return Rep(
            self.P,
            symbol,
            {q: [[sigma(a) for a in row] for row in A] for q, A in self.mats.items()},
        )


# --------------------------------------------------------------------------
# Adapting an irreducible block to the principal element
# --------------------------------------------------------------------------


def projector(rep, chi, d):
    """(d / |R-dagger|) sum_q conj(chi(q)) D(q)."""
    F = rep.P.F
    n = rep.dim
    acc = [[F.zero()] * n for _ in range(n)]
    for q, A in rep.mats.items():
        c = F.conj(chi(q))
        if F.is_zero(c):
            continue
        acc = madd(F, acc, mscale(F, A, c))
    return mscale(F, acc, F.rational(Fr(d, len(rep.P.Q))))


def adapt(rep, symbol, chi, d, projective):
    """Cut the isotypic block of chi out of rep and express it in the eigenbasis of D(s(c)).

    Returns (V columns as list of vectors, norms h_i, exponents e_i) with the
    eigenvalue of column i equal to zeta_{2N}^{e_i} (projective) or zeta_N^{e_i}
    (ordinary), exponents descending.
    """
    P = rep.P
    F = P.F
    n = rep.dim
    Pr = projector(rep, chi, d)
    Dc = rep.mats[P.c]
    N = P.N
    order = 2 * N if projective else N
    M = F.conductor
    assert M % order == 0, (M, order)
    exps = range(
        N if projective else N // 2, -N if projective else -((N - 1) // 2) - 1, -1
    )
    columns, norms, found = [], [], []
    for e in exps:
        lam = F.zeta((M // order) * e)
        A = [
            [F.sub(Dc[i][j], lam if i == j else F.zero()) for j in range(n)]
            for i in range(n)
        ]
        ns = nullspace(F, A)
        if not ns:
            continue
        imgs = [
            [osp._sum(F, (F.mul(Pr[i][j], v[j]) for j in range(n))) for i in range(n)]
            for v in ns
        ]
        imgs = [v for v in imgs if any(not F.is_zero(a) for a in v)]
        if not imgs:
            continue
        assert rank(F, imgs) == 1, (symbol, e, "multiplicity > 1")
        v = imgs[0]
        # clear denominators for readability
        den = 1
        for a in v:
            for c in a:
                den = math.lcm(den, c.denominator)
        v = [F.scale(a, Fr(den)) for a in v]
        h = osp._sum(F, (F.mul(F.conj(a), a) for a in v))
        columns.append(v)
        norms.append(h)
        found.append(e)
    assert len(columns) == d, (symbol, len(columns), d, found)
    return columns, norms, found


def restrict(rep, columns, norms):
    """D_irrep(q) = diag(1/h_i) V^dagger D(q) V, exact in F."""
    F = rep.P.F
    n = rep.dim
    d = len(columns)
    hinv = [finv(F, h) for h in norms]
    out = {}
    for q, A in rep.mats.items():
        AV = [
            [
                osp._sum(F, (F.mul(A[i][j], columns[k][j]) for j in range(n)))
                for k in range(d)
            ]
            for i in range(n)
        ]
        B = [
            [
                F.mul(
                    hinv[i],
                    osp._sum(
                        F, (F.mul(F.conj(columns[i][m]), AV[m][k]) for m in range(n))
                    ),
                )
                for k in range(d)
            ]
            for i in range(d)
        ]
        out[q] = B
    return out


def norm_solutions(
    F, targets, coeff_range=(0, 1, -1, 2, -2), dens=(1, 2, 3, 4, 6, 8, 12)
):
    """For each target r (real, in F) find s = x + i y with x, y in the real subfield and |s|^2 = r.

    Returns a list of s (or None where no small solution exists)."""
    basis = real_basis(F)
    if len(basis) > 4:  # degree-16 fields and beyond: keep the search affordable
        coeff_range, dens = (0, 1, -1), (1, 2, 3, 4, 6)
    i_ = F.i()
    squares = {}
    cands = []
    for den in dens:
        for cs in itertools.product(coeff_range, repeat=len(basis)):
            x = F.zero()
            for c, b in zip(cs, basis):
                if c:
                    x = F.add(x, F.scale(b, Fr(c, den)))
            cands.append(x)
    for x in cands:
        squares.setdefault(F.mul(x, x), x)
    out = []
    for r in targets:
        best = None
        for x2, x in squares.items():
            y2 = F.sub(r, x2)
            if y2 in squares:
                y = squares[y2]
                s = F.add(x, F.mul(i_, y))
                key = (
                    sum(abs(c) for c in x) + sum(abs(c) for c in y),  # ty: ignore[not-iterable]
                    not F.is_zero(y),
                    F.complex(x).real < 0,
                    F.complex(y).real < 0,
                    s,
                )
                if best is None or key < best[0]:
                    best = (key, s)
        out.append(best[1] if best else None)
    return out


def freeze(P, rep, symbol, projective, chi, d, conductors):
    """Adapt, then find the smallest conductor with a unitary realisation; return (M', matrices dict, scalars)."""
    t0 = time.process_time()
    columns, norms, exps = adapt(rep, symbol, chi, d, projective)
    t1 = time.process_time()
    mats = restrict(rep, columns, norms)
    t2 = time.process_time()
    print(f"   {symbol}: adapt {t1 - t0:.1f} s, restrict {t2 - t1:.1f} s", flush=True)
    F = P.F
    targets = [F.mul(norms[0], finv(F, h)) for h in norms]  # h_1 / h_i, real
    for M2 in conductors:
        F2 = field(M2)
        t2 = [lift(F, F2, t) for t in targets]
        print(
            f"   {symbol}: targets {[to_grammar(F2, t) for t in t2]} in Q(zeta_{M2})",
            flush=True,
        )
        t3 = time.process_time()
        found = norm_solutions(F2, t2)
        print(f"   {symbol}: norm search {time.process_time() - t3:.1f} s", flush=True)
        sols = [s for s in found if s is not None]
        if len(sols) == len(found):
            sinv = [finv(F2, s) for s in sols]
            out = {}
            for q, A in mats.items():
                out[q] = [
                    [
                        F2.mul(F2.mul(sinv[i], lift(F, F2, A[i][j])), sols[j])
                        for j in range(d)
                    ]
                    for i in range(d)
                ]
            return M2, F2, out, sols, exps
    raise AssertionError(f"{symbol}: no normalisation found in conductors {conductors}")


# --------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------


def verify(P, F2, symbol, mats, chi, projective, exps):
    """Closure from the generators, unitarity, traces, diagonal c, irreducibility."""
    cov = P.cov
    gens = [P.q1, P.q2]
    e = P.Q[0]
    assert all(F2.is_zero(x) for x in P.wxyz(e)[1:]) and P.wxyz(e)[0] == P.F.rational(
        1
    ), "Q[0] is not the identity"
    d = len(mats[e])
    reached = {e: midentity(F2, d)}
    frontier = [e]
    while frontier:
        nxt = []
        for q in frontier:
            for g in gens:
                p = cov.mul(q, g)
                B = mmul(F2, reached[q], mats[g])
                if p in reached:
                    assert meq(reached[p], B), (symbol, "closure inconsistent")
                else:
                    reached[p] = B
                    nxt.append(p)
        frontier = nxt
    assert len(reached) == len(P.Q), (symbol, "generators do not generate")
    for q in P.Q:
        assert meq(reached[q], mats[q]), (
            symbol,
            "closure differs from the constructed matrix",
        )
    for g in gens:
        assert is_unitary(F2, mats[g]), (symbol, "generator not unitary")
    total = F2.zero()
    for q in P.Q:
        tr = mtrace(F2, mats[q])
        assert tr == lift(P.F, F2, chi(q)), (
            symbol,
            "trace differs from the character",
            q,
        )
        total = F2.add(total, F2.mul(F2.conj(tr), tr))
    assert total == F2.rational(len(P.Q)), (symbol, "not irreducible")
    Dc = mats[P.c]
    for i in range(d):
        for j in range(d):
            if i != j:
                assert F2.is_zero(Dc[i][j]), (symbol, "D(c) not diagonal")
    order = 2 * P.N if projective else P.N
    for i, e_ in enumerate(exps):
        assert Dc[i][i] == F2.zeta((F2.conductor // order) * e_), (
            symbol,
            "eigenvalue order",
        )
    ebar = cov.ebar
    assert meq(
        mats[ebar], mscale(F2, midentity(F2, d), F2.rational(-1 if projective else 1))
    ), (symbol, "full turn")


# --------------------------------------------------------------------------
# Value grammar: printing and parsing
# --------------------------------------------------------------------------


def fmt_q(c):
    return str(c.numerator) if c.denominator == 1 else f"{c.numerator}/{c.denominator}"


def nice_bases(F):
    M = F.conductor
    names = ["1"]
    if M % 4 == 0:
        names.append("i")
    for root, req in (("sqrt2", 8), ("sqrt3", 12), ("phi", 5)):
        if M % req == 0:
            names.append(root)
            if M % 4 == 0:
                names.append(f"i*{root}")
    return names


def value_of(F, name):
    if name == "1":
        return F.rational(1)
    if name.startswith("i*"):
        return F.mul(F.i(), const(F, name[2:]))
    return const(F, name)


def to_grammar(F, x):
    """Shortest of: rational, r*zeta^k, combination over the nice basis, power-basis polynomial."""
    if all(a == 0 for a in x[1:]):
        return fmt_q(x[0])
    M = F.conductor
    for r in (1, -1, 2, -2, Fr(1, 2), Fr(-1, 2)):
        for k in range(1, M):
            z = F.zeta(k)
            if x == F.scale(z, Fr(r)):
                head = "" if r == 1 else "-" if r == -1 else fmt_q(Fr(r)) + "*"
                return f"{head}zeta^{k}" if k != 1 else f"{head}zeta"
    names = nice_bases(F)
    vecs = [value_of(F, n) for n in names]
    cs = coords(F, x, vecs)
    if cs is not None:
        return _terms(zip(cs, names))
    return _terms(
        (a, f"zeta^{j}" if j > 1 else "zeta" if j == 1 else "1")
        for j, a in enumerate(x)
    )


def _terms(pairs):
    out = ""
    for c, name in pairs:
        if c == 0:
            continue
        sign = "-" if c < 0 else "+"
        mag = abs(c)
        if name == "1":
            body = fmt_q(mag)
        elif mag == 1:
            body = name
        else:
            body = f"{fmt_q(mag)}*{name}"
        out += (
            (f"{sign}{body}" if not out else f" {sign} {body}")
            if sign == "-" or out
            else body
        )
    return out or "0"


class Parser:
    """Recursive descent for the spec's value grammar, into F."""

    def __init__(self, F, text):
        self.F = F
        self.toks = self.tokenize(text)
        self.pos = 0

    @staticmethod
    def tokenize(text):
        toks, i = [], 0
        while i < len(text):
            ch = text[i]
            if ch.isspace():
                i += 1
            elif ch.isdigit():
                j = i
                while j < len(text) and text[j].isdigit():
                    j += 1
                toks.append(("int", int(text[i:j])))
                i = j
            elif ch.isalpha():
                j = i
                while j < len(text) and (text[j].isalnum()):
                    j += 1
                toks.append(("name", text[i:j]))
                i = j
            else:
                toks.append(("op", ch))
                i += 1
        return toks

    def peek(self):
        return self.toks[self.pos] if self.pos < len(self.toks) else (None, None)

    def take(self):
        t = self.toks[self.pos]
        self.pos += 1
        return t

    def expr(self):
        x = self.term()
        while self.peek() == ("op", "+") or self.peek() == ("op", "-"):
            op = self.take()[1]
            y = self.term()
            x = self.F.add(x, y) if op == "+" else self.F.sub(x, y)
        return x

    def term(self):
        x = self.factor()
        while self.peek() == ("op", "*") or self.peek() == ("op", "/"):
            op = self.take()[1]
            y = self.factor()
            x = self.F.mul(x, y) if op == "*" else self.F.mul(x, finv(self.F, y))
        return x

    def factor(self):
        neg = False
        if self.peek() == ("op", "-"):
            self.take()
            neg = True
        x = self.atom()
        if self.peek() == ("op", "^"):
            self.take()
            k = self.take()
            assert k[0] == "int"
            x = self.F.power(x, k[1])
        return self.F.neg(x) if neg else x

    def atom(self):
        kind, val = self.take()
        if kind == "int":
            return self.F.rational(val)
        if kind == "name":
            if val == "zeta":
                return self.F.zeta(1)
            return const(self.F, val)
        assert (kind, val) == ("op", "("), (kind, val)
        x = self.expr()
        assert self.take() == ("op", ")")
        return x


def parse_value(F, text):
    p = Parser(F, text)
    x = p.expr()
    assert p.pos == len(p.toks), text
    return x


# --------------------------------------------------------------------------
# The plan per group
# --------------------------------------------------------------------------


def build_group(kind, conductors):
    P = Poly(kind, BASE_CONDUCTOR[kind])
    F = P.F
    results = {}  # symbol -> (M', F2, mats, projective, scalars, exps, path)

    def chi_of(symbol):
        chi, ir = P.character(symbol)
        return chi, ir.dimension, ir.projective

    def do(symbol, rep, path):
        chi, d, projective = chi_of(symbol)
        t0 = time.process_time()
        M2, F2, mats, sols, exps = freeze(
            P, rep, symbol, projective, chi, d, conductors
        )
        t1 = time.process_time()
        verify(P, F2, symbol, mats, chi, projective, exps)
        print(
            f"   built {symbol} in {t1 - t0:.1f} s, verified in {time.process_time() - t1:.1f} s",
            flush=True,
        )
        results[symbol] = (M2, F2, mats, projective, sols, exps, path)
        return (
            Rep(
                P,
                symbol,
                {
                    q: [
                        [
                            lift(F2, F, a) if F2.conductor == F.conductor else a
                            for a in row
                        ]
                        for row in A
                    ]
                    for q, A in mats.items()
                },
            )
            if M2 == F.conductor
            else None
        )

    def do_product(symbol, left, right, path):
        """The plain tensor product of two frozen irreps, columns ordered by weight.

        Used where D(s(c)) has a repeated eigenvalue (I5/2 of 2I: -1 twice), so
        rules 1-2 alone do not fix the basis; the product basis of the canonical
        factors does, with the SU(2) weight 2m = e_left + e_right (unreduced)
        as the exponent."""
        chi, d, projective = chi_of(symbol)
        Ml, Fl, ml, pl, _, el, _ = results[left]
        Mr, _, mr, pr, _, er, _ = results[right]
        assert Ml == Mr
        F2 = Fl
        scale_l = 1 if pl else 2  # exponents of an ordinary factor are in zeta_N units
        scale_r = 1 if pr else 2
        weights = [
            (a * scale_l + b * scale_r, ia * len(er) + ib)
            for ia, a in enumerate(el)
            for ib, b in enumerate(er)
        ]
        weights.sort(key=lambda t: -t[0])
        perm = [col for _, col in weights]
        exps = [w for w, _ in weights]
        assert len(set(exps)) == d, (symbol, exps)
        mats = {}
        for q in P.Q:
            K = mkron(F2, ml[q], mr[q])
            mats[q] = [[K[perm[i]][perm[j]] for j in range(d)] for i in range(d)]
        t1 = time.process_time()
        verify(P, F2, symbol, mats, chi, projective, exps)
        print(
            f"   built {symbol} as a product, verified in {time.process_time() - t1:.1f} s",
            flush=True,
        )
        results[symbol] = (Ml, F2, mats, projective, [F2.rational(1)] * d, exps, path)

    spin = Rep.from_function(P, "E1/2", P.spinor)
    if kind == "T":
        rot = Rep.from_function(P, "T", P.rotation)
        do("E1/2", spin, "spinor map psi(q)")
        do("T", rot, "rotation matrix of q, in the eigenbasis of C3")
        prod = spin.tensor(rot)
        do("1E3/2", prod, "block of E1/2 x T")
        do("2E3/2", prod, "block of E1/2 x T")
    elif kind == "O":
        rot = Rep.from_function(P, "T1", P.rotation)
        a2, _, _ = chi_of("A2")
        do("E1/2", spin, "spinor map psi(q)")
        do("T1", rot, "rotation matrix of q, in the eigenbasis of C4")
        do("T2", rot.twist("T2", a2), "A2 x T1")
        do("E5/2", spin.twist("E5/2", a2), "A2 x E1/2")
        tt = rot.tensor(rot)
        e_rep = do("E", tt, "block of T1 x T1")
        do(
            "G3/2",
            spin.tensor(e_rep) if e_rep else spin.tensor(rot),
            "E1/2 x E" if e_rep else "block of E1/2 x T1",
        )
    else:
        rot = Rep.from_function(P, "T1", P.rotation)
        sigma = 13  # zeta20 -> zeta20^13: sqrt5 -> -sqrt5, i fixed
        do("E1/2", spin, "spinor map psi(q)")
        do(
            "E7/2",
            spin.galois("E7/2", sigma),
            "Galois conjugate sqrt5 -> -sqrt5 of E1/2",
        )
        do("T1", rot, "rotation matrix of q, in the eigenbasis of C5")
        rot2 = rot.galois("T2", sigma)
        do("T2", rot2, "Galois conjugate sqrt5 -> -sqrt5 of T1")
        tt = rot.tensor(rot)
        do("H", tt, "block of T1 x T1")
        do("G", rot.tensor(rot2), "block of T1 x T2")
        do("G3/2", spin.tensor(rot), "block of E1/2 x T1")
        do_product(
            "I5/2",
            "E1/2",
            "T2",
            "E1/2 x T2 in the product basis, columns by weight (sum of factor exponents)",
        )
    return P, results


def irrep_order(P):
    return [ir.symbol for ir in P.ctx.irreps if ir.dimension >= 2]


def report(kind, P, results):
    Mmax = max(r[0] for r in results.values())
    print(
        f"== 2{kind}: base conductor {P.F.conductor}, group conductor {Mmax}, {len(P.Q)} elements, N(c) = {P.N}"
    )
    print(
        f"   generators: {P.labels[(P.cov.proj(P.q1), 1)]}, {P.labels[(P.cov.proj(P.q2), 1)]}"
        f"  (lifted: {P.cov.section(P.cov.proj(P.q1)) == P.q1}, {P.cov.section(P.cov.proj(P.q2)) == P.q2})"
    )
    for symbol in irrep_order(P):
        M2, F2, mats, _, sols, exps, path = results[symbol]
        scal = ", ".join(to_grammar(F2, s) for s in sols)
        print(
            f"   {symbol:6} dim {len(mats[P.Q[0]])}  conductor {M2:3}  exponents {list(exps)}  scalars [{scal}]  <- {path}"
        )
    return Mmax


def write_toml(all_results):
    lines = [
        "# Canonical basis convention: the polyhedral generator data frozen by",
        "# prototypes/canonical_basis_prototype.py (wayfinder ticket #16).",
        "# Regenerate with: uv run python prototypes/canonical_basis_prototype.py --write",
        "",
        'name = "canonical"',
        "",
        "[axial]",
        'rule = "canonical"',
        "",
    ]
    for kind in all_results:
        P, results = all_results[kind]
        Mmax = max(r[0] for r in results.values())
        Fg = field(Mmax)
        g1 = P.labels[(P.cov.proj(P.q1), 1)]
        g2 = P.labels[(P.cov.proj(P.q2), 1)]
        lines.append(f"[polyhedral.2{kind}]")
        lines.append(f"conductor = {Mmax}")
        lines.append(f'generators = ["{g1}", "{g2}"]')
        lines.append("")
        for symbol in irrep_order(P):
            _, F2, mats, _, _, _, path = results[symbol]
            key = symbol if "/" not in symbol else f'"{symbol}"'
            lines.append(f"[polyhedral.2{kind}.irreps.{key}]  # {path}")
            lines.append("matrices = [")
            for g in (P.q1, P.q2):
                A = mats[g]
                rows = []
                for row in A:
                    rows.append(
                        "["
                        + ", ".join(
                            '"' + to_grammar(Fg, lift(F2, Fg, a)) + '"' for a in row
                        )
                        + "]"
                    )
                lines.append("  [" + ", ".join(rows) + "],")
            lines.append("]")
            lines.append("")
    return "\n".join(lines)


def roundtrip(path, all_results):
    import tomllib

    data = tomllib.loads(path.read_text())
    for kind in all_results:
        P, results = all_results[kind]
        sec = data["polyhedral"][f"2{kind}"]
        Fg = field(sec["conductor"])
        for symbol, entry in sec["irreps"].items():
            _, F2, mats, _, _, _, _ = results[symbol]
            for g, A in zip((P.q1, P.q2), entry["matrices"]):
                parsed = [[parse_value(Fg, s) for s in row] for row in A]
                assert meq(
                    parsed, [[lift(F2, Fg, a) for a in row] for row in mats[g]]
                ), (kind, symbol, "round trip")
    print(f"round trip of {path} ok")


def main(argv):
    write = "--write" in argv
    t0 = time.process_time()
    all_results = {}
    kinds = [a[7:] for a in argv if a.startswith("--only=")] or ["T", "O", "I"]
    for kind in kinds:
        base = BASE_CONDUCTOR[kind]
        conductors = [base * k for k in (1, 2, 3, 4, 6)]
        P, results = build_group(kind, conductors)
        report(kind, P, results)
        all_results[kind] = (P, results)
    print(f"built and verified in {time.process_time() - t0:.1f} s")
    text = write_toml(all_results)
    if write:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(text + "\n")
        print(f"wrote {OUT}")
        roundtrip(OUT, all_results)
    else:
        tmp = Path(__file__).with_name("canonical_basis_prototype.out.toml")
        tmp.write_text(text + "\n")
        roundtrip(tmp, all_results)
        tmp.unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
