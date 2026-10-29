"""Prototype: Clebsch-Gordan coefficients by projection (wayfinder ticket #50).

Throwaway code for docs/spec/clebsch-gordan.md on top of the library: the
coupling matrix of every pair ``alpha <= beta`` of a group in the canonical
basis (the trial-state projection, the symmetric-first multiplicity rule, the
norm-factor phase rule), its unitary transform into a chain-adapted basis, the
checks ``CG-1`` ... ``CG-5`` and the pair records of the ``cg`` document.

Every projector sum runs in the version-1 field ``Q(zeta_{M_G})`` of the
canonical matrices (the squared norms are rational there); only the
normalisation runs in the ``cg`` field ``Q(zeta_{M_cg})``, and the matrices are
lifted at the output boundary. Nothing here is production code.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from chain_prototype_basis import ChainBasis, components_of, generating_set
from chain_prototype_exact import NormFactorTable, norm_factor

from pointgroup_tables.document import Model
from pointgroup_tables.field import CyclotomicField, Matrix, Value, cyclotomic_field
from pointgroup_tables.irreps import Irrep

Vector = list[Value]


def cg_conductor(model: Model, M_basis: int | None = None) -> int:
    """``lcm(M_basis, 8)`` axial, ``lcm(M_basis, 24)`` polyhedral."""
    M = model.conductor if M_basis is None else M_basis
    return math.lcm(M, 24 if model.spec.polyhedral else 8)


@dataclass
class Copy:
    gamma: int  # global irrep index of the constituent
    r: int  # copy index
    exchange: int | None  # +1 / -1 for alpha == beta, else None
    trial: int  # the trial product state j0
    N: Value  # the squared norm, in Q(zeta_{M_G})
    columns: list[Vector]  # d_gamma columns in product order, in Q(zeta_{M_cg})


@dataclass
class Pair:
    alpha: int
    beta: int
    copies: list[Copy]
    multiplicities: dict[int, tuple[int, int, int]]  # gamma -> (n, n_plus, n_minus)

    @property
    def dimension(self) -> int:
        return len(self.copies[0].columns[0]) if self.copies else 0

    def matrix(self) -> Matrix:
        """The coupling matrix ``U`` (rows product states, columns coupled vectors)."""
        columns = [c for copy in self.copies for c in copy.columns]
        return [[col[i] for col in columns] for i in range(len(columns[0]))]

    def column_records(self, irreps: list[Irrep]) -> list[dict[str, Any]]:
        out = []
        index = 0
        for copy in self.copies:
            for k in range(len(copy.columns)):
                record: dict[str, Any] = {
                    "index": index,
                    "irrep": copy.gamma,
                    "copy": copy.r,
                    "component": k,
                }
                if copy.exchange is not None:
                    record["exchange"] = copy.exchange
                out.append(record)
                index += 1
        return out

    def constituent_records(self, irreps: list[Irrep]) -> list[dict[str, Any]]:
        out = []
        first = 0
        for gamma in sorted(self.multiplicities):
            n = self.multiplicities[gamma][0]
            count = n * irreps[gamma].dimension
            out.append({"irrep": gamma, "multiplicity": n, "columns": [first, count]})
            first += count
        return out


# ---- the algorithm -----------------------------------------------------------------------


def multiplicities(
    model: Model, alpha: Irrep, beta: Irrep, gamma: Irrep
) -> tuple[int, int, int]:
    """``(n, n_plus, n_minus)``; the split is meaningful for ``alpha is beta`` only."""
    F = model.field
    order = model.order
    total = F.zero
    sym = F.zero
    same = alpha is beta
    for g in range(order):
        chi = F.mul(alpha.characters[g], beta.characters[g])
        weight = F.conjugate(gamma.characters[g])
        total = F.add(total, F.mul(chi, weight))
        if same:
            D = alpha.matrices[g]
            d = alpha.dimension
            trace = F.total(F.mul(D[i][k], D[k][i]) for i in range(d) for k in range(d))
            sym = F.add(sym, F.mul(F.scale(F.add(chi, trace), Fraction(1, 2)), weight))
    n = F.scale(total, Fraction(1, order))
    if not F.is_rational(n) or n[0].denominator != 1:
        raise ArithmeticError("non-integral multiplicity")
    count = int(n[0])
    if not same:
        return count, count, 0
    n_plus = F.scale(sym, Fraction(1, order))
    if not F.is_rational(n_plus) or n_plus[0].denominator != 1:
        raise ArithmeticError("non-integral symmetric multiplicity")
    plus = int(n_plus[0])
    return count, plus, count - plus


def _product_column(
    F: CyclotomicField,
    alpha: Irrep,
    beta: Irrep,
    weights: list[Value],
    j: int,
    d_gamma: int,
    order: int,
) -> Vector:
    """``(d_gamma / |G|) sum_g weights[g] (D_alpha(g) x D_beta(g)) e_j``."""
    da, db = alpha.dimension, beta.dimension
    j1, j2 = divmod(j, db)
    out = [F.zero] * (da * db)
    for g in range(order):
        w = weights[g]
        if F.is_zero(w):
            continue
        A, B = alpha.matrices[g], beta.matrices[g]
        for i1 in range(da):
            a = A[i1][j1]
            if F.is_zero(a):
                continue
            wa = F.mul(w, a)
            for i2 in range(db):
                b = B[i2][j2]
                if not F.is_zero(b):
                    idx = i1 * db + i2
                    out[idx] = F.add(out[idx], F.mul(wa, b))
    scale = Fraction(d_gamma, order)
    return [F.scale(x, scale) for x in out]


def _product_apply(
    F: CyclotomicField, alpha: Irrep, beta: Irrep, g: int, w: Vector
) -> Vector:
    """``(D_alpha(g) x D_beta(g)) w`` as ``D_alpha W D_beta^T``."""
    da, db = alpha.dimension, beta.dimension
    A, B = alpha.matrices[g], beta.matrices[g]
    # X = W B^T : X[i1][j] = sum_j' W[i1][j'] B[j][j']
    X = [
        [
            F.total(
                F.mul(w[i1 * db + jp], B[j][jp])
                for jp in range(db)
                if not F.is_zero(B[j][jp])
            )
            for j in range(db)
        ]
        for i1 in range(da)
    ]
    out = [F.zero] * (da * db)
    for i in range(da):
        for j in range(db):
            out[i * db + j] = F.total(
                F.mul(A[i][ip], X[ip][j]) for ip in range(da) if not F.is_zero(A[i][ip])
            )
    return out


def _swap_index(j: int, d: int) -> int:
    i1, i2 = divmod(j, d)
    return i2 * d + i1


def couple(
    model: Model,
    alpha: Irrep,
    beta: Irrep,
    Fcg: CyclotomicField,
    table: NormFactorTable,
) -> Pair:
    """The coupling matrix of ``(alpha, beta)`` in the canonical basis."""
    F = model.field
    order = model.order
    da, db = alpha.dimension, beta.dimension
    same = alpha is beta
    copies: list[Copy] = []
    mults: dict[int, tuple[int, int, int]] = {}
    for gamma in model.irreps:
        if gamma.projective != (alpha.projective != beta.projective):
            continue
        n, n_plus, n_minus = multiplicities(model, alpha, beta, gamma)
        if n == 0:
            continue
        mults[gamma.index] = (n, n_plus, n_minus)
        dg = gamma.dimension
        weights00 = [F.conjugate(gamma.matrices[g][0][0]) for g in range(order)]
        diag_cache: dict[int, Value] = {}

        def P_diag(j: int, weights00=weights00, diag_cache=diag_cache, dg=dg) -> Value:
            if j not in diag_cache:
                j1, j2 = divmod(j, db)
                total = F.total(
                    F.mul(
                        weights00[g],
                        F.mul(alpha.matrices[g][j1][j1], beta.matrices[g][j2][j2]),
                    )
                    for g in range(order)
                    if not F.is_zero(weights00[g])
                )
                diag_cache[j] = F.scale(total, Fraction(dg, order))
            return diag_cache[j]

        column_cache: dict[int, Vector] = {}

        def P_column(
            j: int, weights00=weights00, column_cache=column_cache, dg=dg
        ) -> Vector:
            if j not in column_cache:
                column_cache[j] = _product_column(
                    F, alpha, beta, weights00, j, dg, order
                )
            return column_cache[j]

        parts = [(+1, n_plus), (-1, n_minus)] if same else [(0, n)]
        r = 0
        for sign, count in parts:
            found: list[tuple[Vector, Value]] = []
            for _ in range(count):
                # Q_jj = P_part[j][j] - sum_s |w_s[j]|^2 / N_s; find the first non-zero
                j0 = None
                for j in range(da * db):
                    q = P_diag(j)
                    if sign:
                        q = F.scale(
                            F.add(q, F.scale(P_column(j)[_swap_index(j, da)], sign)),
                            Fraction(1, 2),
                        )
                    for ws, Ns in found:
                        q = F.sub(q, F.divide(F.norm_squared(ws[j]), Ns))
                    if not F.is_zero(q):
                        j0 = j
                        break
                if j0 is None:
                    raise ArithmeticError(
                        f"no trial state for {gamma.symbol.text} in {alpha.symbol.text} x {beta.symbol.text}"
                    )
                w = P_column(j0)
                if sign:
                    swapped = [w[_swap_index(j, da)] for j in range(da * db)]
                    w = [
                        F.scale(F.add(x, F.scale(y, sign)), Fraction(1, 2))
                        for x, y in zip(w, swapped, strict=True)
                    ]
                for ws, Ns in found:
                    factor = F.divide(F.conjugate(ws[j0]), Ns)
                    w = [F.sub(x, F.mul(factor, y)) for x, y in zip(w, ws, strict=True)]
                N = w[j0]
                if N != F.total(F.norm_squared(x) for x in w):
                    raise ArithmeticError(
                        "the designated coordinate is not the squared norm"
                    )
                if any(not F.is_zero(w[j]) for j in range(j0)):
                    raise ArithmeticError(
                        "a coordinate before the trial state is non-zero"
                    )
                found.append((w, N))
                # the other components by the transfer operators P^gamma_{k 0}
                applied = [_product_apply(F, alpha, beta, g, w) for g in range(order)]
                columns_G: list[Vector] = [w]
                for k in range(1, dg):
                    wk = [F.zero] * (da * db)
                    for g in range(order):
                        c = F.conjugate(gamma.matrices[g][k][0])
                        if F.is_zero(c):
                            continue
                        wk = [
                            F.add(x, F.mul(c, y))
                            for x, y in zip(wk, applied[g], strict=True)
                        ]
                    columns_G.append([F.scale(x, Fraction(dg, order)) for x in wk])
                # normalisation in the cg field
                NM = F.lift(N, Fcg)
                mu = norm_factor(Fcg, NM, table)
                if mu is None:
                    raise ArithmeticError(
                        f"no norm factor for N in Q(zeta_{Fcg.conductor})"
                    )
                inv_mu = Fcg.inverse(mu)
                columns = [
                    [Fcg.mul(F.lift(x, Fcg), inv_mu) for x in col] for col in columns_G
                ]
                copies.append(
                    Copy(gamma.index, r, sign if same else None, j0, N, columns)
                )
                r += 1
    return Pair(alpha.index, beta.index, copies, mults)


def transform(
    pair: Pair, Fcg: CyclotomicField, basis: ChainBasis, irreps: list[Irrep]
) -> Pair:
    """The unitary transform of a canonical-basis pair into the chain-adapted basis:
    ``U' = (U_alpha x U_beta)^dagger U U_gamma`` per copy."""
    Fb = basis.field
    Ua = [[Fb.lift(x, Fcg) for x in row] for row in basis.irreps[pair.alpha].W]
    Ub = [[Fb.lift(x, Fcg) for x in row] for row in basis.irreps[pair.beta].W]
    da, db = len(Ua), len(Ub)
    copies = []
    for copy in pair.copies:
        Ug = [[Fb.lift(x, Fcg) for x in row] for row in basis.irreps[copy.gamma].W]
        dg = len(Ug)
        # rows: (Ua x Ub)^dagger applied to every column
        rotated = []
        for col in copy.columns:
            W = [[col[i * db + j] for j in range(db)] for i in range(da)]
            # X = Ua^dagger W conj(Ub)
            X = [
                [
                    Fcg.total(
                        Fcg.mul(Fcg.conjugate(Ua[i][ip]), W[i][j]) for i in range(da)
                    )
                    for j in range(db)
                ]
                for ip in range(da)
            ]
            Y = [
                [
                    Fcg.total(
                        Fcg.mul(X[ip][j], Fcg.conjugate(Ub[j][jp])) for j in range(db)
                    )
                    for jp in range(db)
                ]
                for ip in range(da)
            ]
            rotated.append([Y[ip][jp] for ip in range(da) for jp in range(db)])
        # columns: C U_gamma
        columns = [
            [
                Fcg.total(Fcg.mul(rotated[k][row], Ug[k][kp]) for k in range(dg))
                for row in range(da * db)
            ]
            for kp in range(dg)
        ]
        copies.append(
            Copy(copy.gamma, copy.r, copy.exchange, copy.trial, copy.N, columns)
        )
    return Pair(pair.alpha, pair.beta, copies, pair.multiplicities)


# ---- checks ------------------------------------------------------------------------------


def check_pair(
    model: Model,
    pair: Pair,
    Fcg: CyclotomicField,
    table: NormFactorTable,
    matrices: list[list[Matrix]] | None = None,
) -> list[str]:
    """CG-1 unitarity, CG-2 reduction on a generating set, CG-3 counts, CG-4
    exchange eigenvectors, CG-5 the designated coefficient (canonical basis
    only). ``matrices`` overrides the irrep matrices (chain basis), in ``Fcg``."""
    irreps = model.irreps
    alpha, beta = irreps[pair.alpha], irreps[pair.beta]
    name = f"{model.spec.symbol} {alpha.symbol.text} x {beta.symbol.text}"
    failures: list[str] = []
    U = pair.matrix()
    n = len(U)
    F0 = model.field

    def D(idx: int, g: int) -> Matrix:
        if matrices is not None:
            return matrices[idx][g]
        return [[F0.lift(x, Fcg) for x in row] for row in irreps[idx].matrices[g]]

    # CG-1
    for c1 in range(n):
        for c2 in range(c1, n):
            inner = Fcg.total(
                Fcg.mul(Fcg.conjugate(U[i][c1]), U[i][c2]) for i in range(n)
            )
            if inner != (Fcg.one if c1 == c2 else Fcg.zero):
                failures.append(f"{name}: CG-1 fails at columns {c1}, {c2}")
                break
        else:
            continue
        break
    # CG-3
    if sum(m[0] * irreps[g].dimension for g, m in pair.multiplicities.items()) != n:
        failures.append(f"{name}: CG-3 dimensions")
    # CG-2 on a generating set
    for g in generating_set(model.table):
        Da, Db = D(pair.alpha, g), D(pair.beta, g)
        da, db = len(Da), len(Db)
        col = 0
        for copy in pair.copies:
            Dg = D(copy.gamma, g)
            dg = len(Dg)
            for k in range(dg):
                # (Da x Db) U[:, col + k] == sum_k' U[:, col + k'] Dg[k'][k]
                v = copy.columns[k]
                W = [[v[i * db + j] for j in range(db)] for i in range(da)]
                X = [
                    [
                        Fcg.total(Fcg.mul(W[i][jp], Db[j][jp]) for jp in range(db))
                        for j in range(db)
                    ]
                    for i in range(da)
                ]
                lhs = [
                    Fcg.total(Fcg.mul(Da[i][ip], X[ip][j]) for ip in range(da))
                    for i in range(da)
                    for j in range(db)
                ]
                rhs = [
                    Fcg.total(
                        Fcg.mul(copy.columns[kp][row], Dg[kp][k]) for kp in range(dg)
                    )
                    for row in range(da * db)
                ]
                if lhs != rhs:
                    failures.append(
                        f"{name}: CG-2 fails for {irreps[copy.gamma].symbol.text} copy {copy.r} component {k} at element {g}"
                    )
                    break
            col += dg
    # CG-4 and CG-5
    d = alpha.dimension
    for copy in pair.copies:
        if copy.exchange is not None:
            for v in copy.columns:
                swapped = [v[_swap_index(j, d)] for j in range(n)]
                if swapped != [Fcg.scale(x, copy.exchange) for x in v]:
                    failures.append(
                        f"{name}: CG-4 fails for {irreps[copy.gamma].symbol.text} copy {copy.r}"
                    )
                    break
        if matrices is None:
            v = copy.columns[0]
            first = next(x for x in v if not Fcg.is_zero(x))
            mu = norm_factor(Fcg, Fcg.norm_squared(first), table)
            if mu is None or first != Fcg.conjugate(mu):
                failures.append(
                    f"{name}: CG-5 fails for {irreps[copy.gamma].symbol.text} copy {copy.r}"
                )
    return failures


# ---- the cg table ------------------------------------------------------------------------


def pair_record(
    pair: Pair, irreps: list[Irrep], Fcg: CyclotomicField, forms: tuple[str, ...]
) -> dict[str, Any]:
    from pointgroup_tables.values import numeric, pretty

    U = pair.matrix()
    record: dict[str, Any] = {
        "irreps": [pair.alpha, pair.beta],
        "dimension": pair.dimension,
        "constituents": pair.constituent_records(irreps),
        "columns": pair.column_records(irreps),
    }
    if "exact" in forms:
        record["exact"] = [[Fcg.exact(x) for x in row] for row in U]
    if "numeric" in forms:
        record["numeric"] = [[numeric(Fcg, x) for x in row] for row in U]
    if "pretty" in forms:
        record["pretty"] = [[pretty(Fcg, x) for x in row] for row in U]
    return record


def cg_table(
    model: Model,
    pairs: list[Pair],
    Fcg: CyclotomicField,
    forms: tuple[str, ...],
    basis: ChainBasis | None = None,
) -> dict[str, Any]:
    records = []
    for irrep in model.irreps:
        r = irrep.record()
        if basis is not None:
            r["components"] = components_of(basis, basis.irreps[irrep.index])
        records.append(r)
    return {
        "irreps": records,
        "pairs": [pair_record(p, model.irreps, Fcg, forms) for p in pairs],
    }


def select_pairs(model: Model, which: str) -> list[tuple[Irrep, Irrep]]:
    out = []
    irreps = model.irreps
    for a in range(len(irreps)):
        for b in range(a, len(irreps)):
            alpha, beta = irreps[a], irreps[b]
            projective = alpha.projective or beta.projective
            if which == "ordinary" and projective:
                continue
            if which == "projective" and not projective:
                continue
            out.append((alpha, beta))
    return out


__all__ = [
    "Copy",
    "Pair",
    "cg_conductor",
    "cg_table",
    "check_pair",
    "couple",
    "cyclotomic_field",
    "select_pairs",
    "transform",
]
