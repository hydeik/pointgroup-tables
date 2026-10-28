"""Clebsch-Gordan multiplicity facts: Sym/Alt split, norms under the proposed rule,
sizes of the full CG unitaries.

Run from the repo root: ``uv run python3 <scratch>/cgfacts/cg_multiplicity.py``.
Helpers are imported from the copied research script ``cyclotomic_closure.py``
in the same directory.
"""

from __future__ import annotations

import math
import os
import sys
import time
from dataclasses import dataclass
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cyclotomic_closure import (
    _projector_entries,
    analyse_norm,
    cg_norms,
    embed,
    multiplicity,
)

from pointgroup_tables.document import Model, build_model
from pointgroup_tables.field import CyclotomicField, Matrix, Value
from pointgroup_tables.irreps import Irrep
from pointgroup_tables.symbols import parse_group
from pointgroup_tables.values import pretty

MULT_GROUPS = ["T", "O", "I", "D6h", "C3v"]
SIZE_GROUPS = ["C3v", "D6h", "Oh", "Ih"]


# ---- helpers -------------------------------------------------------------------------


def compatible(alpha: Irrep, beta: Irrep, gamma: Irrep) -> bool:
    return gamma.projective == (alpha.projective != beta.projective)


def trace_of_square(F: CyclotomicField, D: Matrix) -> Value:
    d = len(D)
    return F.total(
        F.mul(D[i][j], D[j][i])
        for i in range(d)
        for j in range(d)
        if not F.is_zero(D[i][j]) and not F.is_zero(D[j][i])
    )


def rational_int(F: CyclotomicField, x: Value, what: str) -> int:
    if not F.is_rational(x) or x[0].denominator != 1:
        raise ArithmeticError(f"{what}: non-integral value {x}")
    return int(x[0])


def sym_alt_multiplicities(model: Model, alpha: Irrep, gamma: Irrep) -> tuple[int, int]:
    """``(n_sym, n_alt)`` of ``gamma`` in ``Sym^2`` and ``Alt^2`` of ``alpha``, from
    ``tr(D(g))^2`` and ``tr(D(g)^2)`` of the irrep matrices (valid for projective
    irreps too, the sum over ``G`` with the canonical section)."""
    F = model.field
    s_total = F.zero
    a_total = F.zero
    for g in range(model.order):
        chi = alpha.characters[g]
        chi2 = F.mul(chi, chi)
        t2 = trace_of_square(F, alpha.matrices[g])
        cg = F.conjugate(gamma.characters[g])
        if F.is_zero(cg):
            continue
        s_total = F.add(s_total, F.mul(F.add(chi2, t2), cg))
        a_total = F.add(a_total, F.mul(F.sub(chi2, t2), cg))
    scale = Fraction(1, 2 * model.order)
    n_sym = rational_int(F, F.scale(s_total, scale), "n_sym")
    n_alt = rational_int(F, F.scale(a_total, scale), "n_alt")
    return n_sym, n_alt


def swap_index(j: tuple[int, int]) -> tuple[int, int]:
    return (j[1], j[0])


@dataclass
class CGVector:
    alpha: Irrep
    beta: Irrep
    gamma: Irrep
    copy: int  # 1-based, over all parts
    part: str  # "sym", "alt" or "" (alpha != beta)
    component: tuple[int, int]
    vector: list[Value]  # unnormalised partner-0 vector, designated component = N
    norm: Value


def cg_vectors_pair(
    model: Model, alpha: Irrep, beta: Irrep, gamma: Irrep, m: int
) -> list[CGVector]:
    """The unnormalised partner-0 CG vectors of ``gamma`` in ``alpha x beta`` under
    the proposed rule: Sym then Alt for ``alpha == beta`` (Gram-Schmidt within a
    part, row-major trial order), plain Gram-Schmidt in row-major order otherwise.
    Designated component of ``gamma``: index 0."""
    F = model.field
    da, db = alpha.dimension, beta.dimension
    components = [(j1, j2) for j1 in range(da) for j2 in range(db)]
    out: list[CGVector] = []
    if alpha is beta:
        n_sym, n_alt = sym_alt_multiplicities(model, alpha, gamma)
        if n_sym + n_alt != m:
            raise ArithmeticError("n_sym + n_alt != n_gamma")
        parts = [("sym", n_sym, 1), ("alt", n_alt, -1)]
    else:
        parts = [("", m, 0)]
    half = Fraction(1, 2)
    for part, count, sign in parts:
        found: list[tuple[list[Value], Value]] = []
        for j in components:
            if len(found) == count:
                break
            if sign == 0:
                diagonal = _projector_entries(model, alpha, beta, gamma, [j], j)[0]
            else:
                pj, ps = _projector_entries(
                    model, alpha, beta, gamma, [j, swap_index(j)], j
                )
                diagonal = F.scale(F.add(pj, ps) if sign > 0 else F.sub(pj, ps), half)
            if F.is_zero(diagonal):
                continue
            w = _projector_entries(model, alpha, beta, gamma, components, j)
            if sign != 0:
                w = [
                    F.scale(
                        F.add(w[components.index(swap_index(c))], x)
                        if sign > 0
                        else F.sub(x, w[components.index(swap_index(c))]),
                        half,
                    )
                    for c, x in zip(components, w, strict=True)
                ]
            for u, nu in found:
                overlap = F.total(
                    F.mul(F.conjugate(x), y) for x, y in zip(u, w, strict=True)
                )
                if F.is_zero(overlap):
                    continue
                factor = F.divide(overlap, nu)
                w = [F.sub(y, F.mul(factor, x)) for x, y in zip(u, w, strict=True)]
            if all(F.is_zero(x) for x in w):
                continue
            norm = w[components.index(j)]
            if norm != F.total(F.norm_squared(x) for x in w):
                raise ArithmeticError("designated component is not the squared norm")
            found.append((w, norm))
            out.append(CGVector(alpha, beta, gamma, len(out) + 1, part, j, w, norm))
        if len(found) != count:
            raise ArithmeticError(
                f"UNEXPECTED: found {len(found)} of {count} copies of "
                f"{gamma.symbol.text} in {alpha.symbol.text} x {beta.symbol.text} "
                f"part {part!r}: zero diagonal for every trial state"
            )
    return out


def apply_product(
    F: CyclotomicField, Da: Matrix, Db: Matrix, w: list[Value], da: int, db: int
) -> list[Value]:
    """``(Da x Db) w`` with ``w`` indexed row-major ``(i1, i2)``."""
    # W[j1][j2] -> (Da W Db^T)[i1][i2]
    W = [[w[j1 * db + j2] for j2 in range(db)] for j1 in range(da)]
    left = [
        [
            F.total(
                F.mul(Da[i1][j1], W[j1][j2])
                for j1 in range(da)
                if not F.is_zero(Da[i1][j1]) and not F.is_zero(W[j1][j2])
            )
            for j2 in range(db)
        ]
        for i1 in range(da)
    ]
    out: list[Value] = []
    for i1 in range(da):
        for i2 in range(db):
            out.append(
                F.total(
                    F.mul(left[i1][j2], Db[i2][j2])
                    for j2 in range(db)
                    if not F.is_zero(left[i1][j2]) and not F.is_zero(Db[i2][j2])
                )
            )
    return out


def matrix_inverse(F: CyclotomicField, A: Matrix) -> Matrix:
    n = len(A)
    rows = [
        list(A[i]) + [F.one if i == j else F.zero for j in range(n)] for i in range(n)
    ]
    for col in range(n):
        pivot = next(r for r in range(col, n) if not F.is_zero(rows[r][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        inv = F.inverse(rows[col][col])
        rows[col] = [F.mul(inv, x) for x in rows[col]]
        for r in range(n):
            if r != col and not F.is_zero(rows[r][col]):
                f = rows[r][col]
                rows[r] = [
                    F.sub(x, F.mul(f, y))
                    for x, y in zip(rows[r], rows[col], strict=True)
                ]
    return [row[n:] for row in rows]


_transfer_cache: dict[tuple[int, int], tuple[list[int], Matrix]] = {}


def transfer_data(model: Model, gamma: Irrep) -> tuple[list[int], Matrix]:
    """Elements ``g_i`` with ``[D_gamma_k0(g_i)]`` invertible and the inverse
    ``Minv`` with ``w_k = sum_i Minv[k][i] (D x D)(g_i) w_0``."""
    key = (id(model), gamma.index if gamma.index >= 0 else id(gamma))
    if key in _transfer_cache:
        return _transfer_cache[key]
    F = model.field
    d = gamma.dimension
    chosen: list[int] = []
    rows: list[list[Value]] = []
    for g in range(model.order):
        candidate = rows + [[gamma.matrices[g][k][0] for k in range(d)]]
        if _rank(F, candidate) == len(candidate):
            rows = candidate
            chosen.append(g)
            if len(chosen) == d:
                break
    if len(chosen) != d:
        raise ArithmeticError("no invertible transfer system")
    M = rows  # M[i][k] = D_k0(g_i)
    Minv = matrix_inverse(F, M)  # Minv[k][i]
    _transfer_cache[key] = (chosen, Minv)
    return chosen, Minv


def _rank(F: CyclotomicField, A: list[list[Value]]) -> int:
    rows = [list(r) for r in A]
    r = 0
    cols = len(A[0]) if A else 0
    for c in range(cols):
        pivot = next(
            (i for i in range(r, len(rows)) if not F.is_zero(rows[i][c])), None
        )
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inv = F.inverse(rows[r][c])
        rows[r] = [F.mul(inv, x) for x in rows[r]]
        for i in range(len(rows)):
            if i != r and not F.is_zero(rows[i][c]):
                f = rows[i][c]
                rows[i] = [
                    F.sub(x, F.mul(f, y)) for x, y in zip(rows[i], rows[r], strict=True)
                ]
        r += 1
    return r


def partners(model: Model, v: CGVector) -> list[list[Value]]:
    """All partner vectors ``w_k = P_k0 w_0`` (``k = 0 .. d_gamma - 1``), via the
    transfer relation ``(D x D)(g) w_0 = sum_k D_k0(g) w_k``."""
    F = model.field
    gamma = v.gamma
    d = gamma.dimension
    if d == 1:
        return [v.vector]
    chosen, Minv = transfer_data(model, gamma)
    da, db = v.alpha.dimension, v.beta.dimension
    images = [
        apply_product(F, v.alpha.matrices[g], v.beta.matrices[g], v.vector, da, db)
        for g in chosen
    ]
    n = len(v.vector)
    out = []
    for k in range(d):
        wk = [
            F.total(
                F.mul(Minv[k][i], images[i][r])
                for i in range(d)
                if not F.is_zero(Minv[k][i]) and not F.is_zero(images[i][r])
            )
            for r in range(n)
        ]
        out.append(wk)
    return out


def partner_direct(model: Model, v: CGVector, k: int) -> list[Value]:
    """``P_k0 w_0`` by the group sum (for verification)."""
    F = model.field
    alpha, beta, gamma = v.alpha, v.beta, v.gamma
    da, db = alpha.dimension, beta.dimension
    n = da * db
    total = [F.zero] * n
    for g in range(model.order):
        weight = F.conjugate(gamma.matrices[g][k][0])
        if F.is_zero(weight):
            continue
        image = apply_product(F, alpha.matrices[g], beta.matrices[g], v.vector, da, db)
        total = [F.add(t, F.mul(weight, x)) for t, x in zip(total, image, strict=True)]
    scale = Fraction(gamma.dimension, model.order)
    return [F.scale(x, scale) for x in total]


# ---- per-group analyses ---------------------------------------------------------------


def all_pairs(model: Model):
    irreps = model.irreps
    for ia, alpha in enumerate(irreps):
        for beta in irreps[ia:]:
            yield alpha, beta


def triples(model: Model):
    for alpha, beta in all_pairs(model):
        for gamma in model.irreps:
            if not compatible(alpha, beta, gamma):
                continue
            m = multiplicity(model, alpha, beta, gamma)
            if m:
                yield alpha, beta, gamma, m


def norm_report(
    K: CyclotomicField, grouped: dict[Value, list[str]], label: str
) -> None:
    closure = K.conductor
    closure_free = K.conductor
    in_field = 0
    factors = 0
    sevens: list[str] = []
    for value in sorted(grouped, key=lambda v: embed(K, v).real):
        a = analyse_norm(K.conductor, value)
        examples = grouped[value]
        tag = (
            "[sqrt in K]"
            if a.in_field
            else (
                f"[sqrt in Q(zeta_{a.conductor}), not in K]"
                if a.conductor
                else "[not cyclotomic / beyond bound]"
            )
        )
        mu = "mu in K found" if a.hermitian_factor is not None else "NO mu in K found"
        print(
            f"    N = {pretty(K, value):<28} {tag:<36} {mu:<16} {len(examples)}x  e.g. {examples[0]}"
        )
        if a.conductor is None:
            closure = 0
        else:
            closure = math.lcm(closure, a.conductor) if closure else 0
            if a.conductor % 7 == 0:
                sevens.extend(examples)
            free = [e for e in examples if e.endswith(" [mult-free]")]
            if free and closure_free:
                closure_free = math.lcm(closure_free, a.conductor)
        in_field += a.in_field
        factors += a.hermitian_factor is not None
    print(
        f"  => {label}: {len(grouped)} distinct N; sqrt in K for {in_field}; "
        f"Hermitian factor mu in K found for {factors} of {len(grouped)}; "
        f"closure conductor (all triples) = {closure or 'none'}; "
        f"closure conductor (multiplicity-free triples only) = {closure_free or 'none'}"
    )
    if sevens:
        print(f"  => conductor divisible by 7 arises from: {sorted(set(sevens))}")
    else:
        print("  => no conductor divisible by 7")


def analyse_group(name: str, want_sizes: bool) -> None:
    started = time.perf_counter()
    model = build_model(parse_group(name))
    F = model.field
    irreps = model.irreps
    print(f"\n{'=' * 90}\n## {name}  (|G| = {model.order}, K = Q(zeta_{F.conductor}))")
    print(
        "irreps (canonical order): "
        + ", ".join(
            f"{i.symbol.text}({i.dimension}{'p' if i.projective else ''})"
            for i in irreps
        )
    )

    # ---- item 1: multiplicities > 1
    all_triples = list(triples(model))
    multi = [(a, b, g, m) for a, b, g, m in all_triples if m > 1]
    print(
        f"\n### 1. Triples with n_gamma > 1 ({len(multi)} of {len(all_triples)} triples)"
    )
    for a, b, g, m in multi:
        print(
            f"  {a.symbol.text} x {b.symbol.text} > {g.symbol.text}: n = {m}  (alpha == beta: {a is b})"
        )
    if not multi:
        print("  (none)")

    # ---- item 2: Sym/Alt split
    print("\n### 2. Sym^2 / Alt^2 split for alpha == beta")
    counts = {"sym": 0, "alt": 0, "mixed": 0}
    unresolved = []
    full_lines = []
    for a, b, g, m in all_triples:
        if a is not b:
            continue
        ns, na = sym_alt_multiplicities(model, a, g)
        if ns + na != m:
            raise ArithmeticError("n_sym + n_alt != n_gamma")
        where = "Sym" if na == 0 else ("Alt" if ns == 0 else "both")
        if where == "both":
            counts["mixed"] += 1
        else:
            counts[where.lower()] += 1
        full_lines.append(
            f"  {a.symbol.text} x {a.symbol.text} > {g.symbol.text}: n = {m}, n_sym = {ns}, n_alt = {na}  -> {where}"
        )
        if max(ns, na) >= 2:
            unresolved.append((a.symbol.text, g.symbol.text, m, ns, na))
    if name in ("T", "O", "D6h", "C3v"):
        print("\n".join(full_lines))
    else:
        for line in full_lines:
            if "n = 1," not in line:
                print(line)
        print(
            f"  (n = 1 triples: {counts['sym']} purely Sym, {counts['alt']} purely Alt listed in counts below)"
        )
    print(
        f"  summary: {counts['sym']} triples purely in Sym, {counts['alt']} purely in Alt, {counts['mixed']} in both"
    )
    multi_eq = [(a, g, m) for a, b, g, m in multi if a is b]
    if multi_eq:
        if unresolved:
            print(
                "  NOT fully separated by the split (a part still has multiplicity >= 2):"
            )
            for a, g, m, ns, na in unresolved:
                print(f"    {a} x {a} > {g}: n = {m}, n_sym = {ns}, n_alt = {na}")
        else:
            print(
                "  every alpha == beta triple with n > 1 is fully separated by the Sym/Alt split (max part multiplicity 1)"
            )
    else:
        print("  (no alpha == beta triple with n > 1)")

    # ---- item 5 for O
    if name == "O":
        t1 = next(i for i in irreps if i.symbol.text == "T1")
        print("\n### 5. Sanity check O: T1 x T1")
        for g in irreps:
            if compatible(t1, t1, g) and multiplicity(model, t1, t1, g):
                ns, na = sym_alt_multiplicities(model, t1, g)
                print(f"  {g.symbol.text}: n_sym = {ns}, n_alt = {na}")

    # ---- item 3: norms under the proposed rule
    print(
        "\n### 3. Squared norms under the proposed rule (Sym first, then Alt; designated gamma component 0)"
    )
    vectors: dict[tuple[int, int], list[CGVector]] = {}
    grouped: dict[Value, list[str]] = {}
    t0 = time.perf_counter()
    for a, b, g, m in all_triples:
        vs = cg_vectors_pair(model, a, b, g, m)
        vectors.setdefault((a.index, b.index), []).extend(vs)
        for v in vs:
            tag = f"{a.symbol.text} x {b.symbol.text} > {g.symbol.text}"
            if v.part:
                tag += f" [{v.part} #{v.copy}]"
            elif m > 1:
                tag += f" [#{v.copy}]"
            if m == 1:
                tag += " [mult-free]"
            grouped.setdefault(v.norm, []).append(tag)
    nvec = sum(len(v) for v in vectors.values())
    print(
        f"  {nvec} partner-0 vectors over {len(vectors)} pairs ({time.perf_counter() - t0:.1f} s)"
    )
    norm_report(F, grouped, "proposed rule")
    # the multiplicity > 1 vectors in detail
    print("  norms of the multiplicity > 1 triples under the proposed rule:")
    for a, b, g, m in multi:
        vs = [v for v in vectors[(a.index, b.index)] if v.gamma is g]
        for v in vs:
            print(
                f"    {a.symbol.text} x {b.symbol.text} > {g.symbol.text} copy {v.copy} {v.part or '-':<4} e_j = {v.component}: N = {pretty(F, v.norm)}"
            )

    # original rule for comparison
    print(
        "  original rule of cyclotomic_closure.cg_norms (Gram-Schmidt in product order, no split):"
    )
    old = cg_norms(model)
    grouped_old: dict[Value, list[str]] = {}
    mults = {
        (a.symbol.text, b.symbol.text, g.symbol.text): m for a, b, g, m in all_triples
    }
    for r in old:
        tag = f"{r.alpha} x {r.beta} > {r.gamma}"
        if r.copy > 1 or mults[(r.alpha, r.beta, r.gamma)] > 1:
            tag += f" [#{r.copy}]"
        else:
            tag += " [mult-free]"
        grouped_old.setdefault(r.norm, []).append(tag)
    norm_report(F, grouped_old, "original rule")
    only_new = set(grouped) - set(grouped_old)
    only_old = set(grouped_old) - set(grouped)
    print(
        f"  N values only under the proposed rule: {[pretty(F, v) for v in only_new]}"
    )
    print(
        f"  N values only under the original rule: {[pretty(F, v) for v in only_old]}"
    )

    # ---- item 4: sizes and non-zero counts
    if want_sizes:
        print("\n### 4. Sizes of the full CG unitaries")
        t0 = time.perf_counter()
        n_irreps = len(irreps)
        n_pairs = n_irreps * (n_irreps + 1) // 2
        n_ordered = n_irreps * n_irreps
        entries = 0
        nonzero = 0
        entries_ord = 0
        nonzero_ord = 0
        pairs_ord = 0
        verified_pairs = 0
        for a, b in all_pairs(model):
            da, db = a.dimension, b.dimension
            n = da * db
            vs = vectors.get((a.index, b.index), [])
            columns: list[list[Value]] = []
            for v in vs:
                columns.extend(partners(model, v))
            if len(columns) != n:
                raise ArithmeticError(
                    f"{a.symbol.text} x {b.symbol.text}: {len(columns)} columns for dimension {n}"
                )
            nz = sum(1 for col in columns for x in col if not F.is_zero(x))
            entries += n * n
            nonzero += nz
            if not a.projective and not b.projective:
                pairs_ord += 1
                entries_ord += n * n
                nonzero_ord += nz
            # exact check: columns pairwise orthogonal, norms = N of their vector
            if n <= 16 or name != "Ih":
                norms = [v.norm for v in vs for _ in range(v.gamma.dimension)]
                for i, ci in enumerate(columns):
                    for k in range(i, len(columns)):
                        ck = columns[k]
                        ip = F.total(
                            F.mul(F.conjugate(x), y)
                            for x, y in zip(ci, ck, strict=True)
                            if not F.is_zero(x) and not F.is_zero(y)
                        )
                        expected = norms[i] if i == k else F.zero
                        if ip != expected:
                            raise ArithmeticError(
                                f"{a.symbol.text} x {b.symbol.text}: Gram check failed at ({i}, {k})"
                            )
                verified_pairs += 1
        print(
            f"  irreps: {n_irreps}; unordered pairs (alpha <= beta): {n_pairs}; ordered pairs: {n_ordered}"
        )
        print(
            f"  all pairs: sum (d_a d_b)^2 = {entries}; non-zero entries = {nonzero}  ({Fraction(nonzero, entries)} = {nonzero / entries:.4f})"
        )
        print(
            f"  ordinary x ordinary: pairs = {pairs_ord}; sum (d_a d_b)^2 = {entries_ord}; non-zero entries = {nonzero_ord}  ({Fraction(nonzero_ord, entries_ord)} = {nonzero_ord / entries_ord:.4f})"
        )
        print(
            f"  Gram (orthogonality + norm) check passed for {verified_pairs} of {n_pairs} pairs ({time.perf_counter() - t0:.1f} s)"
        )

    print(f"\n({name} done in {time.perf_counter() - started:.1f} s)")


def verify_transfer(name: str, limit: int) -> None:
    """Check the transfer-relation partners against the direct projector sum."""
    model = build_model(parse_group(name))
    checked = 0
    for a, b, g, m in triples(model):
        if g.dimension == 1:
            continue
        for v in cg_vectors_pair(model, a, b, g, m):
            ws = partners(model, v)
            for k in range(g.dimension):
                if ws[k] != partner_direct(model, v, k):
                    raise ArithmeticError(
                        f"transfer mismatch {name}: {a.symbol.text} x {b.symbol.text} > {g.symbol.text} k={k}"
                    )
            checked += 1
            if checked >= limit:
                print(
                    f"  transfer relation verified against P_k0 for {checked} vectors of {name}"
                )
                return
    print(f"  transfer relation verified against P_k0 for {checked} vectors of {name}")


def main() -> int:
    print("# CG multiplicity facts (cg_multiplicity.py)")
    print(
        "Designated component of gamma in cyclotomic_closure._projector_entries: index 0 "
        "(weights = conj(D_gamma[g][0][0])), the same as the proposed rule; no change needed."
    )
    print("\n### Transfer-operator verification")
    verify_transfer("T", 50)
    verify_transfer("O", 50)
    verify_transfer("I", 12)
    for name in MULT_GROUPS:
        analyse_group(name, want_sizes=name in SIZE_GROUPS)
    for name in SIZE_GROUPS:
        if name not in MULT_GROUPS:
            analyse_group(name, want_sizes=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
