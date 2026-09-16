"""A basis convention applied to a group (docs/spec/basis-conventions.md):
the conductor of the document field, the coverage, the generator images of
the axial irreps under the convention's ``[axial]`` rule, the lifted
generating set of the polyhedral data, and the intertwiner with the canonical
realisation that the verification computes and never stores.

The canonical basis and every reference basis convention go through the
same functions; ``conventions.py`` reads the files, ``irreps.py`` closes the
generator images over the binary cover.
"""

from __future__ import annotations

from typing import Any

from . import conventions
from .abstract import PointGroup, point_group
from .double import BinaryRotationGroup
from .errors import CONVENTIONS_HINT, UnknownConventionError
from .field import (
    CyclotomicField,
    Matrix,
    Value,
    canonical_conductor,
    matrix_identity,
    matrix_mul,
)
from .labels import canonical_labels
from .symbols import GroupSpec, parse_group

GeneratorImages = tuple[Matrix, Matrix]  # D(a), D(b) of a two-dimensional axial irrep


# ---- the convention for a group ---------------------------------------------------------


def cover_name(G: PointGroup) -> str | None:
    """``2T``, ``2O`` or ``2I`` for a polyhedral group, ``None`` for an axial one."""
    R = G.rotation_group
    return f"2{R.name}" if R.polyhedral else None


def _uncovered(name: str, G: PointGroup) -> UnknownConventionError:
    return UnknownConventionError(
        f"basis convention '{name}' does not cover {G.spec.symbol}", CONVENTIONS_HINT
    )


def rules_for(name: str, G: PointGroup) -> conventions.BasisConvention:
    """The rules of a basis convention, provided it realises the group; exit code 4
    otherwise (a convention never falls back to the canonical basis silently)."""
    rules = conventions.basis_convention_rules(name)
    family_n = None if G.spec.polyhedral else G.spec.family_n
    if not rules.covers_group(family_n, cover_name(G)):
        raise _uncovered(name, G)
    return rules


def document_conductor(name: str, G: PointGroup) -> int:
    """The conductor of a document of ``G`` under the convention: the declared one
    of the polyhedral data, or the canonical conductor for the axial families
    (the canonical rule and the spin-harmonic rule stay in ``Q(zeta_lcm(2N, 4))``)."""
    rules_for(name, G)
    cover = cover_name(G)
    if cover is None:
        return canonical_conductor(G.rotation_group)
    return conventions.polyhedral_data(name, cover).conductor


def metadata(name: str, G: PointGroup) -> dict[str, Any]:
    """The ``basis_convention`` block: ``{ name }`` for the canonical basis,
    ``{ name, source, coverage }`` for a reference convention."""
    block = conventions.basis_convention(name)
    rules = rules_for(name, G)
    if rules.reference:
        family_n = None if G.spec.polyhedral else G.spec.family_n
        block["coverage"] = "reference" if rules.verified(family_n) else "extrapolated"
    return block


# ---- the spin-harmonic rule -------------------------------------------------------------
#
# Altmann & Herzig's harmonic bases <|j m1>, +-|j -m1>| in closed form
# (docs/spec/basis-conventions.md, "The spin-harmonic rule"; ADR-0008). Exponents
# are kept in halves: ``k2 = 2k`` for an ordinary ``E_k``, ``k2 = k`` for a
# projective ``E_{k/2}``; ``m2 = 2 m1`` is the signed first partner.


def _residue(m2: int, L: int) -> int:
    """The representative of ``m2/2`` modulo ``L`` in ``(-L/2, L/2]``, in halves."""
    r = m2 % (2 * L)
    return r - 2 * L if r > L else r


def _smallest_prime_factor(n: int) -> int:
    p = 2
    while p * p <= n:
        if n % p == 0:
            return p
        p += 1
    return n


def chain(N: int) -> list[int]:
    """The levels of the chain: ``N`` divided by its smallest prime factor, the
    quotient again, and so on; the quotients above ``1``, ascending."""
    levels = []
    L = N
    while L > 1:
        levels.append(L)
        L //= _smallest_prime_factor(L)
    return levels[::-1]


def first_partner(N: int, k2: int) -> tuple[int, int]:
    """``(m2, m~2)`` in halves: the signed first partner ``m1`` of the irrep with
    exponent ``k2/2`` and its residue at the deciding level, the first level at
    which the residues of ``+k`` and ``-k`` differ (they always differ at ``N``);
    ``m1`` is the sign whose residue is positive, except that the ordinary
    ``E_{N/4}`` with ``N/4`` odd and greater than one takes ``-k``."""
    for L in chain(N):
        plus, minus = _residue(k2, L), _residue(-k2, L)
        if plus != minus:
            m2 = k2 if plus > 0 else -k2
            if L == N and k2 % 2 == 0 and N % 4 == 0 and N > 4 and (N // 4) % 2 == 1:
                m2 = -k2
            return m2, _residue(m2, L)
    raise ValueError(f"no deciding level for N = {N}, k2 = {k2}")


def spin_harmonic_images(
    F: CyclotomicField, N: int, k2: int, alpha_quarter_turns: int
) -> GeneratorImages:
    """``D(a), D(b)`` of the two-dimensional irrep with exponent ``k2/2`` of ``D_N``
    in the spin-harmonic basis: ``D(a) = diag(u^-m2, u^m2)`` with ``u = zeta_2N``
    and ``D(b) = exp(-i pi m~) [[0, y], [conj(y), 0]]`` with ``y = exp(-2 i m1
    alpha)``, ``alpha`` the azimuth of the anchored secondary pole in quarter
    turns (``0`` or ``1``)."""
    M = F.conductor
    if M % (2 * N) or M % 4:
        raise ValueError(f"Q(zeta_{M}) does not hold the spin-harmonic images of D_{N}")
    m2, mt2 = first_partner(N, k2)
    zero = F.zero
    Da = [
        [F.zeta(-m2 * M // (2 * N) % M), zero],
        [zero, F.zeta(m2 * M // (2 * N) % M)],
    ]
    t = F.zeta(-mt2 * M // 4 % M)
    y = F.zeta(-m2 * alpha_quarter_turns * M // 4 % M)
    Db = [[zero, F.mul(t, y)], [F.mul(t, F.conjugate(y)), zero]]
    return Da, Db


def secondary_azimuth(anchors: str, spec: GroupSpec) -> int:
    """The azimuth, in quarter turns, of the pole the reference assigns to the
    element its label convention anchors to ``b = (0, 1)``.

    Under the ``altmann-herzig`` anchoring the pole is on ``x`` (``C2'1`` of
    ``Dn``, ``Dnh``, ``Dnd``; ``sigma_v1`` of ``Cnv`` with ``4 | n``) except for
    ``Cnv`` with ``4`` not dividing ``n``, whose ``sigma_v1`` the book puts on
    ``y`` (the twist of ADR-0007 for ``n = 2 (mod 4)``, the book's own numbering
    of T n.1 for odd ``n``). The test is a rule in ``n``, so it extrapolates
    beyond the convention's range like the labels do.
    """
    conventions.label_convention(anchors)
    return 1 if spec.family == "Cnv" and spec.family_n % 4 != 0 else 0


def canonical_images(F: CyclotomicField, N: int, k2: int) -> GeneratorImages:
    """``D(a), D(b)`` of the canonical basis: ``diag(u^k2, u^-k2)`` on ``a`` and
    ``[[0, 1], [1, 0]]`` (ordinary) or ``[[0, i], [i, 0]]`` (projective) on ``b``."""
    zero = F.zero
    Da = [
        [F.root_of_unity(2 * N, k2), zero],
        [zero, F.root_of_unity(2 * N, -k2)],
    ]
    off = F.constant("i") if k2 % 2 else F.one
    return Da, [[zero, off], [off, zero]]


def axial_images(
    name: str, G: PointGroup, F: CyclotomicField, k2: int
) -> GeneratorImages:
    """``D(a), D(b)`` of the two-dimensional irrep with exponent ``k2/2`` of the
    dihedral rotation group of ``G`` under the convention's axial rule."""
    R = G.rotation_group
    N = R.axis_order
    assert N is not None and R.name.startswith("D")
    rules = rules_for(name, G)
    if rules.axial_rule == "spin-harmonic":
        assert rules.anchors is not None
        alpha = secondary_azimuth(rules.anchors, G.spec)
        # the ordinary irreps of Cnv with odd n carry the D_n matrices of C2'(k)
        if G.spec.family == "Cnv" and G.spec.family_n % 2 == 1 and k2 % 2 == 0:
            alpha = 0
        return spin_harmonic_images(F, N, k2, alpha)
    return canonical_images(F, N, k2)


# ---- the polyhedral generating set ------------------------------------------------------


def lifted_generators(
    name: str, cover: BinaryRotationGroup, R_name: str
) -> tuple[tuple[Any, ...], dict[str, tuple[Matrix, ...]]]:
    """The generating set of the convention's polyhedral data for ``R-dagger``,
    lifted by the canonical section (the barred lift for a ``~`` prefix), with
    the generator images per canonical irrep symbol."""
    data = conventions.polyhedral_data(name, cover.name)
    proper = point_group(parse_group(R_name))
    by_label = {
        label.text: r
        for (r, _), label in zip(proper.elements, canonical_labels(proper), strict=True)
    }
    lifted = []
    for label in data.generators:
        barred = label.startswith("~")
        r = by_label[label[1:] if barred else label]
        q = cover.section(r)
        lifted.append(cover.mul(cover.ebar, q) if barred else q)
    return tuple(lifted), data.irreps


# ---- the intertwiner --------------------------------------------------------------------


def conjugate_transpose(F: CyclotomicField, A: Matrix) -> Matrix:
    return [[F.conjugate(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def intertwiner(
    F: CyclotomicField, reference: list[Matrix], canonical: list[Matrix]
) -> Matrix | None:
    """A non-zero ``U`` with ``U D_can(g) = D_ref(g) U`` for every ``g``, by Schur
    averaging ``sum_g D_ref(g) X D_can(g)^dagger`` over the elementary matrices
    ``X`` (both realisations unitary and evaluated at the same section, so the
    factor system cancels); ``None`` when the two are inequivalent."""
    d = len(canonical[0])
    adjoints = [conjugate_transpose(F, D) for D in canonical]
    for i in range(d):
        for j in range(d):
            U = [[F.zero] * d for _ in range(d)]
            for D, A in zip(reference, adjoints, strict=True):
                # D X A with X = E_ij is the outer product of column i of D and row j of A
                for r in range(d):
                    if F.is_zero(D[r][i]):
                        continue
                    for c in range(d):
                        U[r][c] = F.add(U[r][c], F.mul(D[r][i], A[j][c]))
            if any(not F.is_zero(x) for row in U for x in row):
                return U
    return None


def unitary_scale(F: CyclotomicField, U: Matrix) -> Value | None:
    """The non-zero real ``c`` with ``U^dagger U = c 1``, or ``None``: an intertwiner
    between two unitary irreducible realisations is unitary up to such a scale
    (``U^dagger U`` commutes with the irrep, so it is a scalar, and it is positive
    definite); the scale need not be rational, nor a square in the field."""
    d = len(U)
    product = matrix_mul(F, conjugate_transpose(F, U), U)
    c = product[0][0]
    if F.is_zero(c) or not F.is_real(c):
        return None
    scaled = [[F.mul(c, x) for x in row] for row in matrix_identity(F, d)]
    return c if product == scaled else None
