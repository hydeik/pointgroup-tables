"""Prototype: comparison with the papers' printed values (wayfinder ticket #50).

The chain-adapted matrices of ``T > C3`` ([CF98a] eqs 37-38, as restated in
``docs/research/eigenfunction-method.md``), of ``O > C4`` ([FC99] eqs 32-33,
the entries whose extraction is certain) and the ``T > C3`` Clebsch-Gordan
table of [CF98c] Table III, transcribed by hand from the page text, compared
with the prototype's values in the three layers of docs/spec/clebsch-gordan.md
("Relation to the papers' conventions"): a label map from the eigenvalues of
the principal element (``m = -e``), a monomial intertwiner for the matrices,
the isotypic projectors (exact equality) and one phase per multiplicity-free
coupled vector. Nothing here is production code; the transcription is a
working copy, not the reference data of tickets #53 and #54.
"""

from __future__ import annotations

import math
from fractions import Fraction
from itertools import product
from typing import Any

from chain_prototype_basis import ChainBasis
from chain_prototype_cg import Pair

from pointgroup_tables.document import Model
from pointgroup_tables.field import CyclotomicField, Matrix, Value, cyclotomic_field
from pointgroup_tables.values import parse, pretty

Weight = Fraction

# ---- [CF98a]: T > C3 matrices of C2z (components by weight m, rho_m = exp(-2 pi i m / 3)) --

T_C3_C2Z: dict[str, tuple[list[Weight], list[list[str]]]] = {
    "F": (
        [Fraction(1), Fraction(0), Fraction(-1)],
        [["-1/3", "2/3", "2/3"], ["2/3", "-1/3", "2/3"], ["2/3", "2/3", "-1/3"]],
    ),
    "E1/2": (
        [Fraction(1, 2), Fraction(-1, 2)],
        [["-1/3*i*sqrt3", "1/3*i*sqrt2*sqrt3"], ["1/3*i*sqrt2*sqrt3", "1/3*i*sqrt3"]],
    ),
    "E3/2": (
        [Fraction(3, 2), Fraction(1, 2)],
        [["-1/3*i*sqrt3", "1/3*i*sqrt2*sqrt3"], ["1/3*i*sqrt2*sqrt3", "1/3*i*sqrt3"]],
    ),
    "E-3/2": (
        [Fraction(-3, 2), Fraction(-1, 2)],
        [["1/3*i*sqrt3", "-1/3*i*sqrt2*sqrt3"], ["-1/3*i*sqrt2*sqrt3", "-1/3*i*sqrt3"]],
    ),
}

# ---- [FC99]: O > C4 matrices (weights mu, rho_mu = exp(-2 pi i mu / 4)); C2y a coordinate half
# turn perpendicular to C4z, C2c a diagonal half turn. Entries marked uncertain in the survey are
# included and reported separately.

O_C4: dict[str, tuple[list[Weight], dict[str, list[list[str]]], bool]] = {
    "E": (
        [Fraction(0), Fraction(2)],
        {
            "C2y": [["1", "0"], ["0", "1"]],
            "C2c": [["-1/2", "1/2*sqrt3"], ["1/2*sqrt3", "1/2"]],
        },
        True,
    ),
    "F0": (
        [Fraction(1), Fraction(-1), Fraction(0)],
        {"C2y": [["0", "1", "0"], ["1", "0", "0"], ["0", "0", "-1"]]},
        True,
    ),
    "F2": (
        [Fraction(1), Fraction(-1), Fraction(2)],
        {"C2y": [["0", "1", "0"], ["1", "0", "0"], ["0", "0", "-1"]]},
        True,
    ),
    "E1/2": (
        [Fraction(1, 2), Fraction(-1, 2)],
        {"C2y": [["0", "-1"], ["1", "0"]]},
        True,
    ),
    "E3/2": (
        [Fraction(3, 2), Fraction(-3, 2)],
        {"C2y": [["0", "-1"], ["1", "0"]]},
        True,
    ),
    "G3/2": (
        [Fraction(-3, 2), Fraction(-1, 2), Fraction(1, 2), Fraction(3, 2)],
        {
            "C2y": [
                ["0", "0", "0", "1"],
                ["0", "0", "1", "0"],
                ["0", "-1", "0", "0"],
                ["-1", "0", "0", "0"],
            ],
            "C2c": [
                [
                    "-1/4*i*sqrt2",
                    "-1/4*i*sqrt2*sqrt3",
                    "1/4*i*sqrt2",
                    "-1/4*i*sqrt2*sqrt3",
                ],
                [
                    "-1/4*i*sqrt2*sqrt3",
                    "1/4*i*sqrt2",
                    "-1/4*i*sqrt2*sqrt3",
                    "-1/4*i*sqrt2",
                ],
                [
                    "1/4*i*sqrt2",
                    "-1/4*i*sqrt2*sqrt3",
                    "1/4*i*sqrt2",
                    "1/4*i*sqrt2*sqrt3",
                ],
                [
                    "-1/4*i*sqrt2*sqrt3",
                    "-1/4*i*sqrt2",
                    "1/4*i*sqrt2*sqrt3",
                    "-1/4*i*sqrt2",
                ],
            ],
        },
        False,
    ),
}
O_C4_UNCERTAIN = {"G3/2": ["C2c"]}

# ---- [CF98c] Table III: T > C3 Clebsch-Gordan coefficients ---------------------------------
# (alpha, beta): list of (gamma, copy, {m: [(coefficient, (m1, m2)), ...]})

H = Fraction(1, 2)
R3, R23, R6, R2 = "1/3*sqrt3", "1/3*sqrt2*sqrt3", "1/6*sqrt2*sqrt3", "1/2*sqrt2"


def _w(x: Any) -> Fraction:
    return Fraction(x)


T_C3_CG: dict[
    tuple[str, str],
    list[tuple[str, int, dict[Fraction, list[tuple[str, tuple[Fraction, Fraction]]]]]],
] = {
    ("F", "F"): [
        (
            "A0",
            0,
            {
                _w(0): [
                    (R3, (_w(0), _w(0))),
                    (R3, (_w(1), _w(-1))),
                    (R3, (_w(-1), _w(1))),
                ]
            },
        ),
        (
            "A1",
            0,
            {
                _w(1): [
                    (R3, (_w(0), _w(1))),
                    (R3, (_w(1), _w(0))),
                    (R3, (_w(-1), _w(-1))),
                ]
            },
        ),
        (
            "A-1",
            0,
            {
                _w(-1): [
                    (R3, (_w(0), _w(-1))),
                    (R3, (_w(-1), _w(0))),
                    (R3, (_w(1), _w(1))),
                ]
            },
        ),
        (
            "F",
            0,
            {
                _w(0): [
                    (R23, (_w(0), _w(0))),
                    ("-" + R6, (_w(1), _w(-1))),
                    ("-" + R6, (_w(-1), _w(1))),
                ],
                _w(1): [
                    (R23, (_w(-1), _w(-1))),
                    ("-" + R6, (_w(1), _w(0))),
                    ("-" + R6, (_w(0), _w(1))),
                ],
                _w(-1): [
                    (R23, (_w(1), _w(1))),
                    ("-" + R6, (_w(-1), _w(0))),
                    ("-" + R6, (_w(0), _w(-1))),
                ],
            },
        ),
        (
            "F",
            1,
            {
                _w(0): [(R2, (_w(1), _w(-1))), ("-" + R2, (_w(-1), _w(1)))],
                _w(1): [(R2, (_w(0), _w(1))), ("-" + R2, (_w(1), _w(0)))],
                _w(-1): [(R2, (_w(-1), _w(0))), ("-" + R2, (_w(0), _w(-1)))],
            },
        ),
    ],
    ("F", "E1/2"): [
        (
            "E1/2",
            0,
            {
                H: [(R3, (_w(0), H)), ("-" + R23, (_w(1), -H))],
                -H: [("-" + R23, (_w(-1), H)), ("-" + R3, (_w(0), -H))],
            },
        ),
        (
            "E3/2",
            0,
            {
                H: [(R23, (_w(0), H)), (R3, (_w(1), -H))],
                3 * H: [("-" + R3, (_w(1), H)), (R23, (_w(-1), -H))],
            },
        ),
        (
            "E-3/2",
            0,
            {
                -H: [(R23, (_w(0), -H)), ("-" + R3, (_w(-1), H))],
                -3 * H: [("-" + R23, (_w(1), H)), ("-" + R3, (_w(-1), -H))],
            },
        ),
    ],
    ("F", "E3/2"): [
        (
            "E1/2",
            0,
            {
                H: [(R23, (_w(0), H)), ("-" + R3, (_w(-1), 3 * H))],
                -H: [(R3, (_w(-1), H)), (R23, (_w(1), 3 * H))],
            },
        ),
        (
            "E3/2",
            0,
            {
                H: [(R3, (_w(0), H)), (R23, (_w(-1), 3 * H))],
                3 * H: [(R23, (_w(1), H)), ("-" + R3, (_w(0), 3 * H))],
            },
        ),
        (
            "E-3/2",
            0,
            {
                -H: [("-" + R23, (_w(-1), H)), (R3, (_w(1), 3 * H))],
                -3 * H: [(R23, (_w(0), 3 * H)), (R3, (_w(1), H))],
            },
        ),
    ],
    ("F", "E-3/2"): [
        (
            "E1/2",
            0,
            {
                H: [("-" + R3, (_w(1), -H)), ("-" + R23, (_w(-1), -3 * H))],
                -H: [(R23, (_w(0), -H)), ("-" + R3, (_w(1), -3 * H))],
            },
        ),
        (
            "E3/2",
            0,
            {
                H: [(R3, (_w(-1), -3 * H)), ("-" + R23, (_w(1), -H))],
                3 * H: [(R23, (_w(0), -3 * H)), (R3, (_w(-1), -H))],
            },
        ),
        (
            "E-3/2",
            0,
            {
                -H: [("-" + R23, (_w(1), -3 * H)), ("-" + R3, (_w(0), -H))],
                -3 * H: [(R3, (_w(0), -3 * H)), ("-" + R23, (_w(-1), -H))],
            },
        ),
    ],
    ("E1/2", "E1/2"): [
        ("A0", 0, {_w(0): [(R2, (H, -H)), ("-" + R2, (-H, H))]}),
        (
            "F",
            0,
            {
                _w(0): [(R2, (H, -H)), (R2, (-H, H))],
                _w(1): [("1", (H, H))],
                _w(-1): [("-1", (-H, -H))],
            },
        ),
    ],
    ("E1/2", "E3/2"): [
        ("A1", 0, {_w(1): [(R2, (H, H)), ("-" + R2, (-H, 3 * H))]}),
        (
            "F",
            0,
            {
                _w(0): [("-1", (-H, H))],
                _w(1): [(R2, (H, H)), (R2, (-H, 3 * H))],
                _w(-1): [("1", (H, 3 * H))],
            },
        ),
    ],
    ("E1/2", "E-3/2"): [
        ("A-1", 0, {_w(-1): [(R2, (-H, -H)), (R2, (H, -3 * H))]}),
        (
            "F",
            0,
            {
                _w(0): [("1", (H, -H))],
                _w(1): [("1", (-H, -3 * H))],
                _w(-1): [(R2, (-H, -H)), ("-" + R2, (H, -3 * H))],
            },
        ),
    ],
    ("E3/2", "E3/2"): [
        ("A-1", 0, {_w(-1): [(R2, (H, 3 * H)), ("-" + R2, (3 * H, H))]}),
        (
            "F",
            0,
            {
                _w(0): [("-1", (3 * H, 3 * H))],
                _w(1): [("1", (H, H))],
                _w(-1): [("-" + R2, (H, 3 * H)), ("-" + R2, (3 * H, H))],
            },
        ),
    ],
    ("E3/2", "E-3/2"): [
        ("A0", 0, {_w(0): [(R2, (H, -H)), (R2, (3 * H, -3 * H))]}),
        (
            "F",
            0,
            {
                _w(0): [(R2, (H, -H)), ("-" + R2, (3 * H, -3 * H))],
                _w(1): [("1", (3 * H, -H))],
                _w(-1): [("1", (H, -3 * H))],
            },
        ),
    ],
    ("E-3/2", "E-3/2"): [
        ("A1", 0, {_w(1): [(R2, (-H, -3 * H)), ("-" + R2, (-3 * H, -H))]}),
        (
            "F",
            0,
            {
                _w(0): [("-1", (-3 * H, -3 * H))],
                _w(1): [("-" + R2, (-3 * H, -H)), ("-" + R2, (-H, -3 * H))],
                _w(-1): [("1", (-H, -H))],
            },
        ),
    ],
}


# ---- label maps ---------------------------------------------------------------------------


def rho(F: CyclotomicField, m: Weight, n: int) -> Value:
    """``exp(-2 pi i m / n)`` as an element of ``F`` (``M`` a multiple of ``2n``)."""
    M = F.conductor
    k = -m * M // n
    assert k.denominator == 1
    return F.zeta(int(k))


def component_map(
    F: CyclotomicField, D: Matrix, weights: list[Weight], n: int, projective: bool
) -> tuple[list[int], int] | None:
    """For a diagonal ``D`` (the principal element in the chain basis), the
    component index of every paper weight ``m`` with eigenvalue ``rho_m``; the
    sign is the lift (``-1`` when the barred preimage matches a projective irrep)."""
    d = len(D)
    for sign in (1,) if not projective else (1, -1):
        out: list[int] = []
        for m in weights:
            target = rho(F, m, n)
            if sign < 0:
                target = F.neg(target)
            hits = [k for k in range(d) if D[k][k] == target]
            if len(hits) != 1 or hits[0] in out:
                break
            out.append(hits[0])
        else:
            return out, sign
    return None


def monomial_intertwiner(
    F: CyclotomicField,
    A: Matrix,
    B: Matrix,
    perm: list[int],
    phases: list[Value | None],
) -> bool:
    """Extend ``phases`` (``lambda_r``) so that ``A[r][c] = lambda_r conj(lambda_c) B[perm[r]][perm[c]]``;
    ``False`` on an inconsistency."""
    d = len(A)
    changed = True
    while changed:
        changed = False
        for r in range(d):
            for c in range(d):
                a, b = A[r][c], B[perm[r]][perm[c]]
                if F.is_zero(a) != F.is_zero(b):
                    return False
                if F.is_zero(a):
                    continue
                ratio = F.divide(a, b)
                pr, pc = phases[r], phases[c]
                if pr is not None and pc is not None:
                    if F.mul(pr, F.conjugate(pc)) != ratio:
                        return False
                elif pc is not None:
                    phases[r] = F.mul(ratio, pc)
                    changed = True
                elif pr is not None:
                    phases[c] = F.conjugate(F.divide(ratio, pr))
                    changed = True
    return True


Intertwiners = dict[str, tuple[list[int], list[Value]]]


def compare_matrices(
    basis: ChainBasis,
    paper: dict[str, tuple[list[Weight], dict[str, list[list[str]]], bool]],
    name_map: dict[str, str],
    principal: int,
    element_candidates: dict[str, list[int]],
    n: int,
    uncertain: dict[str, list[str]] | None = None,
    found: Intertwiners | None = None,
) -> list[str]:
    """One line per paper irrep: the component map, the candidate elements that
    match up to a monomial intertwiner, and the phases found. Matrices listed in
    ``uncertain`` are compared after the certain ones, with the phases those fixed.
    ``found`` collects ``(perm, phases)`` per paper irrep for the coupling layer."""
    F = cyclotomic_field(math.lcm(basis.conductor, 24))  # the papers' sqrt2, sqrt3
    G = basis.chain.G
    lines = []
    uncertain = uncertain or {}
    for paper_name, (weights, matrices, certain) in paper.items():
        doubtful = uncertain.get(paper_name, [])
        ours = next(
            ci for ci in basis.irreps if ci.irrep.symbol.text == name_map[paper_name]
        )
        lifted = [
            [[basis.field.lift(x, F) for x in row] for row in D] for D in ours.matrices
        ]
        cm = component_map(F, lifted[principal], weights, n, ours.irrep.projective)
        if cm is None:
            lines.append(
                f"{paper_name} -> {name_map[paper_name]}: no component map from the eigenvalues of the principal element"
            )
            continue
        perm, sign = cm
        names = [e for e in matrices if e not in doubtful]
        matches = []
        for choice in product(*(element_candidates[e] for e in names)):
            phases: list[Value | None] = [F.one] + [None] * (len(weights) - 1)
            ok = True
            for e, g in zip(names, choice, strict=True):
                A = [[parse(F, x) for x in row] for row in matrices[e]]
                B = lifted[g]
                # a global sign on D(g) of a half turn is the lift: try both
                if not monomial_intertwiner(F, A, B, perm, phases):
                    Bn = [[F.neg(x) for x in row] for row in B]
                    phases = [F.one] + [None] * (len(weights) - 1)
                    if not monomial_intertwiner(F, A, Bn, perm, phases):
                        ok = False
                        break
            if ok:
                matches.append((choice, phases))
        labels = {e: [G.labels[g].text for g in element_candidates[e]] for e in names}
        if matches:
            choice, phases = matches[0]
            found_phases = [pretty(F, p) if p is not None else "?" for p in phases]
            if found is not None:
                found[paper_name] = (
                    perm,
                    [p if p is not None else F.one for p in phases],
                )
            for e in doubtful:
                A = [[parse(F, x) for x in row] for row in matrices[e]]
                hits = []
                for g in element_candidates[e]:
                    for B in (
                        lifted[g],
                        [[F.neg(x) for x in row] for row in lifted[g]],
                    ):
                        trial = list(phases)
                        if monomial_intertwiner(F, A, B, perm, trial):
                            hits.append(G.labels[g].text)
                            break
                lines.append(
                    f"{paper_name}: the uncertain {e} matches {hits or 'NO element'} with the phases of the certain matrices"
                )
            lines.append(
                f"{paper_name} -> {name_map[paper_name]}: components {perm} (lift sign {sign}); "
                f"matches with {', '.join(f'{e} = {G.labels[g].text}' for e, g in zip(names, choice, strict=True))} "
                f"up to the diagonal unitary {found_phases}{'' if certain else ' [uncertain extraction]'}"
                + (
                    f"; {len(matches)} element choices match"
                    if len(matches) > 1
                    else ""
                )
            )
        else:
            lines.append(
                f"{paper_name} -> {name_map[paper_name]}: NO monomial intertwiner over {labels}{'' if certain else ' [uncertain extraction]'}"
            )
    return lines


# ---- the T > C3 Clebsch-Gordan comparison -------------------------------------------------


def compare_t_cg(
    model: Model,
    basis: ChainBasis,
    pairs: dict[tuple[int, int], Pair],
    Fcg: CyclotomicField,
    principal: int,
    intertwiners: Intertwiners,
) -> list[str]:
    """Isotypic-projector equality after the basis layer (the component map and
    the diagonal phases of the matrix comparison), and the phase of every
    multiplicity-free coupled vector, for every triple of [CF98c] Table III."""
    F = Fcg
    names = {
        "A0": "A",
        "A1": "2E",
        "A-1": "1E",
        "F": "T",
        "E1/2": "E1/2",
        "E3/2": "2E3/2",
        "E-3/2": "1E3/2",
    }
    weights_of = {
        "A0": [Fraction(0)],
        "A1": [Fraction(1)],
        "A-1": [Fraction(-1)],
        "F": [Fraction(1), Fraction(0), Fraction(-1)],
        "E1/2": [Fraction(1, 2), Fraction(-1, 2)],
        "E3/2": [Fraction(3, 2), Fraction(1, 2)],
        "E-3/2": [Fraction(-3, 2), Fraction(-1, 2)],
    }
    index_of = {irrep.symbol.text: irrep.index for irrep in model.irreps}
    lines = []

    def cmap(paper_name: str) -> tuple[dict[Weight, int], dict[Weight, Value]]:
        ci = basis.irreps[index_of[names[paper_name]]]
        if paper_name in intertwiners:
            perm, phases = intertwiners[paper_name]
            return dict(zip(weights_of[paper_name], perm, strict=True)), dict(
                zip(weights_of[paper_name], phases, strict=True)
            )
        D = [[basis.field.lift(x, F) for x in row] for row in ci.matrices[principal]]
        cm = component_map(F, D, weights_of[paper_name], 3, ci.irrep.projective)
        assert cm is not None, paper_name
        return dict(zip(weights_of[paper_name], cm[0], strict=True)), {
            m: F.one for m in weights_of[paper_name]
        }

    for (pa, pb), triples in T_C3_CG.items():
        a, b = index_of[names[pa]], index_of[names[pb]]
        pair = pairs[(min(a, b), max(a, b))]
        swapped = a > b
        (ma, la), (mb, lb) = cmap(pa), cmap(pb)
        da, db = model.irreps[a].dimension, model.irreps[b].dimension
        dim = da * db

        def row_index(
            m1: Weight, m2: Weight, ma=ma, mb=mb, swapped=swapped, da=da, db=db
        ) -> int:
            i, j = ma[m1], mb[m2]
            if swapped:
                return j * da + i
            return i * db + j

        by_gamma: dict[str, list[list[Value]]] = {}
        for gamma, _copy, comps in triples:
            for terms in comps.values():
                v = [F.zero] * dim
                for coeff, (m1, m2) in terms:
                    # |m>_paper = conj(lambda_m) |perm m>_ours (D_paper = Lambda D_ours Lambda^dagger)
                    v[row_index(m1, m2)] = F.mul(
                        parse(F, coeff), F.conjugate(F.mul(la[m1], lb[m2]))
                    )
                by_gamma.setdefault(gamma, []).append(v)
        for gamma, vectors in by_gamma.items():
            g = index_of[names[gamma]]
            ours = [c for copy in pair.copies if copy.gamma == g for c in copy.columns]

            def projector(vs: list[list[Value]], dim=dim) -> Matrix:
                return [
                    [
                        F.total(F.mul(v[r], F.conjugate(v[c])) for v in vs)
                        for c in range(dim)
                    ]
                    for r in range(dim)
                ]

            equal = projector(vectors) == projector(ours)
            phase = ""
            if len(vectors) == model.irreps[g].dimension and equal:
                # one phase per multiplicity-free triple after the column layer: for the
                # paper's component m (our component perm m with phase lambda_m),
                # ours = phase * conj(lambda_m) * paper's
                mg, lg = cmap(gamma)
                phases = set()
                for _gamma, _copy, comps in [t for t in triples if t[0] == gamma]:
                    for m, terms in comps.items():
                        v = [F.zero] * dim
                        for coeff, (m1, m2) in terms:
                            v[row_index(m1, m2)] = F.mul(
                                parse(F, coeff), F.conjugate(F.mul(la[m1], lb[m2]))
                            )
                        w = ours[mg[m]]
                        inner = F.total(
                            F.mul(F.conjugate(x), y) for x, y in zip(v, w, strict=True)
                        )
                        phases.add(pretty(F, F.mul(inner, F.conjugate(lg[m]))))
                phase = f"; phase {sorted(phases)}"
            lines.append(
                f"{pa} x {pb} > {gamma}: projector {'equal' if equal else 'DIFFERENT'} ({len(vectors)} vectors){phase}"
            )
    return lines


__all__ = [
    "O_C4",
    "O_C4_UNCERTAIN",
    "T_C3_C2Z",
    "T_C3_CG",
    "Intertwiners",
    "compare_matrices",
    "compare_t_cg",
    "cyclotomic_field",
]
