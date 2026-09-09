# Basis conventions

Decided in wayfinder ticket
[Basis conventions for irrep matrices: canonical form and reference conventions](https://github.com/hydeik/pointgroup-tables/issues/10).
Builds on [Abstract groups](./abstract-groups.md) (rotation groups and their
generators), [Double groups](./double-groups.md) (the quaternion models of
`2T`, `2O`, `2I`, the generator pairs `q₁, q₂` and the canonical section),
[Irrep labels](./irrep-labels.md) (the principal element `c` of every family
and the canonical Mulliken symbols), [Values](./values.md) (the document
field) and the research note `docs/research/irrep-algorithms.md`.
Rationale in ADR-0004.

Irrep matrices are defined only up to unitary equivalence, so a matrix table
is meaningless without naming its basis. This chapter fixes the **canonical
basis** the tool uses by default and the framework by which a **reference
basis convention** (Altmann & Herzig first) replaces it.

## Two layers

Every irrep matrix table is written in exactly one basis:

- the **canonical basis**, convention-free, fixed by the rules below for the
  axial families and by shipped generator data for the polyhedral families;
  the default, named `canonical`;
- or a **reference basis convention**, a shipped data file that realises
  every irrep as a reference work writes it, selected with
  `--basis-convention <name>`.

Both layers only concern irreps of dimension `d ≥ 2`; a one-dimensional
irrep has no basis freedom and its matrix is its character.

Ordinary and projective irreps are treated alike: a projective irrep is a
representation of `G†`, tabulated on `G` through the canonical section with
the factor system `ω` ([Double groups](./double-groups.md)); the matrix of
the barred element is `−D(g)`.

## The spinor map

One fixed homomorphism `ψ` from the unit quaternions to SU(2) underlies
every two-dimensional projective irrep:

```
ψ(w, x, y, z) = [[ w + i·z,  −y + i·x ],
                 [ y + i·x,   w − i·z ]]
```

so that a rotation by `θ` about the `z` axis (quaternion
`(cos θ/2, 0, 0, sin θ/2)`) is `diag(e^{+iθ/2}, e^{−iθ/2})`, the **positive
exponent first**, the binary rotation about `x` (quaternion `(0, 1, 0, 0)`)
is `[[0, i], [i, 0]]`, and the binary rotation about `y` is
`[[0, −1], [1, 0]]`. `ψ` is a homomorphism (`ψ(i)ψ(j) = ψ(k)`), unitary, and
exact over any field containing the coordinates and `i`.

## Canonical basis: the principle

The canonical basis of every irrep of dimension `d ≥ 2` is **adapted to the
principal element** `c` of the family ([Irrep labels](./irrep-labels.md#named-elements):
`(1, 0)` of `R` for the axial families; `C4^1(1)` for `O`, `Oh`, `Td`;
`C3^1(1)` for `T`, `Th`; `C5^1(1)` for `I`, `Ih`), lifted by the canonical
section for projective irreps:

1. `D(c)` (or `D(s(c))`) is **diagonal**;
2. its eigenvalues `ζ^{e_1}, ..., ζ^{e_d}` are listed with the exponents
   `e_1 > e_2 > ... > e_d` **descending**, the exponent taken in
   `(−N/2, N/2]` for an ordinary irrep (`ζ = ζ_N`) and in `(−N, N]` for a
   projective irrep (`ζ = ζ_2N`), `N` the order of `c`.

For every irrep of every point group these eigenvalues are pairwise distinct
(checked for all polyhedral irreps, e.g. `H` of `I` on `C5`:
`1, ζ_5^{±1}, ζ_5^{±2}`; `I5/2` of `I†` on `s(C5)`: `ζ_20^{±1, ±3, ±5}`), so
rules 1 and 2 fix the basis up to one phase per basis vector. The first basis
vector of a two-dimensional axial `E_k` therefore transforms as the
component-`1` member `1E_k` of the corresponding pair, and subduction to the
cyclic subgroup `⟨c⟩` is read off the diagonal.

The remaining phases are fixed differently in the two kinds of family.

## Canonical basis: axial families

Let `R = C_N` or `D_N`, `a = (1, 0) = c` and `b = (0, 1)` the secondary
element, `w = ζ_N`, `u = ζ_2N`. The generator images are, for every `n`:

| Irrep | `D(a)` | `D(b)` (dihedral only) |
| ----- | ------ | ---------------------- |
| ordinary `E_k`, `1 ≤ k < N/2` | `diag(w^k, w^{−k})` | `[[0, 1], [1, 0]]` |
| projective `E_{k/2}`, `k` odd, `1 ≤ k < N` | `diag(u^k, u^{−k})` on `s(a)` | `[[0, i], [i, 0]]` on `s(b)` |

These are `ψ` applied to the standard embedding of `Dic_N` in SU(2) (`a` a
rotation about `z`, `b` the binary rotation about `x`,
[Double groups](./double-groups.md#canonical-section)), so no phase choice
remains. They satisfy `b² = Ē`, `(ab)² = Ē`, `a^N = Ē` (projective) and the
corresponding ordinary relations exactly. Every other element's matrix is
the product along its canonical word `a^k b^s` with the factor system of the
canonical section; for the projective irreps `D(s(a^k b^s)) = ±D(s(a))^k
D(s(b))^s`, the sign being `ω`-determined ([Double groups](./double-groups.md#factor-system)).

Type II groups (`Cnh` even `n`, `S2n` odd `n`, `Dnh` even `n`, `Dnd` odd
`n`) use the same matrices with `D(i·r) = ±D(r)` (`g`: `+`, `u`: `−`). Type
III groups (`Cnv`, `Cnh`/`Dnh` odd `n`, `S2n` even `n`, `Dnd` even `n`) are
`R` itself with parities attached, so their matrices are those of `R`.

The canonical **conductor** of an axial group is fixed by `R` alone:

| `R` | `M` |
| --- | --- |
| `C_N`, `N` odd | `N` |
| `C_N`, `N` even | `2N` |
| `D_N` | `lcm(2N, 4)` |

(`u = ζ_2N` lies in `Q(ζ_N)` for odd `N`; `i` is needed only when a
secondary element exists.)

## Canonical basis: polyhedral families

For `T`, `O`, `I` and their double groups, rules 1 and 2 leave one phase per
basis vector of every irrep of dimension `≥ 2`. These phases are **not**
fixed by a formula: they are chosen once so that every matrix lies in the
smallest field the irrep allows, and frozen as **shipped data**, the
canonical basis-convention file described below. The data gives, for each
of `2T`, `2O`, `2I`, the images of the generator pair `q₁, q₂` of
[Double groups](./double-groups.md#generators-and-the-projection-onto-t-o-i)
in every irrep of `R†` of dimension `≥ 2`, and the conductor of `R†`. The
construction and verification of this data is the prototype ticket
[Canonical polyhedral bases: generator matrices for every irrep of 2T, 2O, 2I](https://github.com/hydeik/pointgroup-tables/issues/16).

The reference constructions the prototype starts from (all exact):

- `E1/2` of `R†`: `U† ψ(q) U` with `U` the unitary diagonalising `ψ(s(c))`;
- the three-dimensional irrep with `χ(c) = 1 + 2 cos(2π/N)` (`T` of `T`,
  `T1` of `O`, `T1` of `I`): the rotation matrix of the quaternion,
  conjugated into the eigenbasis of `c`;
- `T2` of `O`: `A2 ⊗ T1` (sign of the permutation); `T2` of `I`: the
  Galois conjugate `√5 ↦ −√5` of `T1`;
- everything else: irreducible blocks of tensor products of the above.

Every group whose rotation group is `R` shares the data of `R†`: type II
groups (`Th`, `Oh`, `Ih`) use `D(i·r) = ±D(r)`, and `Td` (`R = O`, twisted)
uses `D((r, χ(r))) = D_O(r)`; the canonical Mulliken symbols of `Td` and of
the `g` block of `Oh` coincide with those of `O`, so the same irrep key
serves all three.

## Reference basis conventions

A reference basis convention is one TOML file per reference work, shipped
with the package, in a namespace separate from the label conventions (the
same reference gets the same `name` in both). It declares its coverage: a
rule for the axial families and per-`R†` data for the polyhedral families.
A group the convention does not cover is an error with exit code 4
([Command line](./cli.md#exit-codes-and-errors)); the tool never falls back
to the canonical basis silently.

```toml
name = "altmann-herzig"
source = "S. L. Altmann and P. Herzig, Point-Group Theory Tables, 2nd ed."

[axial]
rule = "transform"                    # "canonical" or "transform"
ordinary = [["0", "1"], ["1", "0"]]   # U for every two-dimensional ordinary irrep, every n
projective = [["1", "0"], ["0", "i"]] # U for every two-dimensional projective irrep, every n

[polyhedral.2O]
conductor = 8
generators = ["C4^1(1)", "C3^1(1)"]   # canonical labels of elements of R; "~" prefix selects the barred lift

[polyhedral.2O.irreps.E]              # canonical Mulliken symbol of the irrep of R (= of the proper group)
matrices = [                          # one matrix per generator, rows then columns
  [["1", "0"], ["0", "-1"]],
  [["-1/2", "(-1 + i*sqrt2)/2"], ["(1 + i*sqrt2)/2", "-1/2"]],   # illustrative shape only
]

[polyhedral.2O.irreps."E1/2"]
matrices = [ ... ]
```

- **Axial rule.** `canonical` means the canonical matrices. `transform`
  gives two constant unitary `2 × 2` matrices `U`, one for the ordinary and
  one for the projective two-dimensional irreps, applied for every `n` and
  `k` as `D_ref = U† D_can U`. Basis-order reversal, complex conjugation of
  the `E_k` block and phase changes are all of this form. A rule depending
  on `n` or `k` is out of the version 1 vocabulary.
- **Polyhedral data**, keyed by `R†` (`2T`, `2O`, `2I`), shared by every
  group with that rotation group exactly as for the canonical data. A
  convention may give **any generating set** of `R`, named by canonical
  labels; the tool checks that it generates and lifts each named element by
  the canonical section (`~` selects the barred lift). Irreps are keyed by
  their canonical Mulliken symbol; one-dimensional irreps are omitted.
- **Conductor** per `R†`, declared; every entry must lie in `Q(ζ_M)`.

The canonical basis is itself such a file (`name = "canonical"`, no
`source`, `axial.rule = "canonical"`, polyhedral data from the prototype),
so one loader serves both layers.

### Value grammar

Matrix entries are strings in a closed grammar, parsed exactly into the
declared field:

```
expr   := term (("+" | "-") term)*
term   := factor (("*" | "/") factor)*
factor := "-"? atom ("^" integer)?
atom   := integer | "i" | "zeta" | "sqrt2" | "sqrt3" | "sqrt5" | "phi" | "(" expr ")"
```

`zeta` is `ζ_M` for the declared conductor `M`; the other symbols are
embedded by fixed formulas, never by symbolic conversion:

| Symbol | Embedding | Requires |
| ------ | --------- | -------- |
| `i` | `ζ_M^{M/4}` | `4 │ M` |
| `sqrt2` | `ζ_8 + ζ_8^{−1}` with `ζ_8 = ζ_M^{M/8}` | `8 │ M` |
| `sqrt3` | `ζ_12 + ζ_12^{−1}` | `12 │ M` |
| `sqrt5` | `ζ_5 + ζ_5^4 − ζ_5^2 − ζ_5^3` | `5 │ M` |
| `phi` | `(1 + sqrt5) / 2` | `5 │ M` |

A symbol whose requirement the declared conductor does not meet is a data
error. No whitespace significance; integers are unbounded.

## Metadata

Every table document records `basis_convention = { name, source }`
(`{ name = "canonical" }` for the default) and the `field` fixed by the
group and the convention ([Values](./values.md#the-field-of-a-document)):
for a reference convention the conductor is the declared one for polyhedral
groups and, for axial groups, the canonical conductor enlarged to contain
the entries of `U`.

## Verification

Exactly, for every group and every basis convention:

- **Homomorphism by closure.** All matrices are generated by breadth-first
  closure from the generator images along the same words as the group
  elements; an element reached by two words must receive the same matrix
  (for projective irreps, the same matrix up to the sign `ω` prescribes).
  This is equivalent to `D(g) D(h) = ω(g, h) D(gh)` for all pairs, which the
  test suite also checks exhaustively for `|G| ≤ 48`.
- **Unitarity** of every generator image (`D† D = 1` with the field's
  conjugation `ζ ↦ ζ^{−1}`); closure preserves it.
- **Identity of the irrep.** The trace of every generated matrix equals
  the character of the irrep named by the key, computed independently by
  the Burnside class-matrix method; character orthogonality and
  `Σ d² = |G|` hold over the full set.
- **Canonical basis only:** `D(c)` is diagonal with the prescribed
  eigenvalue order (rules 1 and 2).
- **Reference conventions:** the intertwiner `U` between the reference and
  the canonical realisation of each irrep exists and is unique up to a
  phase (computed exactly by Schur averaging); `U` is never part of the
  input and is kept only for the verification chapter.
- **Field.** Every entry lies in the declared `Q(ζ_M)`; minimality of `M`
  is a test-suite property, not a runtime check.

## Version 1 scope

The canonical convention ships in full once the prototype ticket freezes
the polyhedral data. The Altmann & Herzig convention is filled by the task
[Transcribe Altmann & Herzig irrep matrices for the polyhedral groups](https://github.com/hydeik/pointgroup-tables/issues/17)
and the decision
[Altmann & Herzig basis-convention data](https://github.com/hydeik/pointgroup-tables/issues/18).
Bradley & Cracknell and Bilbao basis conventions are not yet specified.
