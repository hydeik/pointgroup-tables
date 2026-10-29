# Clebsch–Gordan coefficients

Decided in wayfinder ticket
[Clebsch–Gordan coefficients: definition, coupling multiplicity and the phase rule](https://github.com/hydeik/pointgroup-tables/issues/46),
on the facts of
[Survey of the eigenfunction method in the Chen–Fan–Ping papers](https://github.com/hydeik/pointgroup-tables/issues/39)
and
[Cyclotomic closure of chain-adapted bases and Clebsch–Gordan coefficients](https://github.com/hydeik/pointgroup-tables/issues/40);
the trade-offs are recorded in
[ADR-0012](../adr/0012-clebsch-gordan-convention.md). Builds on
[Basis conventions](./basis-conventions.md) (the irrep matrices the
coefficients are written in), [Double groups](./double-groups.md) (the
factor system), [Irrep labels](./irrep-labels.md) (canonical irrep order)
and [Values](./values.md) (the field of a document). This chapter fixes
**what a Clebsch–Gordan coefficient of this tool is, how the coupling
multiplicity is resolved and how every phase is fixed**, for every pair of
irreps (ordinary and projective) of every group, in whatever basis the
document names. The `cg` document and subcommand are fixed in the
[Document](#document) section, decided in
[The `cg` document and subcommand](https://github.com/hydeik/pointgroup-tables/issues/47);
the verification ledger rows by ticket
[#48](https://github.com/hydeik/pointgroup-tables/issues/48). Nothing here
changes a version-1 document.

## Definition

Let `α`, `β` be two irreps of the group `G` of a document, of dimensions
`d_α`, `d_β`, with matrices `D_α(g)`, `D_β(g)` in the basis the document
names (canonical, a reference basis convention, or a chain-adapted basis),
and `V_α ⊗ V_β` their product space with the **product basis**
`|α i⟩|β j⟩`, `0 ≤ i < d_α`, `0 ≤ j < d_β`, in **row-major product
order**: `i` outer, `j` inner, so that the product index is `i d_β + j`.

The **coupled vectors** of the pair are

```
|γ k; r⟩ = Σ_{i,j} C^{αβγ r}_{i j k} |α i⟩|β j⟩,     0 ≤ k < d_γ,  0 ≤ r < n_γ,
```

one orthonormal set per irrep `γ` of `G` and per **coupled copy** `r`, and
the **Clebsch–Gordan coefficients** `C^{αβγ r}_{i j k}` are their
coordinates. Collected over every `γ`, `r` and `k`, the coefficients form
the **coupling matrix** `U` of the pair, a `d_α d_β × d_α d_β` matrix whose
rows are the product basis in row-major order and whose columns are the
coupled vectors, ordered by `γ` in the canonical irrep order of `G`
([Irrep labels](./irrep-labels.md#canonical-irrep-order), the global index
that survives `--irreps`), then by the copy `r`, then by the component `k`
in the component order of the document's basis. The coupling matrix is
unitary and reduces the product representation:

```
D_α(g) ⊗ D_β(g) = U ( ⊕_{γ, r} D_γ(g) ) U†      for every g ∈ G.
```

This is the direction of the papers (the coupled state expanded in the
product states): `C^{λ a, λ' a'}_{λ'' a'' τ}` of Chen and Ping is
`C^{αβγ r}_{i j k}` with `(λ, λ', λ'') = (α, β, γ)`, `(a, a', a'') = (i, j, k)`
and `τ = r`.

### The product factor system

The product representation `D_α ⊗ D_β` has the factor system `ω_α ω_β`
([Double groups](./double-groups.md#factor-system)): with the one spin
factor system `ω` of the group, `ω² = 1`, so

| `α` | `β` | `α ⊗ β` | constituents `γ` |
| --- | --- | ------- | ---------------- |
| ordinary | ordinary | ordinary | ordinary irreps |
| projective | projective | ordinary | ordinary irreps |
| ordinary | projective | projective | projective irreps |

The number of projective irreps among `α`, `β`, `γ` is always even, so for
the barred element `D(ḡ) = ±D(g)` the signs cancel in every sum below: every
formula is evaluated over the elements of `G` with the matrices of the
canonical section, exactly as the `matrices` table prints them, and gives
the same result as the sum over `G†`.

### Direct-product decomposition

The **coupling multiplicity** of `γ` in `α ⊗ β` is computed from the
characters,

```
n_γ = (1/|G|) Σ_{g ∈ G} χ_α(g) χ_β(g) χ_γ(g)*,
```

a non-negative integer with `Σ_γ n_γ d_γ = d_α d_β`; only the `γ` of the
kind the table above names can occur. The coupling matrix holds `n_γ d_γ`
columns for `γ`, and none when `n_γ = 0`.

## The algorithm

The coefficients are computed by projection in the product space with the
irrep matrices of the document's basis, so the same code serves every
basis; nothing is read from a table. For every `γ` with `n_γ > 0` let

```
P^γ_{k l} = (d_γ / |G|) Σ_{g ∈ G} D_γ(g)_{k l}^*  (D_α(g) ⊗ D_β(g)),
```

the matrix units of the isotypic component of `γ` (`P^γ_{k k}` is a
Hermitian projector of rank `n_γ`, `P^γ_{k l} P^γ_{l m} = P^γ_{k m}`, the
Wigner projection operators of the eigenfunction method). Every value lies
in the document's field.

1. **Designated component.** `k₀ = 0`, the first component of `γ` in the
   component order of the document's basis.
2. **Trial product state.** For the copy `r`, let `Q` be the projector of
   the part of `P^γ_{k₀ k₀}` not yet exhausted (the [multiplicity rule](#the-multiplicity-rule)
   says which: `P^γ_{k₀ k₀}` itself, or its symmetric or antisymmetric part,
   minus the projectors `|w_s⟩⟨w_s| / N_s` of the copies `s < r` already
   taken in that part). The **trial product state** of the copy is the
   first product basis vector `e_{j₀}` in row-major product order with
   `Q_{j₀ j₀} ≠ 0`. The unnormalised coupled vector of the component `k₀` is

   ```
   w = Q e_{j₀},      N = ⟨w|w⟩ = Q_{j₀ j₀} = w_{j₀},
   ```

   the squared norm `N` a totally positive element of the field, equal to
   the designated coordinate of `w` (`Q` is a Hermitian idempotent). Since
   `Q` is positive semidefinite, `Q_{jj} = 0` forces `Q_{j j₀} = 0`: the
   coordinates of `w` on the product vectors `j < j₀` vanish, so the trial
   state is also the first product state with a non-zero coefficient in the
   designated component.
3. **The other components** follow from the transfer operators, with no
   further choice:

   ```
   w_k = P^γ_{k k₀} w      (k ≠ k₀),      ⟨w_k|w_k⟩ = N.
   ```

4. **Normalisation and phase.** Every component is divided by the
   **canonical norm factor** `μ` of `N` (`μ μ̄ = N`, the primitive of
   [Group chains](./group-chains.md) defined with the exact arithmetic of
   ticket [#44](https://github.com/hydeik/pointgroup-tables/issues/44)):

   ```
   |γ k; r⟩ = w_k / μ,      so that   C^{αβγ r}_{j₀ k₀} = N / μ = μ̄.
   ```

   The **designated coefficient** `C_{j₀ k₀}` is therefore the conjugate
   of the canonical norm factor of its own squared modulus; whenever the
   real square root `√N` lies in the document's field, `μ = √N` and the
   designated coefficient is **real and positive**. With the
   [conductor rule](#the-field-of-a-cg-document) this is the case for every
   coupled copy of every axial group and of the `T`- and `O`-type groups,
   and for every copy of the `I`-type groups except the four (fourteen for
   `Ih`) listed below.

The coupling matrix of a pair is thus determined by the irrep matrices of
the three irreps, the component orders, the canonical irrep order and the
canonical norm factor alone, **in the canonical basis**. In a chain-adapted
basis or a reference basis convention the coefficients are the **unitary
transform** of the canonical-basis ones (decided in ticket
[#44](https://github.com/hydeik/pointgroup-tables/issues/44)): with
`U_α`, `U_β`, `U_γ` the matrices whose columns are the basis vectors of
the three irreps in canonical components (`D'(g) = U† D(g) U`),

```
C'^{αβγ r}_{i' j' k'} = Σ_{i j k} conj(U_α)_{i i'} · conj(U_β)_{j j'} · C^{αβγ r}_{i j k} · (U_γ)_{k k'},
```

the same coupled states written in the new components. The coupled copies,
the multiplicity rule and the exchange signs are therefore independent of
the basis; the coefficients change (they are basis-dependent,
[Glossary](../../GLOSSARY.md)), and the designated coefficient of a copy is
real and positive in the canonical basis only. No squared norm beyond those
of the canonical basis arises, so the canonical norm factors a `cg`
document needs are those of the table in
[Exact computation](./group-chains.md#the-canonical-norm-factor).

## The multiplicity rule

When `n_γ > 1` the coupled copies of `γ` are fixed as follows.

- **`α = β`: symmetric copies first, then antisymmetric.** The swap
  `S |α i⟩|α j⟩ = |α j⟩|α i⟩` commutes with `D_α(g) ⊗ D_α(g)` for every
  `g` (also for a projective `α`), so `P_± = (1 ± S)/2` commute with every
  `P^γ_{k l}` and the isotypic component of `γ` splits into its symmetric
  part `P_+ P^γ_{k₀ k₀}` (multiplicity `n⁺_γ`) and its antisymmetric part
  `P_− P^γ_{k₀ k₀}` (multiplicity `n⁻_γ`, `n⁺_γ + n⁻_γ = n_γ`), with

  ```
  n^±_γ = (1/|G|) Σ_{g ∈ G} ½ ( χ_α(g)² ± tr(D_α(g)²) ) χ_γ(g)*
  ```

  (`tr(D_α(g)²)` rather than `χ_α(g²)`, so that the formula holds with the
  factor system). The copies `r = 0 … n⁺_γ − 1` are taken in the symmetric
  part, the copies `r = n⁺_γ … n_γ − 1` in the antisymmetric part; inside a
  part, by step 2 above with the Gram–Schmidt complement, in trial-state
  order. Every copy is an eigenvector of `S`; its eigenvalue is the
  **exchange sign** `ε_r ∈ {+1, −1}` of the copy (`+1` for the symmetric
  copies). A copy with `n_γ = 1` lies in one part and carries its sign too.
- **`α ≠ β`: trial-state order.** The copies are taken by step 2 with the
  Gram–Schmidt complement `Q_r = P^γ_{k₀ k₀} − Σ_{s < r} |w_s⟩⟨w_s| / N_s`,
  i.e. the copy `r` is projected from the first product state the previous
  copies have not exhausted. The copies are orthonormal.

The symmetric/antisymmetric split alone does not resolve every
multiplicity: `G3/2 ⊗ G3/2 ⊃ T1` of `O` has both copies in the symmetric
part, and `H ⊗ H ⊃ H`, `I5/2 ⊗ I5/2 ⊃ T1`, `I5/2 ⊗ I5/2 ⊃ T2` (symmetric
part) and `I5/2 ⊗ I5/2 ⊃ H` (antisymmetric part) of `I` keep a multiplicity
of two inside one part ([verification data](#verification-data)); the
trial-state order settles those. Unlike the papers, whose copies are
labelled by hand-picked trial states and are not orthogonal, every copy of
this tool is orthonormal and carries a definite exchange sign when
`α = β`.

## Pair order and the exchange relation

The algorithm runs on the pairs `(α, β)` with `α < β` or `α = β` in the
canonical irrep order (for a mixed pair the ordinary irrep is `α`). The
coefficients of the reversed pair are **defined** by exchanging the product
indices:

```
C^{βαγ r}_{j i k} = C^{αβγ r}_{i j k}             (α < β),
C^{ααγ r}_{j i k} = ε_r C^{ααγ r}_{i j k}         (α = β),
```

so the exchange relation holds by construction for every pair, with the
sign `+1` for distinct factors and the exchange sign of the copy for equal
ones. The coupling matrix of `(β, α)` is `U` with its rows permuted from the
`(α, β)` row-major order to the `(β, α)` one (and, for `α = β`, its columns
scaled by `ε_r`): it is again unitary and reduces `D_β ⊗ D_α`. The reversed
pair is **not** recomputed by the algorithm with the roles swapped; that
would fix the trial state in the other product order and could change the
copies and their phases. A `cg` document prints the pairs `α ≤ β` only
([Document](#the-document)); a consumer obtains the reversed pair from
these two formulas.

The exchange relation and unitarity are the symmetry relations the tool
promises. No relation under complex conjugation is promised: the conjugate
pair `(ᾱ, β̄)` and the conjugate constituent `γ̄` are coupled by the same
rule, and because the canonical basis reverses the component order of a
complex-conjugate pair of irreps (`1E3/2`, `2E3/2` of `T`), the conjugate of
a coupling matrix is a coupling matrix of the conjugate triple only up to
the intertwiners of the three pairs and a phase per copy; the verification
checks the subspaces, not the coefficients.

## The field of a `cg` document

The squared norms `N` of the algorithm are rational on every polyhedral
group and `1` or `1/2` on every axial group (proven for every `n` in the
cyclotomic-closure note, observed for the polyhedral groups). The conductor
of a `cg` document is the conductor the same options would give to a
`matrices` document (`M_basis`: the version-1 `M` of the canonical basis, the
convention's conductor under `--basis-convention`, the chain conductor
`M_chain` of ticket #44 under `--chain`), enlarged so that the square roots
of the multiplicity-free norms are in the field:

| Family | `M_cg` | Canonical basis |
| ------ | ------ | --------------- |
| axial (`Cn`, `Cnh`, `S2n`, `Cnv`, `Dn`, `Dnd`, `Dnh`) | `lcm(M_basis, 8)` | adds `√2`: `C3v` 24, `D4h` 8, `D5h` 40, `D6h` 24 |
| polyhedral (`T`, `Th`, `Td`, `O`, `Oh`, `I`, `Ih`) | `lcm(M_basis, 24)` | adds `√2`, `√3`: `T`, `Th` 24; `Td`, `O`, `Oh` 24; `I`, `Ih` 120 |

These are the conductors of the `altmann-herzig` polyhedral data. Under
`--chain` the chain conductor is fixed per chain class
([Group chains](./group-chains.md#the-conductor-of-a-chain-document)), so
the `cg` conductor is `lcm(M_chain, 24)` for a polyhedral group (`24`, or
`120` for `I` and `Ih`) and `lcm(M_G, 8)` for an axial group. With this field every coupled copy of every
axial, `T`-type and `O`-type group has `√N` in the field (its designated
coefficient is real and positive), and so has every copy of `I` and `Ih`
except the copies with the prime `7` in `N`: `T1 ⊗ I5/2 ⊃ I5/2` (both
copies, `N = 7/10` and `1/7`) and `I5/2 ⊗ I5/2 ⊃ T1` (both symmetric
copies, `N = 7/20` and `2/7`) of `I`, and their fourteen `g`/`u` versions in
`Ih`. Those designated coefficients are the conjugates of the canonical
norm factors of `Q(ζ_120)` (`7 = |2 + i*sqrt3|²`: `1/7 ↦ (2 + i*sqrt3)/7`,
`2/7 ↦ (2*sqrt2 + i*sqrt6)/7`, `7/10 ↦ (2*sqrt10 + i*sqrt30)/10`,
`7/20 ↦ (2*sqrt5 + i*sqrt15)/10`, shipped entries of
[Exact computation](./group-chains.md#the-canonical-norm-factor)), of
modulus `√N`. The real-positive closure of every norm would be
`Q(ζ_840)` (degree 192, six times the value size) for the sake of these
four copies; ADR-0012 declines it. The `field` metadata of version 1
describes the result without new keys.

If the canonical norm factor of a squared norm does not exist in the
document's field, the tool fails with the internal error of
[Exact computation](./group-chains.md#errors) (exit `1`, never an
approximation); the condition is unreachable on shipped data, every
squared norm of every pair of the seven polyhedral groups being a square
of the field or one of the four entries above.

## Document

Decided in
[The `cg` document and subcommand](https://github.com/hydeik/pointgroup-tables/issues/47).
The machine form is the `cg` table of
[Output schema](./output-schema.md#the-cg-table) and
[`schema.json`](./schema.json); the command line is in
[Command line](./cli.md#table-subcommands).

### The `cg` subcommand

`pointgroup-tables cg <group> [options]` is a **table subcommand**: exactly
one group, named and normalised as for every other table subcommand
([Command line](./cli.md#the-group-argument)), the same exit codes, and
every common option (`--format`, `--output`, `--values`,
`--label-convention`, `--basis-convention`, `--chain`, `--verify`). One
invocation emits one document holding the Clebsch–Gordan coefficients of
**every pair of irreps** of the group that `--irreps` selects, in the basis
the document names; no option selects a single pair, as no option of
version 1 selects a single irrep.

- `--irreps` selects pairs by their **factors**: `ordinary` the pairs of two
  ordinary irreps (whose constituents are ordinary), `projective` the pairs
  with **at least one** projective factor (ordinary ⊗ projective, whose
  constituents are projective, and projective ⊗ projective, whose
  constituents are ordinary), `all` (the default) every pair. The three
  values partition the pairs as they partition the rows of a character
  table: `ordinary` and `projective` are disjoint and their union is `all`.
- `--double-group` is a usage error (exit `2`) on `cg`, as on `characters`
  and `matrices`: every sum of the algorithm runs over the elements of `G`
  with the matrices of the canonical section, and no index of the table
  refers to an element of `G†`. The `elements` table of a `cg` document is
  always that of `G` and `double` is always `false`.
- `--verify` runs the self-consistency checks `CG-1` … `CG-7` and the
  [document checks](#document-checks) on the document before emitting it;
  `CG-8` is a reference comparison and never runs under `--verify`.
- `cg` is **not part of `all`**, for two reasons: `all` is frozen as the
  union of the five version-1 tables, so that a version-1 `all` document is
  reproduced byte for byte; and a `cg` document declares a different field
  from the other tables of the same group and options
  ([The field of a `cg` document](#the-field-of-a-cg-document)), so adding
  the table would change the `field` record and the exact values of every
  other table of the document. Opt-in inclusion (`all --with cg`) is fog of
  the map, decided only if a consumer asks.

### The document

A `cg` document is a table document of [Output schema](./output-schema.md):
the metadata keys of every table document, the `elements` table of `G`,
under `--chain` the `subgroup_chain` table, and the `cg` table:

```
tables = ["elements", "cg"]                      (default)
tables = ["elements", "subgroup_chain", "cg"]    (with --chain)
```

The metadata is that of a `matrices` document with the same options except
for `field`: `irreps` is present (the document holds an irrep-indexed
table), `double` is `false`, `field` declares the `M_cg` of the table
above, and `basis_convention` and `chain` name the basis exactly as they do
for `matrices`. The document does **not** embed the `matrices` table: the
basis is determined by `basis_convention`, `chain` and `generator.version`
(documents are deterministic), and a consumer who needs the irrep matrices
produces the `matrices` document with the same options and matches the
irreps by their global `index`. Embedding them would double the largest
documents and would need a filtering rule of its own, because the
constituents of the pairs `--irreps projective` selects are of both kinds.
Instead the `cg` table carries the irrep records its indices refer to.

```json
"cg": {
  "irreps": [ <irrep record>, ... ],
  "pairs":  [ <pair record>, ... ]
}
```

- `irreps`: the **full** irrep list of `G` in canonical irrep order, never
  filtered by `--irreps` (the constituents of a projective ⊗ projective
  pair are ordinary, so both kinds are referred to under every selection):
  the same records as `characters.irreps` would hold unfiltered, with
  `reference_symbol` under a label convention, and under `--chain` each
  record followed by the `components` array of
  [Output schema](./output-schema.md#the-matrices-table), the chain labels
  of the components the coupling matrices are indexed by.
- `pairs`: one pair record per **unordered pair** `(α, β)` with `α ≤ β` in
  the global irrep index, both factors in the `--irreps` selection, in
  lexicographic order of `(α, β)`. The pair `(β, α)` is never printed; its
  coefficients follow from the
  [exchange relation](#pair-order-and-the-exchange-relation) by exchanging
  the product indices: `U` with its rows permuted from the `(α, β)`
  row-major order to the `(β, α)` one, and for `α = β` the columns of every
  copy scaled by its exchange sign.

**Pair record**, keys in this order:

| Key | Content |
| --- | ------- |
| `irreps` | `[α, β]`, the global irrep indices of the two factors, `α ≤ β` |
| `dimension` | `d_α d_β`, the order of the coupling matrix |
| `constituents` | one constituent record per irrep `γ` with `n_γ > 0`, in canonical irrep order |
| `columns` | one column record per column of the coupling matrix, `dimension` entries |
| `exact`, `numeric`, `pretty` | the coupling matrix `U` in each requested form, a `dimension × dimension` matrix in the layout of [Values](./values.md#layout-in-tables) (`[row][column]`) |

**Constituent record**: `irrep` (the global index of `γ`), `multiplicity`
(`n_γ`) and `columns` (`[first, count]`: the first column of `γ` and the
number of its columns, `count = n_γ d_γ`). The constituents are the
direct-product decomposition of the pair; an irrep with `n_γ = 0` has no
record.

**Column record**: `index` (the column's position), `irrep` (`γ`), `copy`
(`r`, `0 ≤ r < n_γ`), `component` (`k`, `0 ≤ k < d_γ`, in the component
order of the document's basis) and, for `α = β` only, `exchange`, the
exchange sign `ε_r` of the copy (`1` or `-1`; the key is absent for
`α ≠ β`, never `null`). The columns are ordered by `γ` in canonical irrep
order, then `copy`, then `component`, as the [definition](#definition)
fixes.

**Rows** carry no records: row `i d_β + j` of the coupling matrix is the
product state `|α i⟩|β j⟩` in row-major product order (`i` outer), the
components `i`, `j` in the component order of the document's basis. The
entry `U[i d_β + j][c]` of the column `c = (γ, r, k)` is `C^{αβγ r}_{i j k}`.

The exchange sign, the trial product state and the squared norm of a copy
are all readable from the matrix (the sign by applying the swap, the trial
state as the first non-zero row of the column `k = 0`, the squared norm as
the squared modulus of that entry); only the exchange sign is printed,
because it is the one of the three the tool promises as an invariant
(`CG-4`) and the one a consumer coupling identical factors needs first.

Example: the pair `(E, E)` of `C3v` (`E` has index `2`, `d = 2`), in which
`E ⊗ E = A1 ⊕ A2 ⊕ E`:

```json
{ "irreps": [2, 2],
  "dimension": 4,
  "constituents": [ { "irrep": 0, "multiplicity": 1, "columns": [0, 1] },
                    { "irrep": 1, "multiplicity": 1, "columns": [1, 1] },
                    { "irrep": 2, "multiplicity": 1, "columns": [2, 2] } ],
  "columns": [ { "index": 0, "irrep": 0, "copy": 0, "component": 0, "exchange": 1 },
               { "index": 1, "irrep": 1, "copy": 0, "component": 0, "exchange": -1 },
               { "index": 2, "irrep": 2, "copy": 0, "component": 0, "exchange": 1 },
               { "index": 3, "irrep": 2, "copy": 0, "component": 1, "exchange": 1 } ],
  "exact":   [ <4 rows of 4 exact values> ],
  "numeric": [ <4 rows of 4 numeric values> ] }
```

### Index spaces

A `cg` document uses three index spaces, all of
[Output schema](./output-schema.md#index-spaces): the **global irrep index**
(`irreps` of the pair record, `irrep` of the constituent and column
records, into `cg.irreps`), the **component index** of an irrep in the
basis the document names (`component` of a column record and the factors
`i`, `j` of a row; under `--chain` the chain-adapted component order, whose
chain labels are `cg.irreps[*].components`), and the **product index**
`i d_β + j` of a row. No element or class index occurs in the table; the
`elements` table is embedded as in every document.

### What the options change

| Option | Effect on the `cg` document |
| ------ | --------------------------- |
| `--irreps` | the set of pairs ([above](#the-cg-subcommand)); never the `irreps` list of the table |
| `--values` | the forms present in every pair record; `pretty` uses the value grammar with the fallback of [Values](./values.md#pretty-form-optional) |
| `--label-convention` | `reference_symbol` in the records of `cg.irreps` and the `label_convention` metadata; nothing else |
| `--basis-convention` | the coefficients are the unitary transform into the convention's basis ([The algorithm](#the-algorithm)); `basis_convention` metadata with `coverage`; `field` is `lcm(M_conv, 8)` / `lcm(M_conv, 24)`; the convention must cover the group (exit `4` otherwise, as for `matrices`) |
| `--chain` | the `chain` metadata, the `subgroup_chain` table second in `tables`, `components` in every record of `cg.irreps`; the coefficients are the unitary transform into the chain-adapted basis; `field` is `lcm(M_chain, 24)` / `lcm(M_G, 8)` |
| `--format`, `--output`, `--verify` | as for every document |

The `cg` table is a function of the group, the basis the document names
and the `--irreps` selection: two documents of one group with the same
basis and selection hold the same `cg` table form by form (a form present
in both is byte-identical), whatever the other options.

### Sizes

Estimated from the entry counts of the
[verification data](#verification-data) (unordered pairs of all irreps,
`--values both`, JSON): an exact value of degree `φ(M)` takes about
`7 φ(M)` bytes on one line and, when it does not fit in 100 columns (degree
`32` of `Q(ζ_120)`), one coefficient pair per line
([JSON formatting](./output-schema.md#json-formatting)).

| Document | Entries | Field | Estimated size | Prototype (#50) |
| -------- | ------: | ----- | -------------: | --------------: |
| `cg C3v` | 90 | `Q(ζ_24)` | 20 KB | 24 KB |
| `cg D6h` | 1236 | `Q(ζ_24)` | 150 KB | 244 KB |
| `cg Oh` | 5076 | `Q(ζ_24)` | 0.5 MB | 0.7 MB |
| `cg Ih` | 31428 | `Q(ζ_120)` | 15–20 MB | 24 MB (1.1 million lines) |

`cg Ih` is the largest document of the tool, above `matrices Ih` (5–6 MB),
because its field has degree `32` and its exact values break across lines;
`--values numeric` or `--values exact` roughly halves it and
`--irreps ordinary` reduces it to a quarter (`cg Oh --irreps ordinary`:
0.2 MB in the prototype). The prototype column was measured on the sample
documents of ticket [#50](https://github.com/hydeik/pointgroup-tables/issues/50)
(`prototypes/chain_cg_prototype.py`, `--values both`, JSON; the column
records account for about a third of the lines); the implementation
measures the actual sizes and replaces the estimates in the sizes table of
[Output schema](./output-schema.md#sizes).

### Golden documents

The golden documents of `cg` follow the pattern of
[Verification](./verification.md#test-organisation):
`tests/golden/<group>/cg.<format>` for `C3v`, `D6h`, `Oh`, `Ih` with
default options (`GOLD-1`) and `tests/golden/<group>/cg.chain-default.<format>`
under `--chain default` (`GOLD-2`), in the three formats. The `Ih` files
are the largest golden documents of the repository (about 20 MB each, three
formats, two option sets); the precedent is `matrices Ih`.

## Relation to the papers' conventions

The printed Clebsch–Gordan tables of Chen, Fan and Ping (`T ⊃ C3`,
`T ⊃ D2 ⊃ C2`, `T ⊃ C2`, `O ⊃ C4`, `O ⊃ T ⊃ C3` and the algebraic `Dn ⊃ Cn`
formulas; survey of ticket #39, `docs/research/eigenfunction-method.md`)
relate to this tool's coefficients in three layers, which the reference
comparison `REF-11` of [Verification](./verification.md#comparison-of-coupling-matrices)
applies in this order:

1. **Basis.** The papers' components are the chain labels of their chain,
   with the quantum number `m = −e` of our exponents, so the components of
   every irrep labelled by the cyclic quantum number alone are **reversed**
   relative to the canonical order (exceptions as the survey lists them:
   `F_σ` and `G3/2` of `O ⊃ C4`, `G3/2` of `O ⊃ T ⊃ …`, `(x, y, z)` of
   `T ⊃ D2`); their `E3/2` of `T` is our `2E3/2`, their `A_{±1}` our `2E`,
   `1E`. Their matrices are ours in the chain-adapted basis of the same
   chain up to the permutation and the diagonal phase matrix of ticket #43;
   the coupling matrices transform with `U_α ⊗ U_β` on the rows and `U_γ`
   on the columns.
2. **Phase per copy.** The papers fix the coefficient of a hand-chosen trial
   product state real and positive; this tool fixes the first non-zero
   coefficient of the designated component (in product order) to `μ̄`. For a
   multiplicity-free triple the two coupled vectors span the same line, so
   they agree up to one phase, which is `±1` wherever both are real.
3. **Multiplicity mixing.** For `n_γ > 1` the papers' copies are labelled by
   their trial states and are not orthogonal in general (the printed tables
   orthogonalise them); the tool's copies are orthonormal with symmetric
   copies first. The two sets of copies span the same isotypic subspace and
   are related by an `n_γ × n_γ` invertible (unitary, when both are
   orthonormal) matrix.

The invariant of all three layers is the **isotypic projector**
`Π_γ = Σ_{r,k} |γ k; r⟩⟨γ k; r|` of each `γ` in `α ⊗ β`, written in the
product basis: after layer 1 it must agree exactly with the projector built
from the paper's columns. The papers' exchange rule
`(−1)^{Γ₁+Γ₂+Γ+r}` for `T` and the `(−1)^r` of their `F ⊗ F ⊃ F` copies
coincide with the exchange sign of this chapter on `T ⊃ D2 ⊃ C2`, whose two
copies are the symmetric and antisymmetric ones; for the other groups the
exchange sign is this tool's own.

## Verification

### Properties

An implementation must pass, for every pair `(α, β)` of every group of the
`default` and `large` sets, ordinary and projective, in every basis the
document can name:

- **CG-1 unitarity**: `U† U = 1` exactly;
- **CG-2 reduction**: `(D_α(g) ⊗ D_β(g)) U = U (⊕_{γ,r} D_γ(g))` for every
  `g ∈ G`, with the matrices of the document's basis and the canonical
  section (this is the check that proves the factor system of the product);
- **CG-3 multiplicities**: the number of columns of `γ` equals `n_γ d_γ`
  with `n_γ` from the characters, and `Σ_γ n_γ d_γ = d_α d_β`; for `α = β`
  the split `n⁺_γ + n⁻_γ = n_γ` holds and the copies are counted per part;
- **CG-4 exchange**: the coupling matrix of `(β, α)` is the row permutation
  of that of `(α, β)` (scaled by `ε_r` for `α = β`), and for `α = β` every
  copy is an eigenvector of `S` with eigenvalue `ε_r`;
- **CG-5 designated coefficient**: in every column `(γ, r, k₀ = 0)` the first
  non-zero entry in product order is the conjugate of the canonical norm
  factor of its squared modulus, real and positive whenever `√N` lies in
  the field; for the symmetric copies of `α = β` it is the first non-zero
  entry of the symmetrised column (canonical basis; in a chain-adapted or
  reference basis the coupling matrix equals the unitary transform of the
  canonical one, entry by entry);
- **CG-6 field**: every entry lies in the declared field, and the declared
  conductor follows the table above;
- **CG-7 conjugation**: `n_{γ̄}(ᾱ, β̄) = n_γ(α, β)`, and the isotypic
  projector of `γ̄` in `ᾱ ⊗ β̄` equals the conjugate of that of `γ` in
  `α ⊗ β` transported by the intertwiners `D_{ᾱ} = W_α D_α^* W_α†`,
  `D_{β̄} = W_β D_β^* W_β†` (the intertwiners of `MX-6`);
- **CG-8 reference**: the comparison with the papers' printed tables in the
  three layers above: exact equality of the isotypic projectors after
  layer 1, and, for multiplicity-free triples, equality of the coupled
  vectors up to one phase (a labelled comparison), with `known_discrepancy`
  records where a table differs;
- **mutation**: one coefficient's sign flipped is caught by `CG-2`; one
  column multiplied by a phase is caught by `CG-5`; two copies of a
  multiplicity-two triple swapped are caught by `CG-4` or `CG-5`.

These are the rows `CG-1` to `CG-7` and `CG-9` of the
[check ledger](./verification.md#check-ledger) (decided in ticket
[#48](https://github.com/hydeik/pointgroup-tables/issues/48)): `CG-4`
there is the exchange-eigenvalue half, checked on the document, and `CG-9`
the reversed-pair half, checked on the implementation's objects since a
document prints `α ≤ β` only; `CG-8` is the reference row `REF-11`; the
mutations are `MUT-11` to `MUT-13`. `CG-2` runs on the canonical
generating set of `G`: both sides are projective representations with the
factor system `ω_α ω_β`, so equality on the generators is equality on `G`.
All run on `default` and `large`, `slow` for `|G†| > 96`; `--verify` runs
`CG-1` to `CG-7` on every `cg` document, recomputing the irrep matrices of
the document's basis.

### Document checks

An implementation must check, on every `cg` document of the `seed` set,
for every `--irreps` value and every `--values` form, with and without
`--chain default` (rows `OS-9`, `OS-10`, `OS-11`, `GOLD-1` and `GOLD-2` of the
[check ledger](./verification.md#check-ledger)):

- the document validates against `schema.json`; `tables` is
  `["elements", "cg"]` or `["elements", "subgroup_chain", "cg"]`; `double`
  is `false`; `irreps` is present; `field.conductor` is the `M_cg` of the
  [field table](#the-field-of-a-cg-document);
- `cg.irreps` is the full irrep list of `G`, equal to the unfiltered
  `characters.irreps` of the same options, with `components` added under
  `--chain` and absent otherwise;
- `pairs` holds exactly the pairs `α ≤ β` of the selection, in
  lexicographic order; the pair sets of `ordinary` and `projective` are
  disjoint and their union is that of `all`;
- in every pair record `dimension = d_α d_β`, every form is a
  `dimension × dimension` matrix, `columns` has `dimension` records with
  `index` equal to position, ordered by `irrep` (canonical order), `copy`,
  `component`; `constituents` lists exactly the `γ` with `n_γ > 0` from the
  characters, with `columns = [first, n_γ d_γ]` contiguous and in order;
  `exchange` is present in every column record exactly when `α = β`;
- `cg --double-group` gives `2`; the `all` document of the same options is
  unchanged (`tables`, `field` and bytes) by the existence of `cg`;
- the [golden documents](#golden-documents) equal the checked-in files byte
  for byte.

### Verification data

Exact facts of the canonical basis, computed with the research script of
ticket #40 extended by the rules of this chapter (`cg_multiplicity.py`,
2026-10-28), to be reproduced by the implementation:

- **Triples with `n_γ > 1`** (`α ≤ β`): `T`: `T ⊗ T ⊃ T` (2). `O`:
  `T1 ⊗ G3/2 ⊃ G3/2` (2), `T2 ⊗ G3/2 ⊃ G3/2` (2), `G3/2 ⊗ G3/2 ⊃ T1` (2),
  `G3/2 ⊗ G3/2 ⊃ T2` (2). `I`: `T1 ⊗ I5/2 ⊃ I5/2` (2), `T2 ⊗ I5/2 ⊃ I5/2`
  (2), `G ⊗ H ⊃ H` (2), `G ⊗ G3/2 ⊃ I5/2` (2), `G ⊗ I5/2 ⊃ G3/2` (2),
  `G ⊗ I5/2 ⊃ I5/2` (2), `H ⊗ G3/2 ⊃ I5/2` (2), `H ⊗ I5/2 ⊃ G3/2` (2),
  `H ⊗ I5/2 ⊃ I5/2` (3), `G3/2 ⊗ I5/2 ⊃ G` (2), `G3/2 ⊗ I5/2 ⊃ H` (2),
  `H ⊗ H ⊃ G` (2), `H ⊗ H ⊃ H` (2), `I5/2 ⊗ I5/2 ⊃ T1` (2),
  `I5/2 ⊗ I5/2 ⊃ T2` (2), `I5/2 ⊗ I5/2 ⊃ G` (2), `I5/2 ⊗ I5/2 ⊃ H` (3).
  None for `C3v`, `D6h` (every axial group, proven). `Oh` and `Ih` have the
  `g`/`u` versions (14 and 62 triples).
- **Symmetric / antisymmetric split** (`n⁺_γ + n⁻_γ`): `T ⊗ T ⊃ T` of `T`:
  `1 + 1`; `O`: `G3/2 ⊗ G3/2 ⊃ T1`: `2 + 0`, `⊃ T2`: `1 + 1`; `I`:
  `H ⊗ H ⊃ G`: `1 + 1`, `H ⊗ H ⊃ H`: `2 + 0`, `I5/2 ⊗ I5/2 ⊃ T1`, `⊃ T2`:
  `2 + 0`, `⊃ G`: `1 + 1`, `⊃ H`: `1 + 2`. Exchange signs of the
  multiplicity-free `α = β` triples of `T`: antisymmetric are
  `E1/2 ⊗ E1/2 ⊃ A`, `1E3/2 ⊗ 1E3/2 ⊃ 2E`, `2E3/2 ⊗ 2E3/2 ⊃ 1E`, every other
  symmetric; of `O`: antisymmetric are `E ⊗ E ⊃ A2`, `T1 ⊗ T1 ⊃ T1`,
  `T2 ⊗ T2 ⊃ T1`, `E1/2 ⊗ E1/2 ⊃ A1`, `E5/2 ⊗ E5/2 ⊃ A1`,
  `G3/2 ⊗ G3/2 ⊃ A1`, `G3/2 ⊗ G3/2 ⊃ E`, every other symmetric
  (`T1 ⊗ T1`: `A1`, `E`, `T2` symmetric, `T1` antisymmetric, the textbook
  result).
- **Squared norms `N`** of the copies: `T`: `{1, 2/3, 1/2, 1/3, 1/6}`
  (`T ⊗ T ⊃ T`: `1/6` symmetric, `1/2` antisymmetric); `O`, `Oh`:
  `{1, 3/4, 2/3, 1/2, 3/8, 1/3, 1/4, 1/6}` (`T1 ⊗ G3/2 ⊃ G3/2`: `1/2, 1/3`;
  `T2 ⊗ G3/2 ⊃ G3/2`: `2/3, 3/4`; `G3/2 ⊗ G3/2 ⊃ T1`: `3/8, 1/2`;
  `G3/2 ⊗ G3/2 ⊃ T2`: `1/4, 1/4`); `I`: 28 rationals
  `{1, 9/10, 4/5, 3/4, 7/10, 2/3, 3/5, 9/16, 8/15, 1/2, 2/5, 9/25, 7/20, 1/3, 3/10, 2/7, 4/15, 1/4, 1/5, 1/6, 1/7, 2/15, 1/9, 1/10, 1/12, 1/15, 1/20, 1/30}`,
  `Ih` these and `2/9`, `5/12`, `9/20`; every axial group `{1, 1/2}`. A
  canonical norm factor exists in the version-1 field for every one of
  them; the real-positive closure is `Q(ζ_24)` for `T`, `O`, `Oh`, `C3v`,
  `D6h` and `Q(ζ_120)` for `I`, `Ih` without the four `7`-copies,
  `Q(ζ_840)` with them.
- **Sizes** of the coupling matrices over the unordered pairs `α ≤ β` of
  all irreps (entries, non-zero entries): `C3v` 21 pairs, 90, 44; `D6h`
  171, 1236, 468; `Oh` 136, 5076, 1310; `Ih` 171, 31428, 5798; ordinary
  pairs only: `C3v` 6, 27, 13; `D6h` 78, 324, 160; `Oh` 55, 1332, 378;
  `Ih` 55, 8244, 1558.
