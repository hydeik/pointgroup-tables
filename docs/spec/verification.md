# Verification

Decided in wayfinder ticket
[Verification strategy: self-consistency tests and reference-table comparison](https://github.com/hydeik/pointgroup-tables/issues/12).
Every other chapter ends with a "Verification" section listing the
properties its tables must have. This chapter is the index of all of them:
it fixes where each check runs, on which groups, against which oracle, and
what it takes for the specification to count as implementable. A check is
defined once, in the chapter that owns the table; this chapter assigns it an
identifier and a place in the test suite. Rationale for redistributing the
reference data in ADR-0006.

## Two kinds of proof

An implementation proves its tables in two independent ways:

1. **Self-consistency checks** prove that the tables are *mathematically
   correct*: the Cayley table is a group, the factor system is a normalised
   cocycle with values `±1`, the characters are orthogonal, the matrices
   are unitary homomorphisms with the right traces, and so on. They need no
   external data and hold for every group the tool can construct, including
   an axis order no book tabulates. Every self-consistency check is exact
   (ADR-0002); there is no tolerance anywhere.
2. **Reference comparison** proves that the tables *agree with the
   literature*: that the names, the class memberships, the sign conventions
   of the factor system and the bases of the matrices are the ones the
   reference works print. Self-consistency cannot detect a wrongly named
   irrep or a barred element assigned to the wrong class; only a
   transcribed reference can.

Self-consistency is the proof of correctness; reference comparison is the
proof of convention.

## Where checks run

The **test suite** (`pytest`, under `tests/`) is the home of every check.
It is the proof that ships with a version: the implementation is correct
for every group the suite covers.

In addition, every table subcommand accepts the flag **`--verify`**
([Command line](./cli.md#options)). It runs the self-consistency checks
that apply to the document about to be emitted, for the group and options
actually requested, and only emits the document when all of them pass. It
exists because the axis order `n` is unbounded: the test suite covers a
finite set of groups, and a user who generates `D37h` has `--verify` to
prove that document. `--verify` is off by default; its cost is roughly
`O(|G†|² · Σ d³)` and reaches minutes for large `n`.

`--verify` runs:

- the **base checks** for every table subcommand: the group checks
  (`AG-*`), the class checks (`CL-*`), and, whenever the double group is
  involved (`--double-group`, or a document that contains projective
  irreps), the factor-system checks (`DG-*`);
- the **table checks** for each table in the document: `CH-*` for
  `characters`, `MX-*` for `matrices`, `VL-*` for every table that carries
  values. `IL-*` runs with `characters` and `matrices`;
- the **chain checks** for a `--chain` document: `OS-7` and `OS-8` on the
  `subgroup_chain` table and the `components`, and `CB-2` (every step-1
  block equals the sub-chain realisation). `MX-5` does **not** run under
  `--chain`: a chain-adapted basis does not diagonalise the principal
  element, and `OS-8` (block-diagonal on every step, diagonal on the
  terminal) takes its place; `MX-1` to `MX-4`, `MX-7` and `MX-8` run
  unchanged on the chain-adapted matrices;
- the **coupling checks** for a `cg` document: `CG-1` to `CG-7` on every
  pair of the document (`CG-2` on the canonical generating set, see the
  ledger), and the `MX-*` checks on the irrep matrices the coefficients were
  computed from, which the document does not print
  ([Clebsch–Gordan coefficients](./clebsch-gordan.md#document)) and
  `--verify` recomputes for the document's basis.

It never runs a reference comparison (reference data is test data, not
shipped), nor the output-schema checks other than the chain checks `OS-7`
and `OS-8` above, nor the CLI, golden or mutation checks. Two
version-2 checks need an object the document does not serialise and are
therefore test-suite checks on the implementation's own objects, never run
by `--verify`: the basis matrices of a chain-adapted basis (`CB-1`, `CB-3`:
the phase rule is a convention, not a correctness property, and the
document carries chain labels, not canonical coordinates) and the coupling
matrix of a reversed pair (`CG-9`: a `cg` document prints the pairs
`α ≤ β` only).

On failure `--verify` does not stop at the first check: it runs all of them
and reports every failure, one line each, then exits with code `1` and
writes no document ([Command line](./cli.md#exit-codes-and-errors)):

```
pointgroup-tables: error: verification of 'D3' failed (2 checks)
verify: CH-1: rows 1 and 2 are not orthogonal: <D(1), D(2)> = 6
verify: MX-3: D(C3^1) of irrep 2 is not unitary
```

The identifier after `verify:` is the check identifier of the ledger below.
Every checker is a pure function of the computed document model (the
tables as they will be serialised), never of internal state, so the same
function serves `--verify`, the test suite and the mutation tests.

## Group sets

The test suite is parametrised over named sets of groups. Aliases
(`Cs = C1h`, `D1h = C2v`, `S2 = Ci`, ...) are resolved by the group
argument grammar and counted once.

| Set | Contents | Marker |
| --- | -------- | ------ |
| `seed` | the 43 groups of the seed table in [Double groups](./double-groups.md#verification-data): the 32 crystallographic groups, `I`, `Ih`, `C5`, `C5v`, `D5h`, `D5d`, `C7v`, `D7h`, `D8d`, `C8h` | none |
| `sweep` | every axial family `Cn`, `Cnv`, `Cnh`, `S2n`, `Dn`, `Dnd`, `Dnh` for `n = 1 ... 12` | none |
| `large` | every axial family for `n ∈ {13, 16, 17, 24, 30}` | `slow` |
| `default` | `seed ∪ sweep` | none |

`sweep` exhausts the case distinctions of the labelling and basis rules
(odd and even `n`, `n` divisible by 4 or not, the degenerate `n = 1, 2`
families); `large` proves that the conductor formulas and the `E_k`
subscripts are general and not an artefact of the crystallographic range,
with `13` and `17` prime and `16`, `24`, `30` highly composite.

A check that names no set runs on `default`.

## Exhaustive, never sampled

Every check is exhaustive over its domain: all pairs for the group axioms
and the homomorphism check, all triples for the cocycle condition, every
element of every group in the set. The "sampled beyond `|G| ≤ 48`" clauses
of the abstract-group and double-group chapters are superseded by this
chapter: a sampled check is not a proof, and the cost is bounded (the
largest group of `default` is `Ih†` with 240 elements; `D30h†` in `large`
has 240 as well). Checks that are slow on the largest groups are not
weakened but scheduled: they carry the `slow` marker for groups with
`|G†| > 96`.

## Check ledger

The ledger is the complete list of checks. An implementation session turns
it into the parameter table of the test suite one row at a time; each row
names the chapter that defines the check, the oracle it uses, and the
group set and marker it runs under.

Oracles: **self** (a property of the tables alone), **seed** (the
checked-in tables of known values), **reference** (transcribed reference
data under `tests/reference/`), **schema** (`schema.json`), **golden**
(checked-in documents), **formula** (a closed form of the literature
implemented as a test oracle under `tests/`, never shipped).

Markers: none (the default battery), `slow`, `reference` (skipped with a
reason when the reference file it needs is absent). The `--verify` column
says whether the flag runs the check.

### Groups and classes

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `AG-1` | [Abstract groups](./abstract-groups.md#verification) | Cayley table of `G` satisfies closure, identity, inverses, associativity (all triples) | self | `default`, `large` | `slow` for `large` | yes |
| `AG-2` | Abstract groups | `(order, class count, histogram of element orders)` equals the seed table | seed | `seed` | | yes (`seed` groups only) |
| `AG-3` | Abstract groups | canonical element order is proper-then-improper with contiguous classes | self | `default`, `large` | | yes |
| `DG-1` | [Double groups](./double-groups.md#factor-system) | `ω(E, g) = ω(g, E) = 1` and every value is `1` or `−1` | self | `default`, `large` | | yes |
| `DG-2` | Double groups | cocycle condition over all triples | self | `default`, `large` | `slow` for `\|G†\| > 96` | yes |
| `DG-3` | Double groups | `ω(g, g⁻¹) = ω(g⁻¹, g)` | self | `default`, `large` | | yes |
| `DG-4` | Double groups | Cayley table of `G†` satisfies the group axioms; `\|G†\| = 2\|G\|`; `~E` is central of order 2 | self | `default`, `large` | `slow` for `large` | yes |
| `DG-5` | Double groups | `ω` depends on rotation parts only; for a type II group it is the pull-back of `ω` of `G⁺`; `ω(C2, C2) = −1` for every binary rotation | self | `default`, `large` | | yes |
| `DG-6` | Double groups | class counts of `G†` and the number of projective irreps equal the seed table | seed | `seed` | | yes (`seed` groups only) |
| `CL-1` | [Conjugacy classes](./conjugacy-classes.md#verification) | classes partition the elements; sizes divide `\|G\|`; `size × \|centraliser\| = \|G\|` | self | `default`, `large` | | yes |
| `CL-2` | Conjugacy classes | class labels match the grammar and are pairwise distinct | self | `default`, `large` | | yes |
| `CL-3` | Conjugacy classes | `regular` classes are exactly the split ones; their number equals the number of projective irreps | self | `default`, `large` | | yes |
| `CL-4` | Conjugacy classes | `inverse` is an involution; `base_class`, `split`, `barred` agree with the element records | self | `default`, `large` | | yes |
| `CL-5` | Conjugacy classes | polyhedral class labels reproduce the structures listed in [Label conventions](./label-conventions.md) | self | polyhedral groups | | yes |

### Irrep labels and characters

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `IL-1` | [Irrep labels](./irrep-labels.md#verification) | exactly one canonical symbol per irrep, grammar-conformant, pairwise distinct | self | `default`, `large` | | yes |
| `IL-2` | Irrep labels | the symbol is a function of the character (recomputed from the character table it is reproduced) | self | `default`, `large` | | yes |
| `IL-3` | Irrep labels | `dimension = χ(E)`; `frobenius_schur = (1/\|G\|) Σ ω(g, g) χ(g²)`; `partner` present iff the indicator is `0` and the conjugate is a different irrep | self | `default`, `large` | | yes |
| `IL-4` | Irrep labels | canonical irrep order is the one the chapter fixes | self | `default`, `large` | | yes |
| `CH-1` | this chapter | row orthogonality: `Σ_g χ_i(g) χ_j(g)* = \|G\| δ_ij` over `G` (ordinary) and over `G†` (all irreps of `G†`) | self | `default`, `large` | | yes |
| `CH-2` | this chapter | column orthogonality: `Σ_i χ_i(C) χ_i(C')* = \|centraliser\| δ_CC'` | self | `default`, `large` | | yes |
| `CH-3` | this chapter | `Σ d² = \|G\|` over the ordinary irreps and again over the projective irreps; the number of ordinary irreps equals the class count of `G` | self | `default`, `large` | | yes |
| `CH-4` | this chapter | characters are class functions of `G†`; `χ(~g) = +χ(g)` for ordinary and `−χ(g)` for projective irreps; projective characters vanish off regular classes | self | `default`, `large` | | yes |
| `CH-5` | this chapter | Frobenius–Schur consistency: the indicator is `±1` iff the character is real, `0` iff it is not; `Σ_i ν_i d_i = #{x ∈ G† : x² = E}` over all irreps of `G†` | self | `default`, `large` | | yes |
| `CH-6` | this chapter | characters equal those computed independently by the Burnside class-matrix method | self | `default`, `large` | `slow` for `\|G†\| > 96` | no |
| `VL-1` | [Values](./values.md) | every exact value lies in the declared field and its coefficient vector has length `φ(M)` | self | `default`, `large` | | yes |
| `VL-2` | Values | every numeric value is the correctly rounded double of its exact value | self | `default` | | yes |
| `VL-3` | Values | the declared conductor is minimal: no proper divisor `M'` of `M` contains every value of the unfiltered document (`--irreps all`) | self | `default`, `large` | `slow` for `large` | no |
| `VL-4` | [Group chains](./group-chains.md#exact-computation) | every squared norm of every chain label of every irrep of every chain class of `list chains`, and of every Clebsch–Gordan pair, of the seven polyhedral groups is a square of the declared field, equals `1/2` with factor `(1 + i)/2`, or has a shipped norm-factor entry; every entry satisfies `μ · conj(μ) = N` exactly; every printed row of the norm-factor table of the chapter is reproduced by the shipped data | self | `large` | `slow` | no |
| `VL-5` | Group chains | the declared conductor of a `--chain` document equals the chain conductor of its chain class (`M_G` for the axial families at `n ≤ 12`; the chapter's table for the polyhedral groups) and `list chains` prints it as `conductor`; `VL-3` holds for the unfiltered `matrices` document of every chain class | self | `default` (axial, `T`, `O`), `large` (every class) | | no |
| `VL-6` | Group chains | the exact real square root in `Q(ζ_M)⁺` returns the root for `2` in `Q(ζ_8)`, `3` in `Q(ζ_12)`, `(2 + √3)/6` in `Q(ζ_12)`, `5` in `Q(ζ_5)` and `7/20 · 20/7 = 1` in every field of the conductor table of the chapter, proves non-squareness for `2` in `Q(ζ_12)`, `7` in `Q(ζ_120)`, `(3 ± √3)/6` in `Q(ζ_24)`, and agrees with an independent high-precision evaluation (`mpmath`, a test dependency only) on random squares of `Q(ζ_24)`, `Q(ζ_60)`, `Q(ζ_120)` with coefficient denominators up to `10⁶` | self, `mpmath` | `default` | | no |

### Irrep matrices

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `MX-1` | [Basis conventions](./basis-conventions.md#verification) | homomorphism by closure: every element reached by two words gets the same matrix (up to `ω` for projective irreps) | self | `default`, `large` | | yes |
| `MX-2` | Basis conventions | `D(g) D(h) = ω(g, h) D(gh)` for all pairs and all irreps | self | `default`, `large` | `slow` for `\|G†\| > 96` | yes |
| `MX-3` | Basis conventions | every matrix is unitary under `ζ ↦ ζ⁻¹` | self | `default`, `large` | | yes |
| `MX-4` | Basis conventions | trace of every matrix equals the character | self | `default`, `large` | | yes |
| `MX-5` | Basis conventions | canonical basis: `D(c)` is diagonal with the prescribed eigenvalue order | self | `default`, `large` | | yes |
| `MX-6` | Basis conventions | reference basis convention: the intertwiner to the canonical realisation exists and is unique up to a phase | self | groups the convention covers | | yes |
| `MX-7` | Basis conventions | every entry lies in the declared field | self | `default`, `large` | | yes |
| `MX-8` | Basis conventions | irreps are pairwise inequivalent and irreducible: `(1/\|G\|) Σ_g \|χ(g)\|² = 1` | self | `default`, `large` | | yes |

Under `--chain`, `MX-1` to `MX-4`, `MX-7` and `MX-8` run unchanged on the
chain-adapted matrices, in the test suite for every chain of `list chains`
of every group of the `default` set (the group set of `CB-1`) and under
`--verify` on the document's chain; `MX-5` does not run and `OS-8` takes
its place ([Where checks run](#where-checks-run)).

### Group chains and chain-adapted bases

Decided in ticket
[Verification of chain-adapted bases and CG tables](https://github.com/hydeik/pointgroup-tables/issues/48).
"Every chain" means every chain class of `list chains` for the group, in
its canonical spelling.

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `GC-1` | [Group chains](./group-chains.md#verification) | the chain grammar table of `(string, group, exit code)`, including `--chain` with a `--basis-convention` other than `canonical` → `2` and `--chain` on a utility subcommand → `2` | self | as listed (`O`, `T`, `I`, `D6h`, `D3h`) | | no |
| `GC-2` | Group chains | `--chain default` resolves to the default-chain table, group by group for the polyhedral groups and for the axial families at `n = 1 … 12` (the empty chain for the cyclic shape) | self | `default` | | no |
| `GC-3` | Group chains | the embedding examples: the resolved subgroup of `"C2"` on `D6` is generated by `C2'(1)`, of `"C1h"` on `C6v` by `sigma_v(1)`, of `"D4>C2"` on `O` by `C2(2)`, of `"D3[C2'']>C3"` on `D6` by `C3^1` and `C2''(1)`; the resolved `D4` of `O` contains `C4^1(1)` and the resolved `C3v` of `Td` contains `C3^1(1)` | self | as listed | | no |
| `GC-4` | Group chains | `list chains` agrees with the subgroup catalogue of ticket #41 (canonicity for both kinds recomputed from the characters, terminal rule, `default` flag, `conductor`) **modulo the twin chain classes**: one entry per canonical spelling, the catalogue's chain classes collapsing exactly by the chapter's table of [twin chain classes](./group-chains.md#embedding-which-subgroup-a-symbol-picks) (none in the axial families at `n ≤ 12`); every canonical spelling resolves to the class it was derived from or to its reached twin, and every spelling of the table resolves to the reached twin named there; `list chains O` holds `"C4"` first with `default = true` | self | polyhedral groups, `sweep` | | no |
| `GC-5` | Group chains | the metadata of every `--chain` document holds `chain = { input, canonical, symbols }`, and `--chain` with the canonical spelling of a document's chain reproduces its `symbols`, `subgroup_chain` table and matrices byte for byte | self | `seed`, every chain | | no |
| `CB-1` | [Group chains](./group-chains.md#properties) | intertwiner: for every irrep the basis matrix `U` (columns the chain-adapted vectors in canonical components) is unitary and `D' = U† D U` with the canonical `D` (same characters); the label spaces have rank exactly `1` and their ranks sum to the dimension | self | `default`, every chain; `large` | `slow` for `large` | no |
| `CB-2` | Group chains | sub-chain block equality: for every step-1 label `σ`, the block of `D'_Γ(φ_1(h))`, `h ∈ H_1,lib†`, equals `D'_σ(h)` of `H_1,lib` under the pulled-back chain entry by entry (`matrices O --chain "D4>C4"` contains the `E` blocks of `matrices D4 --chain C4`) | self | `default`, every chain of two or more steps; `large` | `slow` for `large` | yes |
| `CB-3` | Group chains | the norm-factor phase rule: every designated component equals the conjugate of the canonical norm factor of its squared modulus, real positive wherever the real square root lies in the document field; within a block the vectors are those making the restriction equal to the sub-chain realisation | self | `default`, every chain; `large` | `slow` for `large` | no |
| `CB-4` | Group chains | chains through the principal element: `matrices` under the default chain equals the canonical table entry by entry, with `components` added, for every axial family (`D6 --chain C6`, `D3h --chain C3h`, `D6h --chain C6h`, …) and for `T>C3`, `Th>S6`, `Td>S4`, `O>C4`, `Oh>C4h`; every chain vector is a canonical vector with phase `1` and every designated component has `N = 1`; `O "D4>C4"`, `Td "D2d>S4"`, `Oh "D4h>C4h"` differ from the canonical vectors by `diag(1, ∓i)` on the `E1/2`/`E3/2` blocks of the spinor irreps, `I "D5>C5"` and `Ih "D5d>S10"` by the phases and the `I5/2` split stated in [Properties](./group-chains.md#properties) (research of ticket #49, order rule of ticket #57) | self | `sweep`, the seven polyhedral groups | | no |
| `CB-5` | Group chains | centro and twisted inheritance: `matrices Oh --chain "D4h>C4h"` equals `matrices O --chain "D4>C4"` with the parity signs (`D'(i·r) = ±D'(r)`), and the chain-adapted matrices of every type III group equal those of its rotation group's chain on the image elements (`D3h --chain C3h` from `D6 --chain C6`, `Td --chain S4` from `O --chain C4`) | self | the centro and twisted groups of `seed`, every chain | | no |

### Clebsch–Gordan coefficients

The checks are defined in
[Clebsch–Gordan coefficients](./clebsch-gordan.md#properties); `CG-8` of
that chapter is the reference row `REF-11` below.

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `CG-1` | [Clebsch–Gordan coefficients](./clebsch-gordan.md#properties) | unitarity: `U† U = 1` exactly for every pair | self | `default`, `large` | `slow` for `\|G†\| > 96` | yes |
| `CG-2` | Clebsch–Gordan coefficients | reduction: `(D_α(g) ⊗ D_β(g)) U = U (⊕_{γ,r} D_γ(g))` for every `g` of the **canonical generating set** of `G`, with the matrices of the document's basis and the canonical section. Both sides are projective representations with the factor system `ω_α ω_β` (`MX-2` proves it for the irrep matrices), so equality on a generating set is equality on `G`: a proof, not a sample | self | `default`, `large` | `slow` for `\|G†\| > 96` | yes |
| `CG-3` | Clebsch–Gordan coefficients | multiplicities: the columns of `γ` number `n_γ d_γ` with `n_γ` from the characters, `Σ_γ n_γ d_γ = d_α d_β`; for `α = β` the split `n⁺_γ + n⁻_γ = n_γ` holds and the copies are counted per part; the triples, splits and squared norms of the [verification data](./clebsch-gordan.md#verification-data) are reproduced | self | `default`, `large` | | yes |
| `CG-4` | Clebsch–Gordan coefficients | exchange eigenvalue: for `α = β` every column is an eigenvector of the swap `S` with the eigenvalue its `exchange` record states, `+1` for the symmetric copies, which come first | self | `default`, `large` | | yes |
| `CG-5` | Clebsch–Gordan coefficients | designated coefficient: in every column `(γ, r, k₀ = 0)` the first non-zero entry in product order is the conjugate of the canonical norm factor of its squared modulus, real positive whenever `√N` lies in the field; for the symmetric copies of `α = β` the first non-zero entry of the symmetrised column (canonical basis; under `--chain` or a reference basis the coupling matrix equals the unitary transform of the canonical one entry by entry) | self | `default`, `large` | | yes |
| `CG-6` | Clebsch–Gordan coefficients | field: every entry lies in the declared field and the declared conductor is the `M_cg` of the chapter's table (`VL-1`, `VL-3` run on the `cg` document as on any table with values) | self | `default`, `large` | | yes |
| `CG-7` | Clebsch–Gordan coefficients | conjugation: `n_{γ̄}(ᾱ, β̄) = n_γ(α, β)`, and the isotypic projector of `γ̄` in `ᾱ ⊗ β̄` equals the conjugate of that of `γ` in `α ⊗ β` transported by the `MX-6` intertwiners | self | `default`, `large` | `slow` for `\|G†\| > 96` | yes |
| `CG-9` | Clebsch–Gordan coefficients | reversed pair: the coupling matrix of `(β, α)` the implementation produces is the row permutation of that of `(α, β)` from the `(α, β)` to the `(β, α)` row-major order, scaled by `ε_r` per column for `α = β`; a `cg` document prints `α ≤ β` only | self | `default`, `large` | | no |

### Documents and command line

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `OS-1` | [Output schema](./output-schema.md#verification) | every document of every subcommand validates against `schema.json`, for every combination of `--double-group`, `--irreps`, `--values`, with and without `--chain default` | schema | `seed` | | no |
| `OS-2` | Output schema | `tables` and `values` list exactly what is present; all forms of a table have the same shape | self | `seed` | | no |
| `OS-3` | Output schema | every index is in range for its index space; `items[i].index = i`, `classes[c].index = c`; canonical irrep order in `characters` and `matrices` | self | `seed` | | no |
| `OS-4` | Output schema | irrep records of `characters` and `matrices` agree in `all`; `characters.classes` equals the classes of the matching `elements` document | self | `seed` | | no |
| `OS-5` | Output schema, [Command line](./cli.md#verification) | `all` is the union of the five single-table documents; every single-table document embeds the same `elements` | self | `seed` | | no |
| `OS-6` | Output schema, Command line | determinism: two runs with the same arguments are byte-identical, every subcommand (`cg` included) and format, with and without `--chain` | self | `seed` | | no |
| `OS-7` | [Group chains](./group-chains.md#verification-of-the-table) | `subgroup_chain` under `--chain`: present with `chain` and second in `tables`, byte-identical across subcommands and `--double-group`, `--irreps`, `--values`; `irreps` is the unfiltered list of `G`; per step `embedding` is an injective double-group homomorphism into the previous group (Cayley table of `G†`), `cosets` partitions the previous group into left cosets of the recorded subgroup with the stated ordering, `double_cosets` partitions `cosets` with `size` and closure under the subgroup on both sides, `branching` equals the subduction recomputed from the characters through `embedding`, without repeats | self | `seed`, every chain of `list chains` | | yes |
| `OS-8` | Group chains | `components` under `--chain`: `dimension` records per `matrices` entry, `index = position`, `labels` in range and in the chapter's component order (lexicographic in the steps, exponent order at a cyclic step), consistent with `branching` step by step; restricted to each step the matrices are block-diagonal by the labels up to that step, and diagonal on the terminal with the characters of the terminal labels; absent without `--chain` | self | `seed`, every chain of `list chains` | | yes |
| `OS-9` | [Clebsch–Gordan coefficients](./clebsch-gordan.md#document-checks) | `cg` document shape: validates against `schema.json`; `tables` is `["elements", "cg"]` or `["elements", "subgroup_chain", "cg"]`; `double = false`; `irreps` present; `field.conductor` is the `M_cg` of the chapter's table; `cg.irreps` is the full irrep list of `G`, equal to the unfiltered `characters.irreps` of the same options, with `components` under `--chain` and absent otherwise | schema, self | `seed`, every `--irreps` and `--values`, with and without `--chain default` | | no |
| `OS-10` | Clebsch–Gordan coefficients | `cg` pairs and columns: `pairs` holds exactly the pairs `α ≤ β` of the selection in lexicographic order, the pair sets of `ordinary` and `projective` are disjoint with union `all`; per pair `dimension = d_α d_β`, every form a `dimension × dimension` matrix, `columns` has `dimension` records with `index` equal to position ordered by `irrep` (canonical order), `copy`, `component`; `constituents` lists exactly the `γ` with `n_γ > 0` with `columns = [first, n_γ d_γ]` contiguous and in order; `exchange` present in every column record exactly when `α = β` | self | `seed`, every `--irreps` and `--values`, with and without `--chain default` | | no |
| `OS-11` | Clebsch–Gordan coefficients, Command line | `cg --double-group` gives `2`; the `all` document of the same options is unchanged (`tables`, `field`, bytes) by the existence of `cg` | self | `seed` | | no |
| `CLI-1` | Command line | the group-argument table of `(exit code, symbol)`, for the version-1 table subcommands and for `cg` | self | as listed | | no |
| `CLI-2` | Command line | unknown convention gives `4`; option on the wrong subcommand gives `2`; on failure the document stream is empty and no file is created | self | | | no |
| `CLI-3` | Command line | `--verify` on a passing group emits the same bytes as without it; on `list`/`schema` gives `2`; on a failing check gives `1`, no document, one `verify:` line per failing check | self, mutation | `C3v`, `Oh` | | no |
| `GOLD-1` | this chapter | golden documents: every table subcommand (`cg` included) × `json`, `yaml`, `toml` for `C3v`, `D6h`, `Oh`, `Ih` with default options equals the checked-in file byte for byte | golden | those four | | no |
| `GOLD-2` | this chapter, [Group chains](./group-chains.md#verification-of-the-table) | chain golden documents: every table subcommand (`cg` included) × `json`, `yaml`, `toml` for `C3v`, `D6h`, `Oh`, `Ih` under `--chain default`, plus `matrices O --chain "D4>C4"`, `matrices D6 --chain "D3[C2'']>C3"`, `cg T --chain C3` and `cg O --chain C4` (the chains of the papers' CG tables) in the three formats, equal the checked-in files byte for byte | golden | those seven | | no |

### Reference comparison

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `REF-1` | this chapter | structural comparison of the character tables of `G` and `G†` with Altmann & Herzig | reference | groups transcribed | `reference` | no |
| `REF-2` | this chapter | structural comparison of the character tables of `G` and `G†` with Bilbao | reference | 32 crystallographic groups | `reference` | no |
| `REF-3` | this chapter | labelled comparison of the character tables with Altmann & Herzig under the `altmann-herzig` label convention | reference | groups transcribed | `reference` | no |
| `REF-4` | [Double groups](./double-groups.md#verification-data) | `ω(g, h) = −1` exactly where the Altmann & Herzig multiplication table (the `[multiplication]` table of the labels file) marks the product with a tilde, after the convention's `barred` swaps | reference | `C2v`, `D3`, `D4`, `D3h`, `T`, `O`, `Ih` | `reference` | no |
| `REF-5` | Irrep labels | canonical symbols after the `altmann-herzig` convention equal the book's symbols row by row | reference | groups transcribed | `reference` | no |
| `REF-6` | Label conventions | the `altmann-herzig` label convention reproduces the transcribed element labels of every group, element by element | reference | groups transcribed | `reference` | no |
| `REF-7` | Basis conventions | under `--basis-convention altmann-herzig`, the matrices of the transcribed generators equal the book's entry by entry (polyhedral groups from the polyhedral transcription, the ten axial tables from the matrix transcription), and the Wigner matrices of every ket the book lists in its symmetrized-bases tables (`tests/reference/altmann-herzig/bases/`) equal the spin-harmonic rule's | reference | groups transcribed | `reference` | no |
| `REF-8` | Irrep labels | the shipped Koster table: each `Γ_i` entry's identifying characters match the computed irrep it is attached to | reference | 32 crystallographic groups | `reference` | no |
| `REF-9` | this chapter | Chen–Fan–Ping label maps: the transcribed multiplication tables with tildes of `T†`, `O†` and `D6†` (Bradley–Cracknell labels, `tests/reference/chen-fan-ping/labels/`) determine an isomorphism onto `G†` and the section signs `ε`, and `ω(g, h) = −1` exactly where the table marks a tilde after `ε`; the `I†` map built from the Euler-angle tables through the Altmann & Herzig identification of the axes is a homomorphism of `G†` | reference | `T`, `O`, `D6`, `I` | `reference` | no |
| `REF-10` | [Group chains](./group-chains.md#properties), this chapter | chain-adapted matrices against the papers ([comparison up to a monomial intertwiner](#comparison-up-to-a-monomial-intertwiner)): the transcribed matrices, read through the label map and `ε`, form a representation with the tool's factor system; the intertwiner to `matrices <group> --chain <chain>` exists and is monomial; its permutation equals the `m = −e` label map with the survey's exceptions; its phases are computed and reported, never stored; `known_discrepancy` records where a table differs | reference | the chains transcribed: `T>C3`, `T>D2>C2`, `O>C4`, `O>T>C3`, `O>T>D2>C2`, `O>D4>C4`, `Td>D2d>C2`, `I>D5>C5` (projective, and the ordinary irreps against the `I ⊃ C5` tables), the dihedral `⊃ C2` tables of `D4`, `C4v`, `D2d`, `D2`, `C2v`, `D6`, `C6v`, `D3h`, `D3`, `C3v` | `reference` | no |
| `REF-11` | [Clebsch–Gordan coefficients](./clebsch-gordan.md#relation-to-the-papers-conventions) | coupling matrices against the papers ([comparison of coupling matrices](#comparison-of-coupling-matrices)): after the basis layer, the isotypic projector of every constituent equals the paper's exactly; every multiplicity-free coupled vector equals the paper's up to one phase (`±1` where both are real), the phase computed and reported; for `n_γ > 1` the `n_γ × n_γ` matrix relating the paper's copies to the tool's is unitary; `known_discrepancy` records where a table differs (this is `CG-8` of the chapter) | reference | `T>C3`, `T>D2>C2`, `O>C4`, `O>T>C3` | `reference` | no |
| `REF-12` | Group chains, Clebsch–Gordan coefficients | the dihedral closed forms of Chen, Fan, McAven & Butler (2000) as a test oracle (`tests/oracle_dihedral.py`: eqs 2.12a–b for the matrices of `Dn ⊃ Cn` and `Dn ⊃ C2x`, eq 2.26 for the `Dn ⊃ Cn` coupling coefficients, the isomorphisms of eqs 3.2, 4.7, 5.4 for `Cnv`, `Dnd`, `Dnh`) agree with `matrices --chain <Cn>` (the components swapped, no phase), with `matrices --chain C2` (monomial intertwiner) and with `cg` (isotypic projectors; multiplicity-free vectors up to one phase) for every axial family and every `n` of the set | formula | `sweep`, `large` | `slow` for `large` | no |

### Mutation tests

| ID | Check | Oracle | Marker |
| -- | ----- | ------ | ------ |
| `MUT-1` | one Cayley entry replaced by another element is caught by `AG-1` | mutation | |
| `MUT-2` | one `ω(g, h)` sign flipped is caught by `DG-2` and by `MX-2` | mutation | |
| `MUT-3` | one matrix multiplied by a phase `ζ` is caught by `MX-1` or `MX-2` and `MX-4` | mutation | |
| `MUT-4` | one character replaced by another value of the field is caught by `CH-1` and `CH-2` | mutation | |
| `MUT-5` | two class labels swapped is caught by `CL-4` or `CL-2`; two irrep symbols swapped by `IL-2` | mutation | |
| `MUT-6` | one element moved to another class is caught by `CL-1` and `CH-4` | mutation | |
| `MUT-7` | one row of a reference character table altered is caught by `REF-1` | mutation | `reference` |
| `MUT-8` | the `labels` of two `components` records of one `matrices` entry swapped is caught by `OS-8` | mutation | |
| `MUT-9` | one `embedding` image of a subgroup record replaced by another element is caught by `OS-7` | mutation | |
| `MUT-10` | one chain-adapted basis vector multiplied by a phase is caught by `CB-3`, and by `CB-2` when the vector lies in a block of dimension `≥ 2` (the `E` block of `O --chain "D4>C4"`) | mutation | |
| `MUT-11` | one Clebsch–Gordan coefficient's sign flipped is caught by `CG-2` (and `CG-1`) | mutation | |
| `MUT-12` | one column of a coupling matrix multiplied by a phase is caught by `CG-5` | mutation | |
| `MUT-13` | two copies of a multiplicity-two triple swapped (`G3/2g ⊗ G3/2g ⊃ T1g` of `Oh`) is caught by `CG-4` or `CG-5` | mutation | |
| `MUT-14` | one entry of a transcribed paper matrix altered is caught by `REF-10`; one transcribed coupling coefficient altered by `REF-11` | mutation | `reference` |

A mutation test builds the document model of `C3v` and `Oh` (one axial, one
polyhedral group, both with projective irreps), applies the mutation to the
model, and runs the checkers on the corrupted model. The version-2
mutations build, in addition, the `--chain default` documents of `C3v` and
`Oh`, `matrices O --chain "D4>C4"` (a two-step chain with a two-dimensional
block) and the `cg` documents of `C3v` and `Oh` (one axial, one with
multiplicity-two triples). Every check with
`--verify = yes` in the ledger must be the target of at least one mutation;
the suite fails otherwise. This is the proof that the proof works.

## Reference comparison

### Structural comparison

A **structural comparison** matches a computed character table with a
reference table without using any label. Every reference class carries the
invariants the reference lets one read off: its `size`, its `parity`
(proper or improper), and, when the reference lists its elements, the
`order` of its elements in the group of the table. The comparison succeeds
when there is a bijection of classes that preserves the invariants and a
bijection of irreps such that the two tables agree entry by entry, exactly.
The search is a backtracking over the small candidate sets the invariants
leave; a match must exist, and the test reports the matched bijection so a
human can inspect it.

The structural comparison is available for every reference and needs no
label convention, so it runs as soon as the reference file exists. It
proves the mathematical content of the transcription (and thus also catches
transcription errors, which is why a second, independent reference is
kept), but not that the tool's names agree with the reference's.

### Labelled comparison

A **labelled comparison** applies the selected label convention to the
computed tables, then matches rows to the reference by irrep label and
columns by class label, and requires the entries to agree exactly. It runs
only for a reference that has a label convention (`altmann-herzig` once the
convention data exists) and is what proves the convention: element labels,
irrep symbols, class membership and the barred assignment. Comparing a
double-group table requires the convention's `barred` swaps to be applied
first ([Conjugacy classes](./conjugacy-classes.md)).

### Comparison up to a monomial intertwiner

Decided in ticket
[Verification of chain-adapted bases and CG tables](https://github.com/hydeik/pointgroup-tables/issues/48).
The chain-adapted matrices of the papers of Chen, Fan, Ping and co-workers
(the survey of ticket #39, `docs/research/eigenfunction-method.md`) are
written in the chain-adapted basis of the same chain with the papers' own
phases and component order. Two chain-adapted realisations of one irrep for
one canonical chain differ by a **monomial intertwiner**: a permutation of
the components (the label map) times a diagonal phase matrix, because every
label space is one-dimensional. The comparison (`REF-10`) therefore runs
in four steps, per chain and per irrep:

1. **Label map and section.** The paper's element labels (Bradley &
   Cracknell) are mapped to canonical elements by the isomorphism the
   transcribed multiplication table determines (`REF-9`, the analogue of
   `REF-4` and `REF-6` for Altmann & Herzig), with the section signs `ε`
   of [Double groups](./double-groups.md#section-overrides-in-label-conventions) where the
   paper's untilded element is the tool's barred one (the dihedral paper's
   `C_n^k`, `n/2 < k < n`, and the reversed binary axes `C2b`, `C2c`, `C2f`
   of the papers' `O†` table).
2. **Representation.** The transcribed matrices, with `ε`, multiply with
   the tool's factor system on every transcribed product (generators and
   their listed words), and their traces are the characters of the irrep
   the file's `[canonical]` map names.
3. **Intertwiner.** The unique-up-to-phase intertwiner `W` between the
   paper's realisation and `matrices <group> --chain <chain>` of the same
   irrep is computed exactly (the `MX-6` primitive) and must be
   **monomial**; its permutation must equal the label map of the survey
   (`m = −e`, the components reversed, with the exceptions the survey
   lists per chain) composed with the tool's component order (exponent
   order at a cyclic step); the moduli of the entries therefore agree entry by entry after
   the permutation (the label-free content).
4. **Phases.** The diagonal of `W` is computed and reported by the test,
   never stored as data: the papers fix their phases ad hoc, and
   [ADR-0011](../adr/0011-chain-basis-phase-rule.md) declines to follow
   them.

The comparison is exact at every step; a printed table that fails step 2
or 3 is wrong in the print or in the transcription and is resolved by a
`known_discrepancy` record (below), never by weakening the check.

### Comparison of coupling matrices

The papers' Clebsch–Gordan tables relate to the tool's in the three layers
of [Clebsch–Gordan coefficients](./clebsch-gordan.md#relation-to-the-papers-conventions):
basis, one phase per copy, a mixing of copies. The comparison (`REF-11`)
applies them in that order, per pair and per constituent `γ`:

1. **Basis layer.** The paper's product basis is mapped to the tool's with
   the monomial intertwiners `W_α`, `W_β` of step 3 above on the rows
   (`W_α ⊗ W_β`, after the row permutation from the paper's component
   order to row-major product order) and `W_γ` on the columns. This is why
   a CG comparison needs the matrices of the same chain transcribed first.
2. **Isotypic projector.** After layer 1, the projector
   `Π_γ = Σ_{r,k} |γ k; r⟩⟨γ k; r|` built from the paper's columns equals the
   tool's exactly. This is the label-free content: it is independent of
   the phases, of the copy labels and of any orthogonalisation of the
   copies, and holds for `n_γ = 1` and `n_γ > 1` alike.
3. **Multiplicity-free triples.** For `n_γ = 1` the paper's coupled vector
   equals the tool's up to one phase, `±1` wherever both are real; the phase
   is computed and reported (a labelled comparison of the one column).
4. **Multiplicity.** For `n_γ > 1` the `n_γ × n_γ` matrix expressing the
   paper's (orthogonalised) copies in the tool's copies is computed and
   must be unitary; it is reported, never stored. The papers' copies are
   labelled by trial states and the tool's by the symmetric/antisymmetric
   split and trial-state order, so entry-by-entry agreement is neither
   expected nor required.

Printed tables that the tool cannot reach are not compared: the `T ⊃ C2`
tables need a chain that is canonical for the projective irreps only
(`--chain C2` on `T` exits `6`), and the `I ⊃ C5` tables a chain that is
not canonical for `I5/2`; the ordinary irreps of the `I ⊃ C5` tables are
compared against `matrices I --chain "D5>C5"`, whose label spaces for the
ordinary irreps are the `C5` eigenspaces. No paper prints Clebsch–Gordan
coefficients of `I`.

### References and what they check

| Reference | Tables | Groups | Mode |
| --------- | ------ | ------ | ---- |
| Altmann & Herzig, *Point-Group Theory Tables*, 2nd ed. | character tables of `G` and `G†`; element labels; irrep symbols; tilde-marked products; irrep matrices of the polyhedral groups plus axial samples | the 32 crystallographic groups, `I`, `Ih`, `C5v`, `D5h`, `D5d`, `S8`, `D4d` | structural, then labelled |
| Bilbao Crystallographic Server, point-group and double-point-group tables | character tables of `G` and `G†` | the 32 crystallographic groups | structural only (no Bilbao label convention in version 1) |
| Koster, Dimmock, Wheeler & Statz | Bethe symbols | the 32 crystallographic groups | by identifying characters (`REF-8`) |
| Chen, Fan, Ping and co-workers, the papers of the eigenfunction method (keys `[CF98a]`, `[CF98b]`, `[CF98c]`, `[FC99]`, `[F03]`, `[F01]`, `[PC99]`, `[CFP00]`, `[CP99]`, `[FCD99a]`, `[FCD99b]`, `[FCMB99]`, `[CFMB00]` of the survey) | multiplication tables with tildes of `T†`, `O†`, `D6†`; Euler angles of `I†`; chain-adapted irrep matrices of the chains listed under `REF-10`; Clebsch–Gordan tables of `T ⊃ C3`, `T ⊃ D2 ⊃ C2`, `O ⊃ C4`, `O ⊃ T ⊃ C3`; the dihedral closed forms (as a test oracle, not a transcription) | `T`, `O`, `Td`, `I`, `D2`, `C2v`, `D2d`, `D3`, `C3v`, `D3h`, `D4`, `C4v`, `D6`, `C6v`; every axial family through the oracle | up to a monomial intertwiner (`REF-10`); isotypic projectors, then one phase (`REF-11`); formula (`REF-12`) |

Bradley & Cracknell is not transcribed; it appears only in the notes of the
Koster table where its `Γ_i` numbering disagrees, and as the element labels
of the Chen–Fan–Ping papers, whose own multiplication tables fix the map
to the canonical labels (`REF-9`).

The secondary references for Clebsch–Gordan coefficients (the coupling
tables of Koster, Dimmock, Wheeler & Statz, those of Altmann & Herzig where
the book prints them, Butler's tables) are not compared in version 2: each
needs a basis transformation onto the tool's basis that the specification
does not state ([Command line](./cli.md#not-in-version-2)).

## Reference data

Reference data is **test data** and lives under `tests/reference/`, never
in the package. Data the tool ships (label and basis conventions, the
Koster table) lives under `src/pointgroup_tables/data/conventions/` and is
derived from the transcriptions, not identical to them: a convention is
rules plus overrides, a transcription is the book as printed. The
`altmann-herzig` label convention is written by
`tools/ah_convention.py generate` and checked against every transcription
by `tools/ah_convention.py check` (the prototype of `REF-6`).

```
tests/reference/
  README.md                          # sources, citation, licensing note (ADR-0006)
  altmann-herzig/
    labels/<group>.toml              # element labels, the book's identification of each operation,
                                     # and the book's multiplication table with its tilde marks
    characters/<group>.toml          # character tables of G and G†
    matrices/<group>.toml            # generator matrices (polyhedral groups, axial samples)
  bilbao/
    fetch.py                         # the script that produced the files; run once, not by the tests
    characters/<group>.toml          # the Point Group Tables page: classes of G, ordinary irreps
    characters/<group>-double.toml   # the REPRESENTATIONS DPG page: classes of G† listed element by element, every irrep
  koster/
    bethe/<group>.toml               # Γ_i per canonical symbol with identifying characters
  chen-fan-ping/
    README.md                        # paper keys, the chains covered, the label maps
    labels/<group>.toml              # Bradley–Cracknell element labels and the papers' multiplication
                                     # table with tildes (T†, O†, D6†); Euler angles for I†
    matrices/<group>.chain-<spelling>.toml   # chain-adapted matrices of the transcribed generators
    cg/<group>.chain-<spelling>.toml         # Clebsch–Gordan tables, one file per chain
```

Every file is TOML, one per group, named by the normalised Schoenflies
symbol, and starts with the header the label transcription established
(`tests/reference/altmann-herzig/labels/README.md`):

```toml
group = "D3"
table = "T 20.1, p. 176"     # where in the source; a URL for Bilbao
status = "done"              # todo | partial | done
notes = ""                   # anything that could not be read unambiguously
```

The tilde-marked products stay inside the labels file (its
`[multiplication]` table), because the isomorphism between the book's
labels and the canonical ones is derived from that table; `REF-4` reads it
from there. Each kind documents its own body in a `README.md` next to the
files.

Values use the value grammar of [Basis conventions](./basis-conventions.md#value-grammar)
(`i`, `zeta`, `sqrt2`, `sqrt3`, `sqrt5`, `phi`), so characters and matrices
are parsed by the same exact parser as the convention files. A character
table file lists the reference's classes with their invariants and then
one row per irrep:

```toml
double = false

[[classes]]
label = "2C3"           # the reference's own spelling
size = 2
order = 3               # omitted when the reference does not identify its elements
parity = 1

[[irreps]]
label = "E"
characters = ["2", "-1", "0"]
bethe = "Gamma3"        # optional
```

A reference can disagree with itself (a printed table at odds with the
work's own parameter tables or master tables). Such a case is recorded in
the file as a `[[known_discrepancy]]` array of tables, one entry per
disagreement stating the correction and its evidence, never by silently
transcribing a corrected value: either the values stay as printed and the
entry gives the correction (the character tables), or the values are
corrected and the entry keeps the printed value (the matrix tables); in
both cases the entry names the affected entries. The structural comparison
(`REF-1`, `REF-2`) reports a file of the first kind as matching only after
applying the documented correction, and the resolution of the
transcription ticket lists every such case. A character-table entry is
either `swap_projective` (two class columns whose projective entries are
interchanged) or `swap_labels` (two irrep symbols printed on each other's
rows), applied in file order. First instances, both in the Altmann &
Herzig `D5h` table: its projective rows print the `2S5` and `2S5^2`
entries interchanged, and its `E3/2` and `E7/2` rows carry each other's
symbol (the book's bases T 34.6 and matrices T 34.7 decide).

### The Chen–Fan–Ping transcriptions

Decided in ticket
[Verification of chain-adapted bases and CG tables](https://github.com/hydeik/pointgroup-tables/issues/48);
the tables to transcribe are the task tickets of the map. The papers'
tables are transcribed from the user's copies under `PDF/` (gitignored)
by an agent reading the page images, as the Altmann & Herzig matrices
were, with the `transcribed_by` note and `status` of the header; a
`partial` status with `notes` records every entry the survey marked as
uncertain. The files are named by the group and the **canonical spelling
of the chain made file-safe**, the rule of the golden documents
([Test organisation](#test-organisation)): `matrices/T.chain-c3.toml`,
`matrices/T.chain-d2_c2.toml`, `matrices/O.chain-c4.toml`,
`matrices/O.chain-t_c3.toml`, `matrices/I.chain-d5_c5.toml`,
`matrices/D6.chain-c2.toml`, `cg/T.chain-c3.toml`, `cg/O.chain-c4.toml`. A
file holds one chain of one group; the header adds the chain and the
paper:

```toml
group = "O"
chain = "C4"                    # canonical spelling, as `--chain` accepts it
paper = "F03"                   # key of the survey (docs/research/eigenfunction-method.md)
table = "Tables I-III, pp. 545-549"
status = "done"
notes = ""
conductor = 120                 # every entry lies in Q(zeta_conductor); "zeta" is zeta_conductor
```

A `labels` file holds the paper's element labels, the generators the
paper's matrices are given on, and the printed multiplication table with
its tildes in the `[multiplication]` table of the Altmann & Herzig labels
files; for `I†`, whose table is not printed, it holds the Euler angles of
the sixty elements and the identification of each axis with the Altmann &
Herzig axis labels, from which the isomorphism is built (`REF-9`).

A `matrices` file keeps the paper's symbols in a `[book]` section (irrep
symbols in the paper's modified Mulliken notation, the component labels of
every irrep in the paper's order, the element labels of the transcribed
matrices, the generating set) with the paper's constants as printed and in
the value grammar (`[book.constants]`, `[book.exact]`), a `[canonical]`
map from the paper's irrep symbols to canonical symbols fixed by the
traces, and one matrix per transcribed element per irrep in the value
grammar, exactly as the Altmann & Herzig polyhedral files do
(`tests/reference/altmann-herzig/matrices/README.md`).

A `cg` file lists the paper's irrep symbols and component orders in
`[book]` and then one `[[products]]` entry per printed product and
constituent:

```toml
[[products]]
factors = ["F0", "G3/2"]        # the paper's symbols
coupled = "G3/2"
copy = 2                        # the paper's copy label, 1-based; absent when multiplicity-free
trial = ["1", "-1/2"]           # the paper's trial state when printed
rows = [["1", "-3/2"], ["1", "-1/2"], ...]   # the product states of the printed rows, paper labels
columns = { "-3/2" = ["0", "sqrt5*sqrt3/6", ...], "-1/2" = [...], ... }  # one column per coupled component
```

Values are exact in the value grammar (`sqrt(5/6)` is `sqrt5*sqrt3/6`);
the paper's row and column order is kept, and the comparison applies the
label map. A `[[known_discrepancy]]` entry follows the rule of the matrix
tables for both kinds of file: the values are **corrected** before use and
the entry keeps the printed form, names the affected entries and states
the evidence (the single correction, among signs, conjugations,
transpositions and column swaps, under which the printed table satisfies
step 2 of the [monomial comparison](#comparison-up-to-a-monomial-intertwiner)
or the isotypic-projector equality). The papers' phase conventions are
never stored: the diagonal unitary of `REF-10` and the phases and mixing
matrices of `REF-11` are computed by the test and printed in its report.

A Bilbao file is produced by `fetch.py` from the server's HTML and is then
authoritative: the tests read the TOML, never the network, and a change in
the server's pages is not a change in the data. The script stays in the
repository as the record of provenance. The server keeps its character
tables on two pages per group, so each group has a twin: `<group>.toml`
(`double = false`) from the Point Group Tables page, and
`<group>-double.toml` (`double = true`) from the REPRESENTATIONS DPG page,
whose classes are those of `G†` and carry the server's list of operations
(`elements`) and the kind of operations they hold (`barred`); its irreps
carry the server's three labels and its reality. The structural comparison
of `G†` (`tools/bilbao_characters.py check`, the prototype of `REF-2`) may
pair a lifted class of the server with a barred canonical class, since the
server's spinor convention differs from the canonical section
(`tests/reference/bilbao/characters/README.md`).

### Redistribution

The transcribed tables are redistributed with the repository under its
licence (ADR-0006). They are mathematical facts (characters, products,
labels) re-expressed in the project's own format, restricted to what the
tests need, and each file cites its source. `tests/reference/README.md`
carries the citations, the Bilbao terms of use, and a statement that the
files are not reproductions of the books' pages.

## Test organisation

- **Framework**: `pytest`, parametrised over the group sets above. No
  property-based testing library; every domain is a finite set that is
  enumerated.
- **Markers**: `slow` (large groups and the exhaustive pair checks on
  `|G†| > 96`) and `reference` (needs a file under `tests/reference/`;
  skipped with a reason naming the file when it is absent). A plain
  `pytest` runs the default battery without either; continuous integration
  runs `slow` too.
- **Golden documents**: `tests/golden/<group>/<subcommand>.<format>` for
  `C3v`, `D6h`, `Oh`, `Ih`, every table subcommand (`cg` included, its
  `Ih` files the largest at about 20 MB each), all three formats,
  default options. They are the regression and determinism oracle
  (`GOLD-1`). A script under `tools/` regenerates them; a regeneration is
  reviewed as a diff, never accepted blindly. The samples in
  `prototypes/samples/` are their seed. The version-2 chain documents
  (`GOLD-2`) sit next to them as
  `tests/golden/<group>/<subcommand>.chain-<spelling>.<format>`, where
  `<spelling>` is `default` for `--chain default` and otherwise the
  canonical spelling made file-safe (`>` to `_`, `[`, `]`, `'` and `,`
  dropped, lower case): `tests/golden/O/matrices.chain-d4_c4.json`,
  `tests/golden/D6/matrices.chain-d3c2_c3.yaml`. A version-2 tool must
  reproduce the version-1 files unchanged.
- **Test oracles**: `tests/oracle_dihedral.py` implements the closed
  forms of the dihedral paper (`REF-12`) for an arbitrary axis order `n`,
  in the exact arithmetic of the package, and is imported by the tests
  only; nothing under `tests/` is shipped, and the closed forms are never
  the implementation ([Out of scope](https://github.com/hydeik/pointgroup-tables/issues/38)).
- **Mutation tests**: as in the ledger; a fixture applies one mutation to
  the document model and asserts the expected check identifiers fail and
  no other check crashes.
- **`--verify` tests**: `CLI-3` runs the flag end to end on a passing group
  and, through an injected mutation, on a failing one.

## Acceptance criterion

The specification is **implementable** when all of the following hold:

1. **The ledger is complete.** Every item of every chapter's "Verification"
   section has a row in the check ledger with an identifier, an oracle, a
   group set and a marker. An implementation session can transcribe the
   ledger into the test suite without inventing a check.
2. **The schema is covered.** Every field of `schema.json` is defined by
   a chapter of `docs/spec/`.
3. **Nothing is deferred.** No chapter contains a `TODO`, `TBD`, or a
   decision deferred to another ticket. References to the open
   transcription tasks (element labels, character tables and axial matrices,
   polyhedral matrices, Koster table, Bilbao fetch, and the Chen–Fan–Ping
   label maps, matrices and coupling tables of version 2) are allowed: the tests
   that need their data carry the `reference` marker and skip until the
   data exists, so the transcriptions are not a precondition of
   implementation.
4. **The seed documents validate.** For the 43 groups of `seed`, the
   documents of every subcommand validate against `schema.json`, with and
   without `--chain default`.

The transcriptions themselves are not on this list on purpose: they are
human work with the books at hand and may finish after the first
implementation; a `reference` test turns from skipped to green or red the
day its file lands.
