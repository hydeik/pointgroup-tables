# Group chains

Decided in wayfinder ticket
[Chain grammar and subgroup embedding: how `--chain` names a chain and the default chain of every group](https://github.com/hydeik/pointgroup-tables/issues/42),
on the facts of
[Subgroup catalogue: which group chains of every family are canonical](https://github.com/hydeik/pointgroup-tables/issues/41).
Builds on [Abstract groups](./abstract-groups.md) (symbols, aliases,
canonical element order), [Conjugacy classes](./conjugacy-classes.md) (class
labels), [Irrep labels](./irrep-labels.md) (principal elements) and
[Command line](./cli.md). This chapter fixes **what a chain is, how
`--chain` names one, which concrete subgroups the symbols pick, the
default chain of every group, and the chain-adapted basis a chain selects**
(the basis section was decided in
[The chain-adapted basis: the exact eigenfunction method, component order and the phase rule](https://github.com/hydeik/pointgroup-tables/issues/43);
the embedding order rule and the designated component were amended on the
prototype's findings in
[Embedding order rule: whole-sequence or highest-order-first, and the designated component of T ⊃ D2 ⊃ C2](https://github.com/hydeik/pointgroup-tables/issues/58),
and the canonical spelling and the twin chain classes were settled in
[Chain classes that class anchors cannot separate: omit them from list chains, or a finer anchor](https://github.com/hydeik/pointgroup-tables/issues/59)).
The `subgroup_chain` table and the chain labels of `matrices` are fixed in
[The `subgroup_chain` table and the chain labels](#the-subgroup_chain-table-and-the-chain-labels)
(decided in
[The `chain` table and the chain labels of the `matrices` components](https://github.com/hydeik/pointgroup-tables/issues/45));
the normalisation arithmetic and the conductor of a `--chain` document are
fixed in [Exact computation](#exact-computation) (decided in
[Exact eigenvectors over the cyclotomic field: normalisation and the conductor of `--chain` documents](https://github.com/hydeik/pointgroup-tables/issues/44));
nothing here changes a document produced without `--chain`.

## Chains

A **group chain** of the group `G` of a document is a sequence of point
groups

```
G = H₀ ⊃ H₁ ⊃ … ⊃ Hₘ        (m ≥ 0)
```

in which every `Hᵢ` (`i ≥ 1`) is a proper subgroup of `Hᵢ₋₁`, named by a
Schoenflies symbol and resolved to one concrete subgroup of `G` by the
[embedding rule](#embedding-which-subgroup-a-symbol-picks). The groups
`H₁ … Hₘ` are the **steps** of the chain; `Hₘ` is its **terminal**.

A chain is accepted by `--chain` exactly when all three hold:

1. **Every step is a subgroup** of the previous group, of the point-group
   type its symbol names (section [Embedding](#embedding-which-subgroup-a-symbol-picks)).
2. **Every restriction is multiplicity-free for both kinds of irrep**: for
   every `i` and every ordinary and every projective irrep `Γ` of `Hᵢ₋₁`,
   no irrep of `Hᵢ` (of the same kind, with the restricted factor system)
   occurs more than once in `Γ ↓ Hᵢ`. Such a chain is **canonical**. The
   test does not depend on `--irreps`: a chain that is canonical for
   ordinary irreps only (`I ⊃ C5`, `O ⊃ C3`) or for projective irreps only
   (`T ⊃ C2`) is rejected whatever the document holds, so that the chain,
   its table and the basis of `matrices` never depend on a filter.
3. **The terminal is of the cyclic shape**: `Cn`, `Cnh` or `S2n` for some
   `n` (this includes `C1h`, `S2` and `C2`), the groups whose double group
   is abelian. Every ordinary and projective irrep of the terminal is then
   one-dimensional, so the chain labels determine every component of every
   irrep of `G` uniquely. A chain that stops above the cyclic shape
   (`T ⊃ D2`, `O ⊃ D4 ⊃ D2`, `Oh ⊃ O`) is rejected with a hint naming a
   completion; the printed tables of the literature for such chains are
   covered by a completed chain, because for the irreps the shorter chain
   does label, the completed chain's components coincide with the shorter
   chain's up to phase.

Steps below a cyclic-shape group (`O ⊃ C4 ⊃ C2`, `C6h ⊃ C3`) are accepted:
every such restriction is multiplicity-free, the extra labels are merely
redundant. The empty chain (`m = 0`) is accepted only as the
[default chain](#default-chains) of a cyclic-shape group.

Facts behind these rules, proven for every `n` and computed for the
polyhedral groups in the catalogue of ticket #41 (`docs/research/subgroup-chains.md`
on the branch `research/subgroup-chains`):

- every chain of the cyclic shape is canonical;
- for the dihedral shape `D_N` (`Dn`, `Cnv` with `N = n`; `Dnh` for odd `n`
  and `Dnd` for even `n` with `N = 2n`; `Dnh` for even `n` and `Dnd` for odd
  `n` as `Dn × Ci`) a chain is canonical if and only if its only cyclic step
  is `D_d ⊃ C_d` of the same order, or it ends at a secondary `C2` or a
  reflection: `Dn ⊃ Cn`, `Cnv ⊃ Cn`, `Dnh ⊃ Cnh`, `Dnd ⊃ S2n`,
  `Dn ⊃ C2'`, `Cnv ⊃ C1h` and every `D_N ⊃ D_d ⊃ … ⊃ C_d` descent are
  canonical for every `n`; `D_N ⊃ C_d` with `d < N` never is;
- polyhedral: `T ⊃ C3`, `T ⊃ D2 ⊃ C2`; `O ⊃ C4`, `O ⊃ D4 ⊃ C4`,
  `O ⊃ D3 ⊃ C3`, `O ⊃ T ⊃ C3`, `O ⊃ D4 ⊃ D2 ⊃ C2`, `O ⊃ T ⊃ D2 ⊃ C2`;
  `Td ⊃ S4`, `Td ⊃ C3v ⊃ C3`, `Td ⊃ T ⊃ C3`; `Th ⊃ S6`, `Th ⊃ C3`,
  `Th ⊃ T ⊃ C3`; `Oh ⊃ C4h`, `Oh ⊃ C4`, `Oh ⊃ O ⊃ C4`, `Oh ⊃ D4h ⊃ C4h`;
  `I ⊃ D5 ⊃ C5`, `I ⊃ T ⊃ C3`; `Ih ⊃ D5d ⊃ S10`, `Ih ⊃ I ⊃ D5 ⊃ C5` are
  canonical; no one-step chain of `I` or `Ih` is.

## The `--chain` option

`--chain <chain>` is a common option of every table subcommand
([Command line](./cli.md#options)). Its value is one token:

```
chain  := "default" | step (">" step)*
step   := symbol anchor?
anchor := "[" rep ("," rep)* "]"
```

- `symbol` is a Schoenflies symbol under the
  [group-argument grammar](./cli.md#the-group-argument): the same regular
  expression, the same family keys, the same `n` rules, the same aliases,
  matched case-insensitively and normalised the same way (`cs` → `C1h`,
  `v` → `D2`, `d1h` → `C2v`).
- `rep` is the reduced representative of a
  [class label](./conjugacy-classes.md#class-labels) of `G` (the class
  label without its size prefix: `C2`, `C2'`, `C2''`, `sigma_h`, `sigma_v`,
  `sigma_d`, `C4`, `S4`, `C5^2`, `C3^-1`, …), a **class anchor** restricting
  which subgroup the symbol picks ([Embedding](#embedding-which-subgroup-a-symbol-picks)).
- The group `G` itself is **not written**: `--chain "D4>C4"` on the group
  `O` names `O ⊃ D4 ⊃ C4`. A first step equal to `G` is rejected like any
  symbol that is not a proper subgroup.
- No whitespace anywhere. The separator is `>`; the value is written in
  **double quotes** on the command line (`>` is a redirection and `[` a
  glob in every common shell, and a class anchor may contain `'`):
  `--chain "D4>C2[C2']"`.
- `default` (case-insensitive, no anchor, no steps) selects the group's
  [default chain](#default-chains). It is a reserved word of this option
  only; it is never a group symbol.

`--chain` is absent by default, and then no chain applies: the document is
the version-1 document, in the canonical basis. `--chain` combines with
`--basis-convention canonical` only; with any other basis convention it is
a usage error (exit `2`). It combines freely with `--label-convention`,
which renames the irreps of `G` and nothing else: chain labels are always
canonical Mulliken symbols of the steps.

What `--chain` adds to a document (the `chain` metadata key below, the
`subgroup_chain` table, the chain labels of `matrices` components and the
basis switch of `matrices`) is fixed by the later sections of this chapter
and by [Output schema](./output-schema.md); it only adds keys.

## Embedding: which subgroup a symbol picks

A symbol names an abstract point group; `G` may contain that group several
times, as conjugate subgroups (the three `D4` of `O`, the four `C3v` of
`Td`, the `n` subgroups `C2'` of `Dn`) and as non-conjugate subgroups with
the same symbol (the central `C2` and the `C2'`, `C2''` subgroups of
`D4`; the `C1h` of `D6h` generated by `sigma_h`, `sigma_v` or `sigma_d`;
the normal `D2` and the `D2'` of `O`). The symbol is resolved to **one**
subgroup by the rule below, deterministically; the user steers it only
through class anchors. Subgroups are never named by generators, axis
indices or coordinates ([ADR-0010](../adr/0010-chain-embedding-by-rule.md)).

**Admissible isomorphism.** Let `H_lib` be the group the step's symbol
names, with its own canonical element list and double group `H_lib†`
([Abstract groups](./abstract-groups.md), [Double groups](./double-groups.md)),
and let `K` be the group of the previous step (`G` for the first step). An
**admissible isomorphism** is an isomorphism of double groups
`φ: H_lib† → K†` onto the double-group preimage of a subgroup of `K` that

- preserves the parity of every element;
- preserves, for every element, the invariants of its canonical label:
  the kind (rotation, reflection, rotoreflection, inversion, identity), the
  rotation order and the modulus of the exponent;
- sends, for every element that is not a half turn (order `2` in `G`), the
  lifted preimage to the lifted preimage (half turns have a conventional
  section, their two preimages are not distinguished).

Every element of `K` is an element of `G`, so the image of `φ` is a list
of **canonical indices of `G`**.

**Resolution of a step.** Among all admissible isomorphisms for the step:

1. **Class anchors.** Keep those whose image contains, for every `rep` of
   the anchor, at least one element of the class of `G` that `rep` names.
   An anchor whose `rep` is no class label of `G`, or that leaves no
   isomorphism, rejects the chain (exit `6`).
2. **Canonicity.** Keep those for which the restriction `K ⊃ φ(H_lib)` is
   multiplicity-free for both kinds of irrep. If none remains the chain is
   not canonical (exit `6`).
3. **Order.** Of the remaining isomorphisms take the one whose image,
   written as the sequence of canonical indices of `G` of `φ(h)` for `h`
   running through the canonical element list of `H_lib†` **sorted by
   descending element order** (the order of `h` as an element of
   `H_lib†`; elements of equal order in their canonical index order, so
   the lifted principal element of `H_lib` is compared first and the
   half turns last), is **lexicographically smallest**: the elements of
   highest order decide, the lower ones break the ties.

The step is the subgroup `φ(H_lib)` with that `φ`; the chain labels of the
later sections are the Mulliken symbols of `H_lib` pulled back through
`φ`, so `φ` also fixes the orientation of a subgroup that occurs once:
which element of `K` is the principal element of `H_lib` (`C_n^1` rather
than `C_n^-1`, deciding `1E` against `2E`), and for a dihedral `H_lib`
which half turn is its principal element and which its secondary element
(deciding `B1`, `B2`, `B3` of `D2`).

The descending-order comparison was decided in ticket
[#58](https://github.com/hydeik/pointgroup-tables/issues/58), amending ticket #42, which compared the images
in the canonical order of `H_lib†` itself: the prototype of ticket #50
showed that this let the half turns of `T`'s element list, listed before
its `C3`, pick `T` in `I` through `C2(1)` and orient the unique `T` of
`O`, `Td` and `Th` through a non-principal `C3`, so that the first
consequence below was false and the pulled-back sub-chain of `O ⊃ T ⊃ C3`
was `T ⊃ C3^1(k)`, `k ≠ 1`, a chain class `T` itself does not list.

Consequences of the order rule, read off the canonical element order of
`G` (identity first, then the elements about the principal axis, then the
secondary elements in axis order; polyhedral groups by ascending rotation
order and axis index) and the descending-order comparison:

- a symbol picks the conjugate containing the element of its highest
  order with the **smallest canonical index**: in the polyhedral groups the
  conjugate through the polyhedral principal element (`C4^1(1)` for `D4`,
  `C4`, `D4h`, `C4h`, `C4v`; `C3^1(1)` for `C3`, `D3`, `C3v`, `S6`, `D3d`,
  `T`; `S4^-1(1)` for `S4`, `D2d` in `Td`; `C5^1(1)` for `C5`, `D5`, `C5v`,
  `D5d`, `S10`), in the dihedral shape the subgroup containing `C_N^1`, or
  for `D_d` with `d < N` the one containing the secondary element of
  smallest index;
- the **pulled-back sub-chain** of every step is the step group's own
  chain class through its principal element: `O ⊃ T ⊃ C3`, `Td ⊃ T ⊃ C3`,
  `Th ⊃ T ⊃ C3` and `I ⊃ T ⊃ C3` pull back to `T ⊃ C3^1(1)`, the default
  chain `T --chain C3` itself, and `G ⊃ T ⊃ D2 ⊃ C2` to
  `T --chain "D2>C2"`, so the step-1 blocks of the
  [chain-adapted matrices](#properties) are the subgroup's own tables
  and no norm-factor entry is needed for a sub-chain beyond those of the
  subgroup's `list chains` classes;
- an unanchored `C2` under a dihedral-shape group is a **secondary** `C2`
  (the central `C_N^{N/2}` is dropped by canonicity: `D_N ⊃ C2` through the
  central element has multiplicity two for ordinary irreps), the `C2'`
  class before the `C2''` class: `--chain "C2"` on `D6` is the real
  `Dn ⊃ C2x` chain of Chen, Fan, McAven and Butler;
- an unanchored `C1h` under `Cnv` is a `sigma_v` plane, under `Dnh` or `Dnd`
  a `sigma_v` or `sigma_d` plane (`sigma_h` is dropped by canonicity);
- `--chain "D4>C2"` on `O` picks a `C2` perpendicular to the `C4` axis of
  the `D4` (the `C2` on that axis is dropped by canonicity), `C2(2)` before
  `C2'(5)`;
- the two classes of `D_d` in `D_N` at even `N/d` (one through `C2'(1)`,
  one through `C2''(1)`) are both canonical with isomorphic branching
  (`B1 ↔ B2` exchanged); the unanchored symbol picks the `C2'` class, the
  anchor `[C2'']` the other: `--chain "D3[C2'']>C3"` on `D6`;
- `D2` under `O` is the normal `D2` (three `C2`); `O ⊃ D2` is not canonical,
  and `--chain "D2[C2']>C2"` (the `D2'` of the literature, `⟨C2, C2'⟩`) is
  canonical for ordinary irreps only, so both are rejected.

**Fully qualified spelling.** The **class-label signature** of a subgroup
of `G` is the set of class labels of `G` whose classes it meets: everything
a class anchor can see. The **canonical spelling** of an accepted chain is
the spelling in which a step carries a class anchor exactly when another
accepted chain class with the same symbols, and the same signatures at the
steps before, has a different signature at that step, and the anchor is the
**shortest list of reps that excludes every such class**: the fewest reps,
and among lists of that length the first in canonical class order
(`C2v[C2,sigma_h]` under `D4h` in `Oh`, where neither `C2` nor `sigma_h`
alone excludes the `C2v` on a `C2'` axis or the one on `sigma_d` planes).
When no list excludes them all, because the step's signature is contained
in every competitor's, the step stays **bare** and the order rule reaches it:
under `D4` in `O` the normal `D2` meets only the class `C2` and the `D2'`
meets `C2` and `C2'`, so the chains are spelled `"D4>D2>C2"` and
`"D4>D2[C2']>C2[C2]"`. Conjugates within one class are never
distinguished, and the order rule names the representative; a spelling
names one chain class except for the twins below. `list chains` prints
canonical spellings and the `chain` metadata records one, so that a chain
can be reproduced independently of the unanchored defaults; every canonical
spelling resolves to the class it was derived from (ledger row `GC-4`).

**Twin chain classes.** A canonical spelling does not always name one
chain class. Two chain classes that are not conjugate in `G` can have the
same class-label signature at every step, **twin chain classes**: the
elements that would tell them apart lie in one class of `G`, and no element
of `G` that normalises the steps above swaps them (the three `C2` of the
normal `D2` of `O` are one class, but the `D4` above fixes the one on its
axis). Class anchors are class labels of `G` ([ADR-0010](../adr/0010-chain-embedding-by-rule.md)),
so no `--chain` value separates twins: the embedding rule resolves the
shared spelling to one of them, the **reached** twin, and the others are
**not accepted** by `--chain`. `list chains` holds one entry per canonical
spelling, the reached chain, and the `embedding` of the `subgroup_chain`
table shows which subgroup it holds. Decided in ticket
[#59](https://github.com/hydeik/pointgroup-tables/issues/59) against a finer
anchor (a class label of the step above, or an element label) and against
rejecting the shared spelling, because twins differ **only on the
projective irreps**: in every twin pair found, the step at which the twins
differ lies under a group whose ordinary irreps are all one-dimensional
(`D2`, `C2v`, `D2h`) or shares a reflection with the reached subgroup, so
the label spaces of every ordinary irrep coincide (the ordinary blocks are
the same vectors, up to the labels of that step), while the projective label
spaces differ (the `E1/2` of `D2†` is split by the eigenvalue of the `C2` on
the chain's axis in the reached twin and by that of a perpendicular `C2` in
the other). The reached twin is the one through the element of smaller
canonical index: the half turn on the chain's own axis, the first plane.
The twins of the seven polyhedral groups, computed by the prototype of
ticket #50 (`prototypes/chain_prototype_summary.md`, section "Twin chain
classes"); the axial families at `n ≤ 12` have none:

| Group | Spelling | Step | Reached | Unreached |
| ----- | -------- | ---- | ------- | --------- |
| `O` | `D4>D2>C2` | 3 `C2` | `C2(1)`, on the `C4` axis of the `D4` | `C2(2)`, perpendicular to it |
| `Td` | `D2d>D2>C2` | 3 `C2` | `C2(1)`, on the `S4` axis of the `D2d` | `C2(2)` |
| `Th` | `C2v>C1h`, `D2h>C2v>C1h` | last `C1h` | `sigma_h(2)`, the first plane of the `C2v` | `sigma_h(3)` |
| `Oh` | `D4>D2>C2`, `D2d[sigma_d]>D2>C2`, `D4h>D2>C2`, `D4h>D2h>C2`, `D4h>D4>D2>C2`, `D4h>D2h>D2>C2`, `D4h>D2d[sigma_d]>D2>C2`, `O>D4>D2>C2`, `Td>D2d>D2>C2` | last `C2` | `C2(1)`, on the axis of the chain | `C2(2)` |
| `Oh` | `D4h>D2h>C1h` | 3 `C1h` | `sigma_h(1)`, perpendicular to the `D4h` axis | `sigma_h(2)` |
| `Oh` | `D4h>D2h>C2h` | 3 `C2h` | `{C2(1), i, sigma_h(1)}` | `{C2(2), i, sigma_h(2)}` |
| `Oh` | `D4h>C2v[C2,sigma_h]>C1h`, `D4h>D2h>C2v>C1h` | the `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}`, about the `D4h` axis | `{C2(2), sigma_h(1), sigma_h(3)}`, about a perpendicular `C4` axis, with either of its planes as the `C1h` (two unreached classes) |
| `Oh` | `D4h>C2v[C2,sigma_h]>C2`, `D4h>D2h>C2v>C2` | the `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}` | `{C2(2), sigma_h(1), sigma_h(3)}` |
| `Ih` | `Th>C2v>C1h`, `Th>D2h>C2v>C1h` | last `C1h` | `sigma(11)`, the first plane of the `C2v` | `sigma(14)` |

Twenty-one spellings cover forty-four chain classes, twenty-three of them
unreached (`Oh`: 15 spellings, 32 classes; `list chains Oh` holds 143 of
the 160 accepted classes); every other accepted chain class of these groups
has its own spelling.

## Default chains

Every group has one default chain, selected by `--chain default`. It
descends to the **maximal cyclic-shape subgroup of `G` that contains the
principal element** (the [principal element](./irrep-labels.md#named-elements)
of an axial group, the polyhedral principal element of a polyhedral group)
by the shortest canonical route:

| Group | Default chain | Terminal contains |
| ----- | ------------- | ----------------- |
| `Cn`, `Cnh`, `S2n` | empty (`G` alone, no step) | `G` is its own terminal |
| `Dn`, `Cnv` | `Cn` | `C_n` |
| `Dnh` (every `n`) | `Cnh` | `C_n` (even `n`), `S_n^-1` (odd `n`) |
| `Dnd` (every `n`) | `S2n` | `S_2n` (even `n`), `C_n` (odd `n`) |
| `T` | `C3` | `C3^1(1)` |
| `Th` | `S6` | `C3^1(1)` |
| `O` | `C4` | `C4^1(1)` |
| `Oh` | `C4h` | `C4^1(1)` |
| `Td` | `S4` | `S4^-1(1)` |
| `I` | `D5>C5` | `C5^1(1)` |
| `Ih` | `D5d>S10` | `C5^1(1)` |

The axial rows are one step for every `n` and coincide with the chains of
Chen, Fan, McAven and Butler (2000); `T ⊃ C3`, `O ⊃ C4` and `I ⊃ D5 ⊃ C5`
are the chains of the Chen–Fan papers. `I` and `Ih` need two steps because
`I ⊃ C5` is not canonical (`I5/2 ↓ C5†` contains `B5/2` twice). The table is
normative; the rule explains it. The empty default chain of a cyclic-shape
group adds the `chain` metadata and the `subgroup_chain` table (with no
subgroup record) and changes no matrix: the canonical basis of a
one-dimensional irrep is already adapted to every chain.

## Chain-adapted basis

Decided in
[The chain-adapted basis: the exact eigenfunction method, component order and the phase rule](https://github.com/hydeik/pointgroup-tables/issues/43),
on the research notes of
[Survey of the eigenfunction method in the Chen–Fan–Ping papers](https://github.com/hydeik/pointgroup-tables/issues/39)
and
[Cyclotomic closure of chain-adapted bases and Clebsch–Gordan coefficients](https://github.com/hydeik/pointgroup-tables/issues/40).
Rationale in [ADR-0011](../adr/0011-chain-basis-phase-rule.md). This
section defines, for one irrep of `G` and one resolved chain, the basis in
which `matrices` writes the irrep under `--chain`: it is a **change of
basis of the version-1 canonical realisation**, computed exactly, and it
touches neither the characters nor the irrep records of the document.

### Setting

- `Γ` is an irrep of `G` of dimension `d`, ordinary or projective, in its
  canonical realisation `D_Γ` ([Basis conventions](./basis-conventions.md)):
  a representation of `G†` with `D_Γ(ḡ) = −D_Γ(g)` for a projective irrep
  and `D_Γ(ḡ) = D_Γ(g)` for an ordinary one. `e_1, …, e_d` are the canonical
  basis vectors in canonical component order.
- A **resolved chain** is a sequence of subgroups `H_1 ⊃ … ⊃ H_m` of `G`
  with the admissible isomorphisms `φ_i: H_i,lib† → G†` of the
  [embedding rule](#embedding-which-subgroup-a-symbol-picks). The
  definition below applies to any resolved chain, not only to one named by
  `--chain`; it is used recursively on the chains it pulls back.
- The irreps of `H_i,lib` that can label a component of `Γ` are those **of
  the kind of `Γ`**: ordinary for an ordinary `Γ`, projective (with the
  factor system of `H_i,lib†`, which `φ_i` carries to the restriction of
  `G†`'s) for a projective `Γ`. `σ` runs over them in the
  [canonical irrep order](./irrep-labels.md#canonical-irrep-order) of
  `H_i,lib`.

### Label spaces

For every step `i` and every irrep `σ` of `H_i,lib` of the kind of `Γ`, the
**isotypic projector** of `σ` in `Γ` is

```
P^σ_i = (dim σ / |H_i,lib†|) · Σ_{h ∈ H_i,lib†} conj(χ_σ(h)) · D_Γ(φ_i(h)),
```

the sum running over both preimages of every element of `H_i` (so that for
an ordinary `Γ` it equals the sum over `H_i` itself). The `P^σ_i` are
Hermitian idempotents with entries in `Q(ζ_lcm(M_G, M_{H_i}))`, and
projectors of different steps commute. For a sequence `s = (σ_1, …, σ_m)` of
one irrep per step the **label space** is

```
V_s = range(P^{σ_1}_1 · P^{σ_2}_2 · … · P^{σ_m}_m).
```

Because the chain is canonical for both kinds and ends at the cyclic shape
([Chains](#chains)), every `V_s` has dimension `0` or `1` and the non-zero
ones sum to `Γ`. The **chain labels** of `Γ` are the sequences `s` with
`dim V_s = 1`; there are exactly `d` of them. An implementation checks both
facts for every irrep (ranks `≤ 1`, ranks summing to `d`) and treats a
failure as an internal error, never as a usage error: a chain accepted by
the grammar section cannot fail here.

A one-dimensional `Γ` has one chain label, its own restriction, and its
chain-adapted matrix is its character. The empty chain (the default chain of
a cyclic-shape group) selects the canonical basis unchanged.

### Component order

The components of `Γ` in the chain-adapted basis are its chain labels in
**chain order**: by `σ_1` first, then `σ_2`, and so on (lexicographic in
the steps), the labels of one step compared by the **step order** of that
step's group:

- a step of the **dihedral or polyhedral shape** (`D_d`, `T` and their
  centro and twisted relatives) orders its labels by its
  [canonical irrep order](./irrep-labels.md#canonical-irrep-order);
- a step of the **cyclic shape** (`Cn`, `Cnh`, `S2n`, including `C2` and
  `C1h`) orders its labels in **exponent order**: by the exponent `e` of
  the label's eigenvalue at the step group's own principal element `c`
  ([Irrep labels](./irrep-labels.md#named-elements): `C_n` for `Cn`, for
  `Cnh` with even `n` and for `S2n` with odd `n`; `S_2n` for `S2n` with even
  `n`; `S_n^-1` for `Cnh` with odd `n`, so `C_2` for a `C2` step and `σ_h`
  for a `C1h` step, whose embedded images are secondary elements of `G`),
  `χ_σ(c) = ζ_N^e`
  with `e ∈ (−N/2, N/2]` for an ordinary label and `χ_σ(s(c)) = ζ_2N^e`
  with `e ∈ (−N, N]` for a projective label, `N` the order of `c`,
  **descending**; labels with equal exponents (only `Cnh` with even `n`,
  where `Ag` and `Au` share `e = 0`, and `σ_h` is central so that the two
  never label components of one irrep) by canonical irrep order. This is
  the rule by which [Basis conventions](./basis-conventions.md) orders the
  canonical components of version 1, applied to the step's own principal
  element.

Components with the same `σ_1` are therefore contiguous, which is what
makes the restriction to `H_1` block-diagonal with one block per label,
and `1E` precedes `2E` at every cyclic step as in version 1. Examples:
the `C4` step of `O` orders `B, 1E, A, 2E` (`e = 2, 1, 0, −1`), so `T1` has
the components `(1E, A, 2E)` and `E` the components `(B, A)`, as the
canonical basis lists them; the `C3` step of `T` orders `1E, A, 2E`; the
`Cnh` step of an odd `Dnh` orders by the exponent at `S_n^-1`, the
canonical order of `Dnh` itself; a secondary `C2` terminal orders
`B, A` (`e = 1, 0`) and a `C1h` terminal `A'', A'`.

Decided in
[Component order at cyclic-shape steps: subgroup irrep order, or the exponent order that makes every one-step default chain reproduce version 1 entry by entry](https://github.com/hydeik/pointgroup-tables/issues/57)
on the computation of
[Does the version-1 canonical basis coincide with the chain-adapted basis through the principal element](https://github.com/hydeik/pointgroup-tables/issues/49),
amending the first rule of ticket #43 (every step by the subgroup's
canonical irrep order). That rule put `A` before `1E` at the `C4` and `C3`
steps of the polyhedral groups and, because the `1E`/`2E` labels of an odd
`Cnh` are assigned at `c₀ = C_n` while the canonical order of `Dnh` uses
`c = S_n^-1`, swapped the two components of `E_k'` (`k` odd) and `E_k''`
(`k` even) of every odd `Dnh`; the exponent order makes every one-step
chain through the principal element reproduce the canonical component
order ([Properties](#properties)).

### The norm-factor phase rule

Every phase of the basis is fixed by one primitive. A rule **designates**
one non-zero quantity `x` (a vector component or a matrix element) whose
squared modulus `N = x · conj(x)` is a totally positive element of the real
subfield of the document field `F`; the phase is then chosen so that

```
x = conj(μ_F(N)),
```

where `μ_F(N) ∈ F` is the **canonical norm factor** of `N`: an element with
`μ · conj(μ) = N`, equal to the real positive square root `√N` whenever that
root lies in `F`, and otherwise fixed by
[Exact computation](#the-canonical-norm-factor), which also fixes the
conductor of `F` so that the factor exists for every `N` a chain produces. Where `√N ∈ F` the rule is
therefore the familiar "designated quantity real and positive"; where it is
not (`N = (3 ± √3)/6` under `T ⊃ D2 ⊃ C2`, whose square roots lie in no
cyclotomic field), the designated quantity is a canonical element of modulus
`√N` instead, and every value stays in `F`
([Cyclotomic closure](https://github.com/hydeik/pointgroup-tables/issues/40)).
The same primitive fixes the phases of the Clebsch–Gordan coefficients
([Clebsch–Gordan coefficients](./clebsch-gordan.md)).

### The basis

The chain-adapted basis `(w_s)` of `Γ`, indexed by its chain labels in
component order, is defined by recursion on the length `m` of the chain:

0. **`m = 0`**: `w_s = e_s`, the canonical basis in canonical order.
1. **Blocks.** Group the chain labels by their first entry: the **block**
   `B_σ = { s : σ_1 = σ }` for every `σ` that occurs. Realise `σ` in
   `H_1,lib` adapted to the **pulled-back chain**
   `H_1,lib ⊃ φ_1⁻¹(H_2) ⊃ … ⊃ φ_1⁻¹(H_m)` (embeddings `φ_1⁻¹ ∘ φ_i`), by
   this same definition, in the canonical realisation of `σ` from the
   version-1 tables of `H_1,lib`: the result is the matrix family `D'_σ(h)`
   with components indexed by `(σ_2, …, σ_m)`, in the same order as `B_σ`.
   For a one-dimensional `σ` it is the character of `σ`.
2. **Within a block.** Choose unit vectors `w_s ∈ V_s`, `s ∈ B_σ`, such that

   ```
   ⟨w_s | D_Γ(φ_1(h)) | w_{s'}⟩ = D'_σ(h)_{s s'}      for every h ∈ H_1,lib†.
   ```

   Such a family exists because the restriction of `Γ` to `V_σ = ⊕ V_s` is
   equivalent to `σ`, and it is unique up to one common factor of modulus
   one (Schur's lemma: `σ` occurs once). Concretely, `w_t` for the first
   label `t` of the block fixes every other `w_s` through the entries
   `D'_σ(h)_{s t} ≠ 0`.
3. **Between blocks.** Fix the common factor of every block by its first
   label `t` and its **designated component** `e_j`: among the canonical
   indices with `(P_t)_{jj} ≠ 0` (`P_t` the projector onto `V_t`;
   equivalently `⟨e_j | w_t⟩ ≠ 0`), `j` is the smallest one whose
   canonical basis vector is fixed by the principal element `c` of `G`
   (`D_Γ(c) e_j = e_j`: the canonical basis is adapted to `c`, so these
   are the components of exponent `0`), and the smallest one of all when
   the block has no such index (every projective irrep; an ordinary irrep
   whose block misses the exponent-`0` component). Apply the norm-factor
   rule to `x = ⟨e_j | w_t⟩`. Here `N = x · conj(x) = (P_t)_{jj}`, so
   `w_t = P_t e_j / μ_F((P_t)_{jj})` is the explicit form. The preference
   for the fixed component was decided in ticket [#58](https://github.com/hydeik/pointgroup-tables/issues/58)
   (amending ticket #43, which took the smallest index outright): where
   `c` normalises `H_1` it makes `D'_Γ(c)` permute the blocks without
   phases ([Properties](#properties)); it never touches a chain through
   the principal element, where every block is a canonical vector.

The **chain-adapted matrices** are `D'_Γ(g)_{s s'} = ⟨w_s | D_Γ(g) | w_{s'}⟩`
for every `g` of the element list (`G`, or `G†` through the canonical
section), with entries in `F`. Only the designated components need a norm
factor: the other vectors of a block are reached from `w_t` by matrices over
`F`, so `F` must contain the fields of the pulled-back chains (the conductor
rule of the arithmetic section).

No further convention is imposed. In particular the tool does not require
complex-conjugate partner irreps to have conjugate bases, does not reproduce
any paper's phases entry by entry, and does not follow the Condon–Shortley
convention; the papers' matrices are reference data compared with ours up to
a permutation of components and a diagonal unitary
(`REF-10` of the [check ledger](./verification.md#check-ledger), decided in
ticket [#48](https://github.com/hydeik/pointgroup-tables/issues/48)).

### Properties

An implementation must exhibit these consequences of the definition; they
are the checks `CB-1` to `CB-5` of the
[check ledger](./verification.md#check-ledger) for chain-adapted bases
(`CB-1` intertwiner and ranks, `CB-2` block equality with the sub-chain
realisation, `CB-3` the phase rule, `CB-4` chains through the principal
element, `CB-5` centro and twisted inheritance; `OS-8` the block
structure). `--verify` runs `CB-2`, `OS-7` and `OS-8` on a `--chain` document and
skips `MX-5`; `CB-1` and `CB-3` need the basis matrices, which the
document does not carry, and run in the test suite only:

- **Block structure.** For every step `i`, `D'_Γ(φ_i(h))`, `h ∈ H_i,lib†`,
  is block-diagonal in the components grouped by `(σ_1, …, σ_i)`; for
  `i = 1` the block of `σ` **equals** `D'_σ` of step 1 entry by entry; for
  the terminal, `D'_Γ(φ_m(h))` is diagonal with the characters of the
  terminal labels.
- **The principal element across blocks.** Where the principal element
  `c` of `G` normalises `H_1` without lying in it, it permutes the
  one-dimensional blocks of an ordinary irrep (`T ⊃ D2`, `Th ⊃ D2`,
  `Th ⊃ D2h`, `Th ⊃ C2v` is not normalised), and the designated
  component rule makes `D'_Γ(c)` carry the first (and only) vector of a
  block designated at a `c`-fixed component `e_j` onto the vector of the
  image block with phase `1`: `D_Γ(c) P_t e_j = P_{c·t} D_Γ(c) e_j =
  P_{c·t} e_j`, and both blocks share `N = (P_t)_{jj}`. On `T` of
  `T ⊃ D2 ⊃ C2` (every block holds the exponent-`0` component, `N = 1/3`)
  `D'_T(C3^1(1))` is therefore the real cyclic permutation of
  [CF98b] eq 17, and by block equality so is the `T` block of every
  `G ⊃ T ⊃ D2 ⊃ C2`. Projective irreps have no fixed component and keep
  the first-index rule; where `c` does not normalise `H_1` (`O ⊃ D3`,
  `I ⊃ T`) the rule only chooses the designated squared norm (`1/3`
  instead of `(2 − √3)/6` on `T1` of `O ⊃ D3 ⊃ C3`, `1/6` instead of
  `5/12` on `T1` of `O ⊃ D3 ⊃ C2`). Checked through `CB-3` (the rule) and
  reported by `REF-10` (the phases against the papers); no ledger row of
  its own (ticket #58).
- **Chains through the principal element.** For a one-step chain whose
  subgroup is the cyclic-shape group of the principal element (`Dn ⊃ Cn`,
  `Cnv ⊃ Cn`, `Dnh ⊃ Cnh`, `Dnd ⊃ S2n`, `T ⊃ C3`, `O ⊃ C4`, `Td ⊃ S4`,
  `Th ⊃ S6`, `Oh ⊃ C4h`: every default chain but those of `I` and `Ih`) the
  chain-adapted basis **is the canonical basis**: every `V_s` is spanned by
  one `e_k` (the eigenvalues of the principal element are distinct, with
  the single exception `I5/2` of `I†`, which no one-step chain of `I`
  reaches), so the designated component is `e_k` itself with `N = 1` and
  the vector is `e_k` with phase `1`; and the exponent order of the step is
  the canonical component order, since every `V_s` is an eigenvector of
  `D_Γ(c)` with eigenvalue `χ_σ(c)`. `matrices` under the default chain
  therefore **equals the canonical table entry by entry, with chain labels
  added**, for every axial family and for `T`, `Th`, `Td`, `O`, `Oh`
  (computed for every axial family at `n ≤ 12` and the five polyhedral
  groups in ticket [#49](https://github.com/hydeik/pointgroup-tables/issues/49), order rule decided in ticket
  [#57](https://github.com/hydeik/pointgroup-tables/issues/57)).
  Longer chains through the principal element span the same label spaces
  on every irrep but `I5/2`: `O ⊃ D4 ⊃ C4`, `Oh ⊃ D4h ⊃ C4h` and
  `Td ⊃ D2d ⊃ S4` differ from the canonical vectors by the diagonal unitary
  `diag(1, −i)` on every `E1/2` block and `diag(1, i)` on every `E3/2` block
  of the spinor irreps (ordinary irreps: phase `1`; every `N = 1`), forced
  by the within-block rule because the canonical `D4` writes `D(C2'(1))`
  of `E1/2` as `[[0, i], [i, 0]]` where the canonical `O` writes `D(C2(2))`
  as `[[0, −1], [1, 0]]`; `I ⊃ D5 ⊃ C5` and `Ih ⊃ D5d ⊃ S10`, the default
  chains of `I` and `Ih`, differ by a phase on the second component of
  every `E1`, `E2`, `E1/2`, `E3/2` block of `D5` (`−ζ_20^4` on `E1 ⊂ T1`,
  `−ζ_20^8` on `E2 ⊂ T2`, `−ζ_20^3` on `E1/2 ⊂ E1/2`, …; `N = 1`), and on
  `I5/2` the two `−1` eigenvectors of `C5` (canonical components `0` and
  `5`) are split by `D5` into `w = (1 − i)/2 · (e_0 ∓ ζ_20 e_5)` for
  `1E5/2`, `2E5/2`, with `N = 1/2` and `μ = (1 + i)/2` (the `(20, 1/2)` row
  of the norm-factor data). No chain of `I` or `Ih` reproduces the
  canonical basis: the within-block rule ties the `D5` blocks to the
  canonical `D5` tables, which the canonical `I` matrices do not restrict
  to.
- **Centro groups.** For `G = R × Ci` and a chain of subgroups
  `H_i = (S_i, S_i⁰)` ([Subgroup catalogue](https://github.com/hydeik/pointgroup-tables/issues/41)),
  the projectors of `Γ_g` and `Γ_u` on `H_i` equal those of `Γ` on `S_i`,
  so the chain-adapted basis of `Γ_g` and of `Γ_u` is that of `Γ` for the
  rotation chain `R ⊃ S_1 ⊃ …`, and `D'(i·r) = ±D'(r)`.
- **Twisted groups.** A type III group is its rotation group with parities
  attached, so its chain-adapted matrices are those of the rotation chain on
  the image elements (`D3h ⊃ C3h` from `D6 ⊃ C6`, `Td ⊃ S4` from `O ⊃ C4`).
- **Determinism.** The basis depends on the abstract group, the resolved
  chain, the canonical realisation of `G` and of every `H_i,lib`, and the
  canonical norm factor; the same document is emitted for the same
  arguments and tool version.

## Exact computation

Decided in wayfinder ticket
[Exact eigenvectors over the cyclotomic field: normalisation and the conductor of `--chain` documents](https://github.com/hydeik/pointgroup-tables/issues/44),
on the research note of
[Cyclotomic closure of chain-adapted bases and Clebsch–Gordan coefficients](https://github.com/hydeik/pointgroup-tables/issues/40)
and on an exact survey of candidate norm factors run for the ticket.
Rationale in [ADR-0013](../adr/0013-chain-conductor-and-norm-factor.md).
This section fixes **how the chain-adapted basis is computed exactly, the
exact real square root, the canonical norm factor `μ_F(N)`, the conductor
of a `--chain` document, the error raised when a value cannot be placed in
the field, and the cost**. Every value of a `--chain` document is an
element of one cyclotomic field, as in version 1 ([Values](./values.md));
nothing here extends the value model.

### What is computed

The definition of [the basis](#the-basis) needs no eigenvalue problem, no
nullspace and no Gram–Schmidt step: every basis vector is the image of a
canonical basis vector under exact projectors.

1. **Chain field.** The isotypic projectors `P^σ_i` have entries in the
   **chain field** `Q(ζ_L)`, `L = lcm(M_G, M_{H_1}, …, M_{H_m})`, with
   `M_G` and `M_{H_i}` the version-1 conductors of `G` and of the step
   groups `H_i,lib`. Every projector, rank and transfer computation below
   runs in `Q(ζ_L)` with the arithmetic of version 1 (reduced power basis,
   rational coefficients); only the normalisation of step 3 runs in the
   document field `F = Q(ζ_{M_chain})`, `L | M_chain`, and every value is
   lifted to `F` at the output boundary (`ζ_L = ζ_M^{M/L}`), as version 1
   lifts values computed in a smaller field.
2. **Label spaces.** `P_s = P^{σ_1}_1 · … · P^{σ_m}_m` for every sequence
   `s`; its rank is computed by exact Gaussian elimination over `Q(ζ_L)`
   and checked to be `0` or `1`, the ranks summing to `d`
   ([Label spaces](#label-spaces)). This exact rank is the only
   linear-algebra primitive beyond matrix products: no linear system is
   solved.
3. **First vector of a block.** For the first label `t` of a block with
   designated canonical index `j` (the first `j` with `(P_t)_{jj} ≠ 0`
   among the components fixed by the principal element, else the first of
   all, [The basis](#the-basis)),

   ```
   w_t = P_t e_j / μ_F(N),      N = (P_t)_{jj}:
   ```

   the unnormalised vector is the `j`-th column of `P_t`, its squared norm
   is that column's diagonal entry (`P_t` is a Hermitian idempotent), and
   `N` is a totally positive element of the real subfield of `Q(ζ_L)`.
4. **Other vectors of a block.** `w_s = P^σ_{s t} w_t` with the
   matrix-element projector of the step-1 group,

   ```
   P^σ_{s t} = (dim σ / |H_1,lib†|) · Σ_{h ∈ H_1,lib†} conj(D'_σ(h)_{s t}) · D_Γ(φ_1(h)),
   ```

   which carries `w_t` to the unique unit vector `w_s ∈ V_s` satisfying the
   [within-block condition](#the-basis) (the orthogonality relations of
   `D'_σ`); `⟨w_s|w_s⟩ = 1` exactly, no further normalisation is needed.
   `D'_σ` is the chain-adapted realisation of `σ` for the pulled-back
   chain, computed by these same steps in the canonical realisation of
   `H_1,lib`; its field is the chain conductor of that sub-chain, which
   `M_chain` contains by the
   [conductor rule](#the-conductor-of-a-chain-document).
5. **Matrices.** `D'_Γ(g)_{s s'} = ⟨w_s | D_Γ(g) | w_{s'}⟩`, one exact
   product per element of the element list.

The eigenfunction method of the papers (common eigenvectors of the class
operators of the chain) and this projector form select the same label
spaces; the projector form is used because the characters have already
diagonalised it and no root finding is involved.

### Exact real square root

The arithmetic gains one decision procedure, **`real_sqrt_F(a)`**: for an
element `a` of the real subfield `F⁺` of `F = Q(ζ_M)`, the positive real
square root `√a ∈ F⁺`, or the proof that `a` is not a square in `F`. It
never approximates:

- `a` is checked to be totally positive (real, and positive in every
  embedding); a totally positive `a` has a square root in `F` if and only
  if it has a real one;
- the denominator `D` of `a` is cleared: `A = D² a` is an algebraic
  integer, hence so is its root `x = √A`, hence `x` has integer
  coefficients over the power basis (`Z[ζ_M]` is the ring of integers);
- a prime `p ≡ 1 (mod M)` with `p > 2B` is chosen, `B` a bound on the
  absolute value of the coefficients of `x` (from `|σ_k(x)|² = σ_k(A)` and
  the inverse Vandermonde matrix of the `M`-th roots of unity, both bounded
  rationally); modulo `p` the field splits completely, `ζ ↦ r^k` for a
  primitive `M`-th root of unity `r` modulo `p`, so each `σ_k(x)` is a
  square root of `σ_k(A)` modulo `p` and the only freedom is one sign per
  pair `{k, −k}` (a real `x` has `σ_k(x) = σ_{−k}(x)`); for every sign
  pattern (`2^{φ(M)/2 − 1}`, the identity embedding positive) the
  coefficients are reconstructed modulo `p` by the inverse Vandermonde,
  lifted to symmetric residues in `(−p/2, p/2)`, and the candidate is
  **verified exactly** by `x · x = A`; the first verified pattern gives
  `√a = x / D`, with `x > 0` in the identity embedding;
- no pattern verifying proves that `a` is not a square in `F`: a root
  would have produced its own pattern.

Cost control: the root is taken in the smallest cyclotomic subfield
`Q(ζ_m) ⊆ F` containing `a` and lifted; when it does not exist there it is
sought along the quadratic steps of the tower up to `F`
(`Q(ζ_120)⁺ = Q(ζ_60)⁺(√2)`: `√a ∈ Q(ζ_120)` iff `a` or `2a` is a square
in `Q(ζ_60)`), so the sign search never exceeds `2^7` patterns for the
conductors of this specification. Floating-point arithmetic is used
nowhere in the decision: the floating-point reconstruction of the research
script was found to miss genuine squares with large denominators.

### The canonical norm factor

For a totally positive `N` in the real subfield of the document field
`F = Q(ζ_M)`, the **canonical norm factor** is

```
μ_F(N) = √N                         when N is a square in F (real_sqrt_F),
μ_F(N) = (1 + i)/2                  when N = 1/2 and √2 ∉ F   (i = ζ_M^{M/4}, 4 | M),
μ_F(N) = the shipped entry (M, N) ↦ μ    otherwise,
```

and the designated quantity of the
[norm-factor rule](#the-norm-factor-phase-rule) is `conj(μ_F(N))`. The
factor depends on `M` and `N` alone, so the basis is reproducible from the
document field. The three cases are exhaustive for every `N` the tool
meets:

- **Axial families, every `n`.** Every squared norm is `1` or `1/2`
  (proven for all `n` in the cyclotomic-closure note, §2.6): `1` for a
  label space spanned by a canonical basis vector, `1/2` for the two
  labels of an `E_k` or `E_{k/2}` split by a secondary `C2` or a reflection
  (`Dn ⊃ C2'`, `Cnv ⊃ C1h`, the `D_N ⊃ D_d ⊃ …` descents). Under the
  [conductor rule](#the-conductor-of-a-chain-document) the field is
  `Q(ζ_{M_G})` with `4 | M_G` for every group of the dihedral shape, so
  `√2 ∈ F` iff `8 | M_G` (`D4h`, `D8h`, `C8v`, …) and otherwise
  `μ = (1 + i)/2`, designated component `(1 − i)/2`. The cyclic shape has
  `N = 1` only. The chain-adapted **matrices** of the axial families are
  the same under either case, because the two labels of a split `E_k`
  share one `N` and the factor cancels in `⟨w_s|D(g)|w_{s'}⟩`; only the
  basis vectors of the `subgroup_chain` table show the phase.
- **Polyhedral groups.** The squared norms of every irrep of every chain
  class form a finite set; those that are not squares in the declared
  field are **shipped data**, a table `(M, N) ↦ μ` with `μ · conj(μ) = N`
  verified exactly when the table is loaded and by the test suite, in the
  spirit of the frozen polyhedral bases of
  [ADR-0004](../adr/0004-canonical-basis.md). The entries for the chains
  of `T` and `O` (and of `Td`, `Th`, `Oh`, whose chains have the same
  norms), for `I ⊃ D5 ⊃ C5` and for the Clebsch–Gordan coefficients of `I`
  and `Ih` are these, written with the symbols of the
  [value grammar](./basis-conventions.md#value-grammar):

  | `M` | `N` | where it occurs | `μ_F(N)` |
  | --- | --- | --------------- | -------- |
  | 24 | `(3 - sqrt3)/6` | `T ⊃ D2 ⊃ C2`, `O ⊃ T ⊃ C3`, `O ⊃ D3 ⊃ C3` (`E1/2`) | `(3 - sqrt3 + i*sqrt6)/6` |
  | 24 | `(3 + sqrt3)/6` | the same | `(3 + sqrt3 + i*sqrt6)/6` |
  | 24 | `(3 - sqrt3)/12` | `O ⊃ T ⊃ C3`, `O ⊃ D3 ⊃ C3` (`G3/2`) | `(3*sqrt2 - sqrt6 + 2*i*sqrt3)/12` |
  | 24 | `(3 + sqrt3)/12` | the same | `(3*sqrt2 + sqrt6 + 2*i*sqrt3)/12` |
  | 20 | `1/2` | `I ⊃ D5 ⊃ C5` (`I5/2`) | `(1 + i)/2` (the closed form above) |
  | 120 | `1/7` | `cg`: `T1 ⊗ I5/2 ⊃ I5/2` | `(2 + i*sqrt3)/7` |
  | 120 | `2/7` | `cg`: `I5/2 ⊗ I5/2 ⊃ T1` | `(2*sqrt2 + i*sqrt6)/7` |
  | 120 | `7/10` | `cg`: `T1 ⊗ I5/2 ⊃ I5/2` | `(2*sqrt10 + i*sqrt30)/10` |
  | 120 | `7/20` | `cg`: `I5/2 ⊗ I5/2 ⊃ T1` | `(2*sqrt5 + i*sqrt15)/10` |

  Every other squared norm of these chains and pairs is a square in the
  declared field (`1`, `1/2`, `1/3`, `1/4`, `1/6`, `2/3`, `3/8`, `3/4`,
  `(2 ± sqrt3)/6`, … with `sqrt2`, `sqrt3`, `sqrt6` in `Q(ζ_24)`), so the
  designated quantity is real and positive there. The entries for the
  chains of `I` through `T` (`I ⊃ T ⊃ C3`: nine values of degree up to
  eight in `Q(ζ_60)⁺`, factors with up to fifteen non-zero coefficients;
  `I ⊃ T ⊃ D2 ⊃ C2`: five values in `Q(ζ_120)⁺`, factors with up to
  twenty coefficients; the `Ih` counterparts share them) and for every
  chain class of `list chains` not named above are **generated by the
  prototype ticket [#50](https://github.com/hydeik/pointgroup-tables/issues/50)**
  with the procedure below and frozen in the implementation as a data
  file (`prototypes/chain_prototype_norm_factors.json`, regenerated
  2026-10-28 under the rules of ticket #58, each entry with its exact,
  pretty and numeric form and the chains and irreps that need it; the
  entry count and the uses of every entry are in the generated
  `prototypes/chain_prototype_summary.md`); this chapter does not print
  them. The
  file also holds the entries of the `O` chains through `D3 ⊃ C2` and
  `C3v ⊃ C1h` (`(24, 1/4 ± √2/8)`, `(24, 1/2 ± √2/4)`) and the `T1`, `T2`
  entries of `I ⊃ T ⊃ C3` (`(60, 1/3 ± √5/15)`). Three entries the rules
  of tickets #42 and #43 needed are gone: `(12, 2/3)` for the `T ⊃ C3`
  that the whole-sequence order rule pulled back through a non-principal
  `C3` of `T` (every pulled-back sub-chain is now the subgroup's own chain
  class), and `(24, 5/12)` on `T1` of `O ⊃ D3 ⊃ C2` and an `H` entry of
  `I ⊃ T ⊃ C3` under the first-index designated component. A `--chain`
  document never depends on how an entry was found, only on the entry.

**Provenance of the entries (not normative).** The printed entries were
found by an exact search over candidates `ν = u + i·v` with `u`, `v` short
real sums `a + Σ_j (ζ^j + ζ^{−j})` (`a ∈ {0, 1, 2}`, at most two cosine
terms), filtered by quadratic characters at prime ideals of degree one
and verified by `real_sqrt_F(N / (ν · conj(ν))) = c`, `μ = c · ν`; the
printed form is the hit with both real sums positive. The entries of the
chains of `I` through `T` need products of such candidates (the search of
the research script `docs/research/scripts/cyclotomic_closure.py`, branch
`research/cyclotomic-closure`, exhibits them). No enumeration is part of
the specification: the table is. A search that finds another `μ'` with
`μ' · conj(μ') = N` has found another valid phase, not an error; the
shipped entry stays.

### The conductor of a `--chain` document

A `--chain` document declares one field, as every document does
([Values](./values.md#the-field-of-a-document)), in the version-1 `field`
record without new keys. Its conductor `M_chain` is fixed by the group and
the chain class alone, never by `--irreps`, `--double-group` or the table
subcommand, under the version-1 principle that the declared field is the
smallest field holding every value of the unfiltered document:

- the **chain field** `Q(ζ_L)`, `L = lcm(M_G, M_{H_1}, …, M_{H_m})`, is
  necessary: the chain-adapted matrices contain the characters of every
  step (cyclotomic-closure note, §2.5);
- `M_chain` is the **smallest multiple `M'` of `L` dividing `lcm(L, 24)`**
  such that `M'` is a multiple of the chain conductor of every pulled-back
  sub-chain (`H_1,lib ⊃ φ_1⁻¹(H_2) ⊃ …`) and every squared norm `N` of
  every chain label of every irrep of `G`, ordinary and projective, is a
  Hermitian norm from `Q(ζ_{M'})`: has a canonical norm factor there.

The factor `24` (`√2` and `√3`) is the smallest enlargement any surveyed
chain needed; `T ⊃ D2 ⊃ C2` proves that `L` alone does not suffice:
`(3 ± √3)/6` is no Hermitian norm from `Q(ζ_12)`. The rule evaluates to:

- **Axial families, every chain and every `n`: `M_chain = M_G`.** Every
  step conductor divides `M_G` (`C_d ⊂ C_N` or `D_N` with `d | N`, `D_d`
  with `d | N`, a secondary `C2` or reflection with conductor `4`;
  centro and twisted groups inherit the rotation chain), and
  `N ∈ {1, 1/2}` with `(1 + i)/2 ∈ Q(ζ_{M_G})`. The default-chain document
  of an axial group therefore has the version-1 field, and since its
  matrices equal the canonical table entry by entry
  ([Properties](#properties)), the same exact values byte for byte.
- **Polyhedral groups: one row per chain class.** A type II or III group
  takes the row of its rotation chain (the same squared norms and the
  same `M_G`):

  | Group | Chain class | `L` | `M_chain` | Why |
  | ----- | ----------- | --: | --------: | --- |
  | `T` | `C3` | 12 | 12 | every `N = 1` |
  | `T` | `D2>C2` | 12 | 24 | `(3 ± √3)/6`, no factor in `Q(ζ_12)` |
  | `O` | `C4`, `D4>C4` | 8 | 8 | every `N = 1` |
  | `O` | `T>C3`, `D3>C3`, `T>D2>C2` | 24 | 24 | `L` itself: the shipped entries, or squares |
  | `I` | `D5>C5`, `D5>C2` | 20 | 20 | `N ∈ {1, 1/2, 1/4}`, `1/2 ↦ (1 + i)/2` |
  | `I` | `T>C3` | 60 | 60 | every `N` a Hermitian norm from `Q(ζ_60)` (exhibited exactly; nine shipped entries) |
  | `I` | `T>D2>C2` | 60 | 120 | the pulled-back sub-chain `T ⊃ D2 ⊃ C2` has conductor `24`, and `24 ∤ 60`: the smallest multiple of `L` holding the sub-chain realisation is `120` (five shipped entries in `Q(ζ_120)`; ticket #50 corrected this row from `60`) |
  | `Td` | `S4` ≙ `O>C4`; `C3v>C3` ≙ `O>D3>C3`; `T>C3` ≙ `O>T>C3` | | 8; 24; 24 | |
  | `Th` | `S6`, `C3`, `T>C3` ≙ `T>C3` | | 12 | |
  | `Oh` | `C4h`, `C4`, `O>C4`, `D4h>C4h` ≙ `O>C4`, `O>D4>C4` | | 8 | |
  | `Ih` | `D5d>S10`, `I>D5>C5` ≙ `I>D5>C5` | | 20 | |

  The rows of the chain classes of `list chains` not in this table (for
  example `O ⊃ D4 ⊃ D2 ⊃ C2`) are computed with the same rule by the
  prototype ticket #50 and frozen with the norm-factor data; a chain with
  redundant steps below the cyclic shape takes the row of its class. The
  prototype's evaluation of the rule for every accepted chain class is
  `prototypes/chain_prototype_conductors.json`, generated with
  `uv run python3 prototypes/chain_cg_prototype.py chains`; the number of
  classes per group and the distribution of their conductors (`Td` and
  `O` 8 for the classes through `S4` / `C4` / `D2d` / `D4` and 24 for
  those through `C3v` / `D3` / `T`, `Th` 24 for the classes through `C2v`,
  `D2`, `D2h` or `T>D2>C2`, every class of `D6h` 12) are in the generated
  `prototypes/chain_prototype_summary.md`; this chapter states no counts,
  so that a regeneration cannot leave it stale. Every row the table lists
  is reproduced there.

Consequences:

- `list chains` prints the row as `conductor` ([`list chains`](#list-chains)),
  so a consumer knows the field before reading a document.
- A `--chain` document of a polyhedral group may declare a conductor
  different from the version-1 one (`O --chain "T>C3"`: `24` instead of
  `8`); a consumer reads `field.conductor` per document, as version 1
  already requires. The default chains of `T`, `O`, `Td`, `Th`, `Oh` keep
  the version-1 conductor; those of `I` and `Ih` keep `20`.
- A `cg` document under `--chain` declares `lcm(M_chain, 24)` for the
  polyhedral groups (`24`, or `120` for `I` and `Ih`) and `lcm(M_G, 8)` for
  the axial families, by the rule of
  [Clebsch–Gordan coefficients](./clebsch-gordan.md#the-field-of-a-cg-document);
  its coefficients are the unitary transform of the canonical-basis ones,
  so no squared norm beyond those of the canonical basis arises.
- The minimality check `VL-3` of [Verification](./verification.md) holds
  for every unfiltered `--chain` document by construction (the `large`
  tier runs it for every chain class).

### Errors

A squared norm without a canonical norm factor in the declared field, a
label space of rank other than `0` or `1`, or ranks not summing to the
dimension is an **internal error**, exit code `1`
([Command line](./cli.md#exit-codes-and-errors)): nothing is written to
`-o`, no approximation and no value outside the field is ever emitted. The
message names the quantity and where it arose:

```
pointgroup-tables: error: internal error: no canonical norm factor for '(3 - sqrt3)/6' in Q(zeta_12) (group 'T', chain 'D2>C2', irrep 'E1/2', chain label 'B1 > 1E1/2')
```

The condition is unreachable on shipped data: proven for the axial
families at every `n`, and checked for every chain class and every
Clebsch–Gordan pair of the seven polyhedral groups by the verification
below. A chain accepted by the [grammar section](#chains) can therefore
not fail here except through a defect of the shipped tables.

### Cost

For one irrep of dimension `d` and a chain of `m` steps the work is `m`
projector sums of `|H_i,lib†|` matrices `d × d` over `Q(ζ_L)`, at most `d`
rank computations of `d × d` matrices, one `real_sqrt_F` or table lookup
per block, and `|G†| · d²` matrix elements for the output. `d ≤ 6`
(`I5/2` of `I†`) and the problems never leave the irrep: nothing of the
size of the group algebra is formed. `n` stays unbounded: for the axial
families `L = M_chain = M_G`, the version-1 conductor, so the cost grows
with `n` exactly as a version-1 `matrices` document does; for the
polyhedral groups `M_chain ≤ 60` (`120` for `cg` of `I` and `Ih`, the
field version 1 already uses for the `altmann-herzig` convention of `I`).
The sign search of `real_sqrt_F` is at most `2^7` patterns through the
subfield tower, and the norm-factor table is a lookup.

## Metadata

A document produced with `--chain` carries the top-level metadata key
`chain`, placed after `basis_convention`:

```json
"chain": { "input": "D4>C2", "canonical": "D4>C2[C2]", "symbols": ["O", "D4", "C2"] }
```

| Key | Content |
| --- | ------- |
| `input` | the option value as typed (`"default"` when the default chain was selected) |
| `canonical` | the canonical spelling of the resolved chain (`""` for the empty chain) |
| `symbols` | the normalised symbols of the chain, `G` first, one per group |

The resolved subgroups themselves (their element lists as canonical
indices of `G`, the embedding isomorphisms, the coset decompositions, their
irreps and the branching rules) belong to the
[`subgroup_chain` table](#the-subgroup_chain-table-and-the-chain-labels),
not to the metadata. The metadata key is named `chain`; because metadata
keys and tables share the top-level namespace of a document
([Output schema](./output-schema.md#metadata)), the table is named
`subgroup_chain`.

## The `subgroup_chain` table and the chain labels

Decided in
[The `chain` table and the chain labels of the `matrices` components](https://github.com/hydeik/pointgroup-tables/issues/45).
This section fixes what a `--chain` document adds beyond the `chain`
metadata: the **`subgroup_chain` table**, present in every document
produced with `--chain`, and the **chain labels** on the components of the
`matrices` table. The key-by-key layout and the machine schema are in
[Output schema](./output-schema.md#the-subgroup_chain-table); this section
fixes the content.

### What `--chain` adds, subcommand by subcommand

| Document | Added by `--chain` | Changed by `--chain` |
| -------- | ------------------ | -------------------- |
| every table subcommand | the `chain` metadata key, the `subgroup_chain` table (`tables` lists it right after `"elements"`) | nothing else |
| `matrices` (and `all`) | the `components` array of every irrep entry (the chain labels) | the irrep matrices, now in the chain-adapted basis |
| `elements`, `cayley`, `factor`, `characters` | the two items above | nothing: these tables do not depend on a basis |

`all` stays the union of the single-table documents under the same
options: every one of them carries the same `subgroup_chain` table, so
`all` carries it once. The `basis_convention` metadata is unchanged
(`{ name = "canonical" }`, the only value `--chain` combines with); **the
presence of the `chain` key declares that the matrices are in the
chain-adapted basis of that chain**, the version-1 realisation being the
`m = 0` case. There is no separate basis marker.

### Invariance

The `subgroup_chain` table is a function of the group and the resolved
chain alone. It does not depend on `--double-group`, `--irreps` or
`--values`, so the same table appears, byte for byte, in every document of
one group under one chain whatever the other options; `--label-convention`
only adds `reference_symbol` to the irrep records of `G` it carries (they
stay identical to those of `characters`). The table holds indices and
records only, never a value of the document field, which is why it needs no
`values` form.

### The `subgroup_chain` table

```json
"subgroup_chain": {
  "irreps": [ <irrep record of G>, ... ],
  "items":  [ <subgroup record of H_1>, <subgroup record of H_2>, ... ]
}
```

- `irreps` is the **full** irrep list of `G` (ordinary then projective, in
  canonical irrep order, never filtered by `--irreps`), the records the
  first subgroup's `branching` is indexed by, so that the table reads on
  its own in a document without `characters`.
- `items` holds one **subgroup record** per step, in chain order (`H_1`
  first); the empty chain gives an empty array. A subgroup record carries no
  `index`: its position is the step number, and nothing refers to a step by
  index.

A subgroup record, for the step `K ⊃ H` (`K = H_{i−1}`, `G` for the first
step; `H = H_i` with its abstract group `H_lib` and embedding `φ_i`), keys
in this order:

| Key | Content |
| --- | ------- |
| `group` | the group record of `H_lib` as [Output schema](./output-schema.md#metadata) defines it for the document's group, without `input_symbol`: `{ symbol, family, n, order, type, rotation_group }` |
| `embedding` | the isomorphism `φ_i` as an array of `\|H\|` **`G†` indices**: `embedding[k]` is the index of `φ_i(h̃_k)`, where `h_k` is the `k`-th element of the canonical element list of `H_lib` and `h̃_k` its lifted preimage in `H_lib†`. The index is `≥ \|G\|` exactly when `φ_i` sends the lift of a half turn to the barred lift ([Embedding](#embedding-which-subgroup-a-symbol-picks)); `φ_i(bar(h̃_k))` is the barred partner of `embedding[k]` and is not listed |
| `cosets` | the **left cosets** `g H` of `H` in `K`, one array of `G` indices per coset, each with its representative (the smallest canonical index of the coset) first and the other members ascending; cosets in ascending order of their representatives, so the first coset is `H` itself |
| `double_cosets` | the double cosets `H g H` of `H \ K / H`, one record `{ representative, size, cosets }` each: `representative` the smallest canonical index of the double coset, `size` its number of elements, `cosets` the ascending indices into `cosets` of the left cosets it is the union of; records in ascending order of their representatives, so the first is `H` itself |
| `irreps` | the full irrep list of `H_lib`, ordinary then projective in its canonical irrep order, as irrep records with canonical symbols only (no `reference_symbol`, whatever `--label-convention`: a label convention renames the irreps of `G` and nothing else) |
| `branching` | the **branching rule** of the step: `branching[j]` is the ascending list of indices into this record's `irreps` of the irreps of `H_lib` that occur in the restriction of the irrep `j` of `K` (indexed by `K`'s full irrep list: `subgroup_chain.irreps` for the first step, the previous record's `irreps` afterwards), pulled back through `φ_i`. An ordinary irrep restricts to ordinary irreps, a projective one (with `D(bar h̃) = −D(h̃)`) to projective irreps of `H_lib` with its own factor system. No index repeats: the chain is canonical |

The element and coset indices live in the index space of `G`
(`0 ≤ i < |G|`); only `embedding` uses `G†` indices, under the
[barred-index rule](./output-schema.md#index-spaces) that lets a `G`
document name a barred element as `i + |G|`. Both conventions are the
same in a `--double-group` document, which is what makes the table
independent of that flag. Right cosets, the double-coset members and the
coset representatives' products are not listed: a consumer derives them
from the Cayley table of the `elements` registry, which every document
embeds, and the double cosets are the generator set of the double-induced
technique of the eigenfunction method
([Survey](https://github.com/hydeik/pointgroup-tables/issues/39)), the
reason they are recorded at all.

Example, `D3 --chain C3` (the default chain of `D3`; canonical element
order `E, C3^1, C3^-1, C2'(1), C2'(2), C2'(3)`; irreps of `D3`
`A1 A2 E`, `E1/2 1E3/2 2E3/2`; irreps of `C3` `A 1E 2E`,
`B3/2 1E1/2 2E1/2`):

```json
"subgroup_chain": {
  "irreps": [ { "index": 0, "symbol": "A1", ... }, ..., { "index": 5, "symbol": "2E3/2", ... } ],
  "items": [
    { "group": { "symbol": "C3", "family": "Cn", "n": 3, "order": 3, "type": "proper", "rotation_group": "C3" },
      "embedding": [0, 1, 2],
      "cosets": [[0, 1, 2], [3, 4, 5]],
      "double_cosets": [ { "representative": 0, "size": 3, "cosets": [0] },
                         { "representative": 3, "size": 3, "cosets": [1] } ],
      "irreps": [ { "index": 0, "symbol": "A", ... }, ..., { "index": 5, "symbol": "2E1/2", ... } ],
      "branching": [[0], [0], [1, 2], [4, 5], [3], [3]] }
  ]
}
```

`E ↓ C3 = 1E + 2E` gives `branching[2] = [1, 2]`; both `1E3/2` and `2E3/2`
of `D3†` restrict to `B3/2` of `C3†`, so `branching[4] = branching[5] = [3]`.
Under `O --chain "D4>C4"` the second record's `branching` is indexed by the
irreps of `D4` (`A1 A2 B1 B2 E`, `E1/2 E3/2`) and its `cosets` are the two
cosets of `C4` in `D4`, as `G` indices of `O`.

### Chain labels of the `matrices` components

Under `--chain`, every irrep entry of the `matrices` table carries, after
`frobenius_schur` and before the value forms, a `components` array with
one **component record** per component, in the
[component order](#component-order) of the chain-adapted basis:

```json
{ "index": 2, "symbol": "E", "dimension": 2, "projective": false, "letter": "E",
  "frobenius_schur": 1,
  "components": [
    { "index": 0, "labels": [1], "symbols": ["1E"] },
    { "index": 1, "labels": [2], "symbols": ["2E"] }
  ],
  "exact": [ ... ], "numeric": [ ... ] }
```

| Key | Content |
| --- | ------- |
| `index` | the component's position, `0 ≤ index < dimension`: the row and column of the matrices it labels, and the index the Clebsch–Gordan tables refer to ([Clebsch–Gordan coefficients](./clebsch-gordan.md)) |
| `labels` | the chain label `(σ_1, …, σ_m)` as indices: `labels[j]` is the index, in `subgroup_chain.items[j].irreps`, of the irrep of `H_{j+1},lib` that labels the component |
| `symbols` | the same chain label as canonical Mulliken symbols, `symbols[j]` the `symbol` of that irrep record; display only, like every label of a document |

`labels` and `symbols` have one entry per step, so under the empty chain
every component record is `{ "index": k, "labels": [], "symbols": [] }`.
The entries of `labels` are consistent with `branching`: `labels[0]` occurs
in `subgroup_chain.items[0].branching[<irrep index>]`, and `labels[j]` in
`items[j].branching[labels[j−1]]`. The component records are the only
addition to an irrep entry; the irrep record proper (`index` …
`frobenius_schur`) stays identical to that of `characters`, and a document
without `--chain` has no `components` key.

### Formats

The table and the component records contain integers, strings and arrays
only. In YAML and TOML they are dumps of the same object
([Output schema](./output-schema.md#json-formatting)): in TOML the table is
`[subgroup_chain]` with `irreps` and `items` as arrays of inline tables,
and an empty `items` or `labels` is the empty array `[]`, so no `null` is
ever needed.

### Verification of the table

An implementation must check, for every `--chain` document (these are the
rows `OS-7`, `OS-8` and `GOLD-2` of the [check ledger](./verification.md#check-ledger)):

- `subgroup_chain` is present exactly when `chain` is, listed in `tables`
  right after `"elements"`, and is byte-identical across every subcommand,
  `--double-group`, `--irreps` and `--values` for one group and chain;
  `subgroup_chain.irreps` equals the unfiltered irrep list of `G`;
- every subgroup record is consistent with the Cayley table of `G†`:
  `embedding` is injective, maps the canonical element list of `H_lib` to
  elements of the previous group, and is a homomorphism of double groups
  (the product of two images, read in the Cayley table of `G†`, is the
  image of the product in `H_lib†`); `cosets` partitions the previous
  group's elements into `|K| / |H|` left cosets of the recorded `H`, with
  representatives and members ordered as stated; `double_cosets` partitions
  `cosets`, each `size` equals `|H|` times the number of its cosets, and
  every double coset is closed under left and right multiplication by `H`;
- `branching` equals the subduction computed from the characters of `K`
  and `H_lib` through `embedding` (`⟨χ_Γ ∘ φ, χ_σ⟩` over `H_lib†`), has no
  repeated index, and the dimensions add up;
- the `components` of every `matrices` entry are in the
  [component order](#component-order) (lexicographic in the steps,
  exponent order at a cyclic step), consistent with `branching` as stated, and agree with the
  [block structure](#properties): restricted to `φ_j(H_{j},lib†)` the
  matrices are block-diagonal in the components grouped by
  `labels[0 … j−1]`, and on the terminal diagonal with the characters of
  `labels[m−1]`;
- golden documents: every table subcommand × `json`, `yaml`, `toml` for
  `C3v`, `D6h`, `Oh`, `Ih` under `--chain default`, plus `matrices O
  --chain "D4>C4"` and `matrices D6 --chain "D3[C2'']>C3"`, equal the
  checked-in files byte for byte.

## `list chains`

The utility subcommand `list chains <group>` lists the chains the group
accepts, in `--format` (`-f`, `-o` only, like the other `list`
subcommands):

```json
{
  "group": { "symbol": "O", "input_symbol": "o", ... },
  "default": "C4",
  "chains": [
    { "chain": "C4", "symbols": ["O", "C4"], "default": true, "conductor": 8 },
    { "chain": "D3>C3", "symbols": ["O", "D3", "C3"], "default": false, "conductor": 24 },
    ...
  ]
}
```

`group` is the group record of [Output schema](./output-schema.md#metadata);
`chain` is the canonical spelling; the array holds **every canonical
spelling accepted by `--chain`** (canonical for both kinds, cyclic-shape
terminal, one entry per canonical spelling, redundant steps below the cyclic
shape excluded), each entry being the chain the embedding rule resolves the
spelling to (the reached twin where the spelling has
[twin chain classes](#embedding-which-subgroup-a-symbol-picks); an unreached
twin has no entry, as no `--chain` value names it), sorted by number of
steps and then by canonical spelling in ASCII order; `conductor` is the
chain conductor of that chain
([The conductor of a `--chain` document](#the-conductor-of-a-chain-document)).
The group argument follows the grammar and exit codes of the
table subcommands (`2`, `3`, `5`). The document is not covered by
`schema.json`.

## Errors

| Code | Condition | Example message |
| ---- | --------- | --------------- |
| `2` | the value does not match the grammar: empty, whitespace, a bad separator, a malformed anchor, `default` with steps or an anchor; `--chain` with a non-`canonical` `--basis-convention`; `--chain` on a utility subcommand | `pointgroup-tables: error: invalid chain 'D4>': expected a group symbol after '>'` |
| `3` | a step's symbol names no group (the group-argument rule) | `pointgroup-tables: error: unknown group 'Xh' in chain 'Xh>C4'` |
| `5` | a step's symbol has an unsupported `n` | `pointgroup-tables: error: unknown group 'S3' in chain 'S3': S2n takes even subscripts only` |
| `6` | **invalid chain**: a step is not a proper subgroup of the previous group (after anchors); the chain is not canonical for both kinds; the terminal is not of the cyclic shape; an anchor names no class of `G` or no subgroup contains it | `pointgroup-tables: error: invalid chain 'D4>C3' for 'O': 'C3' is not a subgroup of 'D4'` |

Hints follow the [error format](./cli.md#exit-codes-and-errors): a code-`6`
error adds `hint: run 'pointgroup-tables list chains O' for the accepted
chains`, and a chain rejected for stopping above the cyclic shape names a
completion (`hint: 'T>D2' stops above the cyclic shape; try 'T>D2>C2'`).

## Verification

An implementation must pass (rows `GC-1` to `GC-5`, `CB-1` to `CB-5`,
`VL-4` to `VL-6`, `OS-7`, `OS-8`, `GOLD-2`, `REF-9`, `REF-10`, `REF-12`
and `MUT-8` to `MUT-10` of the
[check ledger](./verification.md#check-ledger)):

- the grammar table: at least `"D4>C4"` on `O` → `0`; `"d4>c4"` → `0` with
  `symbols` `["O", "D4", "C4"]`; `"default"` → `0`; `"DEFAULT"` → `0`;
  `"O>D4>C4"` → `6`; `"D4>"` → `2`; `"D4 > C4"` → `2`; `""` → `2`;
  `"default>C4"` → `2`; `"Xh"` → `3`; `"S3"` on `D6h` → `5`; `"D2"` on
  `T` → `6` (stops above the cyclic shape); `"C5"` on `I` → `6` (ordinary
  only); `"C2"` on `T` → `6` (projective only); `"C3[sigma_h]"` on `D3h`
  → `6` (no `C3` contains a `sigma_h`); `"C4"` with
  `--basis-convention altmann-herzig` → `2`;
- the default-chain table above, group by group for the polyhedral groups
  and for the axial families at `n = 1 … 12` (the empty chain for the
  cyclic shape);
- the embedding examples of the [consequences](#embedding-which-subgroup-a-symbol-picks):
  the resolved subgroup of `"C2"` on `D6` is generated by `C2'(1)`, of
  `"C1h"` on `C6v` by `sigma_v(1)`, of `"D4>C2"` on `O` by `C2(2)`, of
  `"D3[C2'']>C3"` on `D6` by `C3^1` and `C2''(1)`; the resolved `D4` of `O`
  contains `C4^1(1)` and the resolved `C3v` of `Td` contains `C3^1(1)`;
- `list chains` agrees with the catalogue of ticket #41 (canonicity for
  both kinds, terminal rule) for the seven polyhedral groups and the axial
  families at `n ≤ 12`, one entry per canonical spelling, the catalogue's
  chain classes collapsing exactly by the table of twin chain classes (none
  in the axial families); every canonical spelling resolves to the class it
  was derived from or to its reached twin, and the spellings of the twin
  table resolve to the reached twin named there; `list chains O` holds
  `"C4"` first with `default = true`;
- the metadata of every `--chain` document holds `input`, `canonical` and
  `symbols` as above, and `--chain` with the canonical spelling of a
  document's chain reproduces that document's `symbols`, `subgroup_chain`
  table and matrices byte for byte;
- the chain-adapted basis: for every irrep of `T†`, `O†`, `I†`, `D6h†` and
  every chain of `list chains`, the label spaces have rank `≤ 1` and ranks
  summing to the dimension, the chain-adapted matrices satisfy the
  multiplication table with the factor system, are unitary, and are
  block-diagonal as the [properties](#properties) state, with the step-1
  blocks equal to the sub-chain realisation (`matrices O --chain "D4>C4"`
  contains the `E` blocks of `matrices D4 --chain C4`); `matrices D6
  --chain C6` and `matrices D6h --chain C6h` equal the canonical tables
  entry by entry; `matrices O --chain C4` equals the canonical table
  entry by entry with the chain labels added (`CB-4`); `matrices Oh --chain "D4h>C4h"` equals
  `matrices O --chain "D4>C4"` with parity signs; every designated
  component equals the conjugate of the canonical norm factor of its
  squared modulus, and is real positive wherever the real square root lies
  in the document field;
- the exact computation ([Exact computation](#exact-computation)):
  `real_sqrt_F` returns the root for `2` in `Q(ζ_8)`, `3` in `Q(ζ_12)`,
  `(2 + sqrt3)/6` in `Q(ζ_12)`, `5` in `Q(ζ_5)` and `7/20 · 20/7 = 1` in
  every field, and proves non-squareness for `2` in `Q(ζ_12)`, `7` in
  `Q(ζ_120)` and `(3 ± sqrt3)/6` in `Q(ζ_24)`, and agrees with an
  independent high-precision evaluation on random squares of `Q(ζ_24)`,
  `Q(ζ_60)` and `Q(ζ_120)` with coefficient denominators up to `10⁶`;
  every shipped norm-factor entry satisfies `μ · conj(μ) = N` exactly and
  every printed entry of the table above is reproduced; every squared
  norm of every chain label of every irrep of every chain class of
  `list chains` of the seven polyhedral groups, and every squared norm of
  every Clebsch–Gordan pair of those groups, is a square in the declared
  field, equals `1/2` with `(1 + i)/2`, or has a shipped entry (`large`
  tier: the internal error of the [Errors](#errors) subsection is
  unreachable); the chain conductor of every chain class equals the table
  above where the table has a row, and is minimal (`VL-3`) for the
  unfiltered `matrices` document of every chain class, with `M_chain = M_G`
  for the axial families at `n ≤ 12`; `list chains` prints that value as
  `conductor`.
