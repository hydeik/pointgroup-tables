---
status: accepted
---

# Chain-adapted phases follow the tool's own rule (norm factor, subduction-exact blocks, designated component), not the papers'

The eigenfunction method fixes a chain-adapted basis only up to one phase
per basis vector, and every Chen–Fan–Ping paper closes the gap with a
different external convention (agreement with Butler's tables, real
symmetry-adapted functions, hermiticity of the projection operators,
associated-irrep relations). A "designated component real and positive"
rule, the natural single rule, is impossible in exact arithmetic: for
`T ⊃ D2 ⊃ C2`, `O ⊃ T ⊃ C3` and `O ⊃ D3 ⊃ C3` the squared norms of the
chain-adapted vectors have square roots in no cyclotomic field (ticket #40).
We fix every phase by three rules of our own: the **norm-factor rule** (a
designated quantity equals the conjugate of the canonical norm factor of
its squared modulus, which is the real positive square root whenever that
root lies in the document field), **subduction-exact blocks** (the
restriction of an irrep to the first subgroup is, block by block, the
subgroup's own chain-adapted matrices for the pulled-back chain), and the
**designated canonical component** (the first non-vanishing canonical
component of the first vector of every block). The papers' matrices become
reference data compared up to a permutation and a diagonal unitary. See
`docs/spec/group-chains.md`, "Chain-adapted basis".

**Amended by ticket #57 (2026-10-28).** The labels of a cyclic-shape step
are ordered in **exponent order** (by the exponent of the label's
eigenvalue at the step's own principal element, descending, as version 1
orders canonical components), not in the subgroup's canonical irrep order
as ticket #43 first decided. The exact survey of ticket #49 showed that the
irrep order permuted the components of every polyhedral default chain and
swapped two components of some `E'`/`E''` irreps of every odd `Dnh`
(the `1E`/`2E` labels of an odd `Cnh` are assigned at `C_n`, the canonical
order at `S_n^-1`); with the exponent order every default chain but those
of `I` and `Ih` reproduces the version-1 `matrices` entry by entry. The
cost, accepted: a cyclic step's order follows its exponents, not the
printed irrep list (`B, 1E, A, 2E` for `C4`, `B, A` for a secondary `C2`),
and it applies uniformly to every cyclic step, pulled-back sub-chains
included. See `docs/spec/group-chains.md`, "Component order".

**Amended by ticket #58 (2026-10-28).** The designated canonical component
of a block is the first non-vanishing one **among the components fixed by
the principal element of `G`** (eigenvalue `1`), and the first
non-vanishing one of all only when the block has none (every projective
irrep). Ticket #43's "first non-vanishing" made `D'(C3)` on `T` of
`T ⊃ D2 ⊃ C2` the cyclic permutation times `ζ_3` on every entry, where
[CF98b] eq 17 prints the real permutation: the three `D2` blocks were all
designated at the `m = −1` eigenvector of `C3`. Designating at the
`C3`-fixed component makes the principal element permute blocks without
phases wherever it normalises the first subgroup; the shipped norm-factor
data shrink (`(24, 5/12)` leaves, `I ⊃ T ⊃ C3` trades an `H` entry for
`T1`/`T2` entries `1/3 ± √5/15`), no conductor changes, and the chains
through the principal element and every projective irrep are untouched.
The cost, accepted: a two-clause rule instead of one. See
`docs/spec/group-chains.md`, "The basis" and "Properties".

## Considered Options

- **Norm-factor rule with subduction-exact blocks and the designated
  canonical component** (chosen). One primitive shared with the
  Clebsch–Gordan phase rule, every value in one cyclotomic field per
  document, the subgroup's tables visible verbatim inside the group's, and
  every one-step chain through the principal element reproducing the
  canonical basis, order included (amendment of ticket #57). The cost: polyhedral chain matrices inherit the frozen canonical
  phases of ADR-0004 (the axial ones are formula-defined), and no paper's
  table is matched entry by entry.
- **Designated component real and positive.** Rejected: provably leaves
  every cyclotomic field on some canonical chains; a quadratic extension of
  the value model for a convention no reference shares is not worth it.
- **The papers' conventions, entry by entry.** Rejected: there is no single
  rule across the papers, several are stated as agreement with another
  author's tables, and none covers the chains the tool accepts beyond the
  printed ones.
- **Intertwining-element rule for the phases between blocks** (the first
  element of the group with a non-zero matrix element between two blocks'
  first vectors, made canonical by the norm-factor rule). It makes the
  matrices independent of the frozen canonical phases, but it re-phases the
  canonical basis even for `Dn ⊃ Cn` (the projective `C2'` matrix becomes
  `[[0, 1], [−1, 0]]`), so the default chain would no longer coincide with
  the canonical basis of version 1. Rejected for that reason.
- **Phases free within each block too** (every vector by its designated
  component). Simpler, but the subgroup's blocks would agree with the
  subgroup's own tables only up to phases, losing the property the papers
  call "subduction built into the matrices".
