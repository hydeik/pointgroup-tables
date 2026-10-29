---
status: accepted
---

# A chain names its subgroups by Schoenflies symbol and class anchor, and a lexicographic rule picks the concrete subgroup

A Schoenflies symbol names an abstract group, but a point group contains
the same abstract group several times: as conjugate subgroups (the three
`D4` of `O`) and as non-conjugate subgroups with the same symbol (the
central `C2` and the `C2'`, `C2''` subgroups of `D4`; the `sigma_h`,
`sigma_v` and `sigma_d` planes of `D6h`). The chain-adapted basis of
`--chain` depends on which one is meant, and the chain labels depend on how
its elements are identified with the subgroup's own (`C_n^1` against
`C_n^-1` decides `1E` against `2E`). We resolve a symbol by a fixed rule:
the admissible isomorphism (parity, label invariants and lifts preserved)
whose image is lexicographically smallest in the canonical indices of `G`,
after dropping the isomorphisms that break canonicity; the only input a
user has is an optional **class anchor** `[C2']`, `[sigma_h]` naming a
class of `G` the subgroup must meet. Conjugate subgroups within one class
are never distinguishable. See `docs/spec/group-chains.md`.

**Amended by ticket #58 (2026-10-28).** The lexicographic comparison runs
over the element list of the symbol's double group **sorted by descending
element order** (ties by canonical index), not over the list in its own
order. The prototype of ticket #50 showed that the original order let the
half turns of `T`'s element list, listed before its `C3`, pick `T` in `I`
through `C2(1)` and orient the unique `T` of `O`, `Td` and `Th` through a
non-principal `C3`: the stated consequence "the conjugate through the
principal element" was false, and the pulled-back sub-chain of
`O ⊃ T ⊃ C3` was not `T ⊃ C3` itself (an extra norm-factor entry, the
subgroup's tables not visible inside the group's). With the elements of
highest order compared first every consequence holds and every pulled-back
sub-chain is the subgroup's own chain class, at no change of any conductor.
The alternative "compare the principal element's image first, then the
rest in canonical order" was not taken: it can choose another conjugate
and orientation of `T` in `I`, and the descending order needs no notion
beyond the element orders.

## Considered Options

- **Symbol plus class anchor, lexicographic rule** (chosen). Every chain
  of the literature is typed as the paper writes it (`"C3"`, `"D4>C4"`,
  `"D5>C5"`, `"C2"` for Chen's real `Dn ⊃ C2x` basis), the anchors reuse
  the class labels every document already prints, and every subgroup
  class that admits a canonical chain is reachable. The cost: the rule is
  the definition, so the picked conjugate and orientation are a consequence
  of the canonical element order rather than a stated geometric choice;
  `list chains` and the `chain` table make the result explicit.
- **Primed group symbols** (`C2''`, `D2'`, `C2v''`), the notation of
  Altmann & Herzig and of the research note. Rejected: the meaning of a
  prime varies by case (in `D6` both `D3` classes contain a secondary
  element, in `O` only `D2'` does, the planes of `C1h` are `h`/`v`/`d`
  rather than primed), so it needs a table per group instead of one rule.
- **Subgroups by generators, axis indices or coordinates** (`D4[C4^1(2)]`,
  a generator list, a rotation matrix). Rejected: it reopens the choice of
  conjugate, which changes nothing but a similarity transform of the
  matrices, and it makes the chain an input that fixes coordinates, which
  the glossary forbids a chain to do.
- **Reject ambiguous strings** and require an anchor whenever a symbol
  admits several classes. Rejected: `"C2"` on `Dn` and `"D4>C2"` on `O`
  are the papers' chains and should work as typed; the canonical spelling
  in the metadata and in `list chains` gives the unambiguous form to
  anyone who needs it.
- **Let the user pick a conjugate.** Rejected for version 2: conjugate
  subgroups give bases related by one irrep matrix, and the `chain` table
  records the representative, so a consumer can transform; the map rules
  subgroups named by coordinates out of scope.
