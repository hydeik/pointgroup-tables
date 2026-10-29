---
status: accepted
---

# Clebsch–Gordan coefficients are fixed by a trial-state projection, a symmetric-first multiplicity rule and the closure field of the multiplicity-free norms

The Clebsch–Gordan coefficients of a point group are defined only up to a
phase per coupled copy and a unitary mixing of the copies of an irrep that
occurs more than once in a product, and the literature fixes both by hand,
table by table (Koster, Butler, the eigenfunction-method papers of Chen,
Fan and Ping). The tool must fix them by a rule that is stated once, works
for every group, every `n`, every basis a document can name, and keeps
every value in one cyclotomic field per document (ADR-0002). We chose: the
coupled copy is the projection of the first product basis vector the
isotypic projector does not annihilate; for `α = β` the symmetric copies
come first, then the antisymmetric ones, the trial-state order deciding
inside a part; the first non-zero coefficient of the designated component
is the conjugate of the canonical norm factor of its squared modulus, which
is the real positive square root whenever that root lies in the document's
field; and a `cg` document enlarges the basis conductor to
`lcm(M, 8)` (axial) or `lcm(M, 24)` (polyhedral), the closure of the
multiplicity-free norms, and no further. See `docs/spec/clebsch-gordan.md`.

## Considered Options

- **Closure of the multiplicity-free norms with the norm-factor fallback**
  (chosen). Every coefficient of every axial, `T`-type and `O`-type group
  and all but four copies of `I` (fourteen of `Ih`) have the familiar
  real positive designated coefficient and read like the printed tables;
  the field is `Q(ζ_24)` / `Q(ζ_120)`, the same as the `altmann-herzig`
  polyhedral data and the `--chain` documents, so a `cg --chain` document
  needs no further enlargement. The cost: the four `I` copies whose norms
  carry the prime `7` (`T1 ⊗ I5/2 ⊃ I5/2`, `I5/2 ⊗ I5/2 ⊃ T1`) have a
  unit-modulus complex designated coefficient, and one rule has two
  faces.
- **Keep the version-1 conductor and always use the norm factor.** The
  smallest documents (`Q(ζ_8)` for `O`, `Q(ζ_20)` for `I`) and the one rule
  that the chain-adapted bases need anyway. Rejected: `1/√3` of `O` would
  print as `(1 − i√2)/3` and `1/√2` of every axial group as a unit multiple
  of `(1 − i)/2`, unlike every reference table, which makes the reference
  comparison a comparison up to phase everywhere and the tables unfamiliar
  to every reader.
- **Real positive everywhere, conductor `840` for `I`-type groups.** Every
  designated coefficient real positive. Rejected: degree 192 instead of
  32, six times the value size of every `I` and `Ih` document, for four
  copies; and the rule is impossible in the chain-adapted bases, where the
  norms have no cyclotomic square root (ticket #40), so it would still need
  the fallback.
- **The papers' multiplicity convention** (copies labelled by hand-chosen
  trial states, not orthogonal). Rejected: it needs a table of trial states
  per group, and the copies are not orthonormal, so the coupling matrix is
  not unitary.
- **Trial-state order without the symmetric/antisymmetric split.** One rule
  for every pair, the rule of the research survey. Rejected: it gives no
  exchange relation for `α = β` (the copies mix symmetric and antisymmetric
  parts), while the split costs nothing (the swap commutes with the
  product representation, also for projective irreps) and yields the
  textbook `[E²] = A1 + E`, `{E²} = A2`.
- **A second class operator (CSCO-III) to separate the copies of `α ≠ β`.**
  Chen's intrinsic operator. Rejected: it ties the multiplicity label to a
  chain and does not resolve every case either (`G3/2 ⊗ G3/2 ⊃ T1` of `O`
  stays degenerate under the symmetric split and needs an order anyway).
