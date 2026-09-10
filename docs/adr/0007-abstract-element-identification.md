---
status: accepted
---

# A label convention identifies elements abstractly, so that the canonical irrep rules reproduce the reference's symbols

A label convention maps canonical labels to a reference's labels through an
isomorphism between the abstract group and the reference's group. That
isomorphism is fixed by anchors (`E`, the principal rotation, the first
secondary element, `i`), and for the groups with two classes of secondary
elements there are two natural choices: the **geometric** one, which sends
canonical direction `0` to the reference's `x` axis, and the
**letter-preserving** one, which sends the canonical `sigma_v(1)` to the
reference's `sigma_v1`. In `C2v` and `C6v` (every `Cnv` with
`n ≡ 2 (mod 4)`) the two differ, because the reference's `sigma_v` planes
contain the master group's `C2'` axes while the canonical `sigma_v(1)` is
perpendicular to direction `0`.

We choose the letter-preserving isomorphism. See
`docs/spec/label-conventions.md` (the twist) and the anchors in
`tools/ah_book.py`.

## Considered Options

- **Letter-preserving isomorphism** (chosen). Canonical `sigma_v(1)` prints
  as `sigma_v1`, and the canonical rule "`B1` is symmetric under
  `sigma_v(1)`" reproduces the book's `B1` (symmetric under `sigma_y` in
  `C2v`, T 50.4; under `sigma_v` in `C6v`, Table 17.6) without an irrep
  override. The cost is one extra swapped lift per group and a rule (the
  twist) in the spec. The map is fixed by the letters both layers share,
  which is what a reader compares.
- **Geometric isomorphism.** Canonical direction `0` is the book's `x`
  axis everywhere, and the direction table applies without a twist; but
  canonical `sigma_v(1)` prints as `sigma_d1` in `C6v` and `sigma_x` in
  `C2v`, and the book's `B1` is the canonical `B2`, so every `Cnv` with
  `n ≡ 2 (mod 4)` needs an `[irreps.overrides]` swap. Two layers that
  disagree on letters would confuse every document that shows both.
- **Change the canonical rule** so that `sigma_v(1)` of `Cnv` contains
  direction `0` instead of being perpendicular to it. Rejected: the
  canonical labels are convention-free and fixed (ADR-0001, ticket #3);
  a reference's habits must not leak into them.
