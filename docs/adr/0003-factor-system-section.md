---
status: accepted
---

# The factor system is fixed by the positive-preimage section on integer models of the binary rotation groups

The double group `G†` is built from the binary cover `R†` of the rotation
group (`Z_2n`, `Dic_n`, or the binary polyhedral groups `2T`, `2O`, `2I` as
unit quaternions with coordinates in `(1/2) Z[α]`, `α ∈ {1, √2, φ}`). The
section that fixes the spin factor system `ω` picks, for every rotation, the
preimage that is *positive*: residue in `(−n/2, n/2]` for the axial groups,
first non-zero quaternion coordinate positive for the polyhedral groups. The
inversion lifts trivially. In SU(2) terms this is the rotation angle taken in
`(−π, π]` with binary-rotation axes in a positive hemisphere, the convention
of Altmann & Herzig and Bilbao. See `docs/spec/double-groups.md`.

## Considered Options

- **Angle in `[0, 2π)`** (the convention of many physics texts): the lift of
  `C_n^k` is `a^k` for `0 ≤ k < n`, giving a "carry" cocycle
  `ω(k, k') = −1 iff k + k' ≥ n`. Rejected because it disagrees with the
  primary reference (Altmann & Herzig) on every rotation by more than `π`,
  and because the canonical labels already reduce exponents to `|j| ≤ m/2`,
  so a `(−π, π]` section is the one that makes `s(C_m^j) = s(C_m)^j`.
- **Angle in `(−π, π]`, positive hemisphere** (chosen): agrees with the
  primary reference, and can be stated as a single sign rule on the integer
  models.
- **Permutation models for `2T`, `2O`, `2I`**: would keep one kind of model
  for every rotation group, but permutations carry no sign, so the section
  could only be "the preimage with the smaller closure index", which no
  reference convention would match. Quaternions were chosen for the covers
  while the permutation models of `T`, `O`, `I` are kept, joined by a
  projection fixed through generator images.
- **Section as part of the basis convention**: a reference's choice of which
  preimage is `g` and which `ḡ` could be absorbed as a sign in the irrep
  basis. Rejected in favour of a `barred` override table in the *label*
  convention, since which element carries the tilde is a labelling matter and
  the factor table must follow it, not only the matrices.

## Consequences

- `ω(C2, C2) = −1` for every binary rotation, `ω(C3, C3) = −1`,
  `ω(C_m^j, C_m^{−j}) = +1`. Consumers used to the `[0, 2π)` convention will
  see a factor system that differs by a coboundary; the spec says so.
- `2T`, `2O`, `2I` are a second kind of integer model (quaternions over a
  quadratic ring) next to the permutation models; the generator pairs and
  the projection are part of the specification, since changing them changes
  which canonical element is barred.
- Reference conventions that bar different elements are supported by a
  per-group `barred` table in the label convention file, which transforms
  the factor table and the projective irreps by the corresponding
  coboundary. Version 1 ships none.
