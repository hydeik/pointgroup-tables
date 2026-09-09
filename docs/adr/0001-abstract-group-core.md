---
status: accepted
---

# Point groups are handled as abstract groups with no fixed coordinate system

Every point group is represented internally by an isomorphic abstract group
(generators and relations, or its Cayley table); the tool never fixes a
Cartesian frame or manipulates rotation matrices. Geometric labels such as
`C3+` or `sigma_v(1)` are assigned by a selectable label convention on top of
the abstract group, so the same tables can be relabelled to match different
reference conventions (Altmann & Herzig, Bradley & Cracknell, Bilbao) without
recomputation.

## Considered Options

- **Matrix-first**: elements as 3x3 orthogonal matrices in a standard
  setting; labels derived from axis and angle. Rejected because it bakes one
  orientation into every table and forces floating-point or heavy symbolic
  linear algebra just to multiply elements.
- **Abstract-first** (chosen): elements are symbols, products come from the
  group presentation, geometry is a labelling layer.

## Consequences

- The abstract group must carry the parity (proper/improper) of every
  element: `Ci` and `C2` are isomorphic as bare groups but have different
  double groups. The model is therefore `G ⊆ R × {±1}` with `R` a rotation
  group (see `docs/spec/abstract-groups.md`).

- Label conventions and basis conventions must be specified in terms of
  generators of the abstract group, not in terms of axes.
- The spin double group is defined as an abstract central extension, and the
  factor system by a chosen section of that extension, rather than through an
  SU(2) lift computed from rotation angles.
