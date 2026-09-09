# pointgroup-tables

A command-line tool that emits machine-readable tables (group elements, Cayley
table, factor system, character tables, irreducible-representation matrices)
for finite point groups and their spin double groups.

## Language

### Groups

**Point group**:
A finite group of orthogonal transformations of three-dimensional space that
fix the origin, identified by its Schoenflies symbol (e.g. `C3v`, `D6h`, `Ih`).
The axis order `n` of the axial families is any positive integer; the group
need not be crystallographic.
_Avoid_: symmetry group, crystal class, Hermann–Mauguin symbol as an identifier

**Family**:
One of the infinite series (Cn, Cnv, Cnh, S2n, Dn, Dnd, Dnh) or the isolated
polyhedral groups (T, Td, Th, O, Oh, I, Ih) into which point groups are sorted.
A family plus an axis order names exactly one point group.
_Avoid_: type, class (reserved for conjugacy class)

**Double group**:
The group of order 2|G| obtained by lifting a point group G to SU(2), so that
a rotation by 2π is a distinct element from the identity. Written G†.
_Avoid_: spin group, spinor group

**Abstract group**:
The group isomorphic to a point group, defined without any coordinate system
as a subgroup of `rotation group × parity`. Every computation in this project
happens on the abstract group; geometry enters only through labels.
_Avoid_: matrix group, coordinate representation

**Rotation group**:
The finite group of proper rotations (Cn, Dn, T, O or I) from which a point
group is built, handled as an exact integer object (residues, pairs, or
permutations).
_Avoid_: proper subgroup, SO(3) subgroup

**Parity**:
The sign +1 (proper) or −1 (improper) attached to every element of a point
group; part of the abstract group, since groups with the same bare structure
but different parities have different double groups.
_Avoid_: determinant, chirality, handedness

**Group type**:
One of proper (Type I, all elements proper), centro (Type II, contains the
inversion, `R × {±1}`) or twisted (Type III, improper elements are the
rotations outside an index-2 subgroup of `R`).
_Avoid_: kind, class

**Element**:
A member of the abstract group of a point group, referred to by a label (e.g.
`E`, `C3+`, `sigma_v(1)`, `i`) supplied by a label convention. Elements are
handled as abstract symbols, never as floating-point matrices.
_Avoid_: operation, operator, symmetry

**Cayley table**:
The multiplication table of a point group: the product of every ordered pair of
elements.
_Avoid_: multiplication table, product table, group table

### Representations

**Factor system**:
The function ω(g, h) on pairs of elements that measures the failure of a
projective representation to be a true representation. In this project ω is
always the spin factor system, taking values ±1, induced by the double group.
_Avoid_: factor table (as a concept; it is the *table* of the factor system),
cocycle, multiplier, Schur multiplier

**Projective representation**:
A double-valued (spinor) representation of a point group: a representation of
its double group that is not a representation of the point group itself, or
equivalently a representation of the point group with the spin factor system.
_Avoid_: ray representation, extra representation, spinor representation

**Ordinary representation**:
A single-valued representation of a point group, i.e. one with a trivial
factor system.
_Avoid_: vector representation, true representation, linear representation

**Irrep**:
An irreducible representation, ordinary or projective, of a point group. Each
irrep has a Mulliken symbol as its name and may carry a Bethe symbol as an
alias.
_Avoid_: irreducible rep, species

**Character table**:
The table of characters of every irrep on every conjugacy class of a point
group (or of its double group, for projective irreps).

**Irrep matrix table**:
The table of explicit matrices of an irrep on every element, given in a basis
fixed by a basis convention.
_Avoid_: representation matrices (without naming the convention)

**Conjugacy class**:
A set of elements conjugate to one another within the group. Characters are
constant on a conjugacy class.
_Avoid_: class (alone, when ambiguous)

**Frobenius–Schur indicator**:
The value +1, 0 or −1 attached to an irrep stating whether it is real,
complex, or quaternionic (pseudo-real).
_Avoid_: reality, realness type

### Conventions

**Label convention**:
A named, selectable rule that assigns a geometric label to every element of an
abstract group (e.g. which generator is called `C3+`, which reflection is
`sigma_v(1)`). A label convention implies an orientation of the group in
space, but the tool itself never fixes coordinates.
_Avoid_: setting convention, orientation, standard setting, naming scheme

**Basis convention**:
The rule fixing the basis in which irrep matrices are written, selectable on
the command line. Irrep matrices are only defined up to unitary equivalence,
so a table is meaningless without naming its basis convention.
_Avoid_: phase convention, gauge

**Mulliken symbol**:
The name of an irrep in the A/B/E/T (and G/H for the icosahedral groups)
notation, extended to arbitrary axis order (E_k for ordinary, E_{k/2} for
projective irreps).

**Bethe symbol**:
The Γ_i name of an irrep, kept only as an alias of the Mulliken symbol.

### Number representation

**Exact value**:
A number written as an element of a cyclotomic field (a symbolic expression
such as `exp(2*pi*I/3)`), never rounded.
_Avoid_: symbolic value, analytic value

**Numeric value**:
The double-precision complex approximation of an exact value, emitted next to
it for consumers that cannot evaluate exact values.
_Avoid_: float, approximate value
