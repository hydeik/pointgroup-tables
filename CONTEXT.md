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
a rotation by 2π is a distinct element from the identity. Written G†. Built
from the abstract group as a central extension; the inversion lifts
trivially.
_Avoid_: spin group, spinor group

**Full turn**:
The central element Ē of a double group covering the identity: the rotation
by 2π. Its canonical label is `~E`.
_Avoid_: bar-E as a concept name, minus identity, −E

**Section**:
The choice, for every element g of a point group, of one of its two
preimages in the double group, written s(g). The factor system is entirely
determined by the section; the canonical section takes the positive
preimage (rotation angle in (−π, π], binary-rotation axis in the positive
hemisphere).
_Avoid_: lift (for the choice itself), gauge, sign convention

**Lifted element** and **Barred element**:
For an element g, the lifted element is s(g) and the barred element is
Ē·s(g), written ḡ and labelled with the prefix `~`. Which of the two
preimages is barred is a labelling matter and may be swapped per element by
a label convention.
_Avoid_: tilde element, extra element, doubled element

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
always the spin factor system, taking values ±1, induced by the double group
through the section: s(g) s(h) = ω(g, h) s(gh).
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
An irreducible representation over the complex numbers, ordinary or
projective, of a point group. Each irrep has a Mulliken symbol as its name
and may carry a Bethe symbol as an alias. Complex-conjugate one-dimensional
irreps are distinct irreps, never merged into a real pair.
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
constant on a conjugacy class. Classes are ordered and represented by their
first member in the canonical element order.
_Avoid_: class (alone, when ambiguous)

**Class label**:
The name of a conjugacy class: its size followed by the reduced canonical
label of its representative (`2C3`, `3sigma_v`, `6C2'`, `4C3^-1`), and for
double groups `2~C3` or `C2+~C2`.
_Avoid_: class symbol, class header

**Representative**:
The first member of a conjugacy class in the canonical element order; for
a class of the double group, the lifted or barred image of the
representative of the class it covers.
_Avoid_: class leader, generator (of a class)

**Regular class**:
A conjugacy class of a point group whose preimage in the double group is
two classes (it **splits**) rather than one (it is **merged**). Projective
characters can be non-zero only on regular classes, and the number of
regular classes equals the number of projective irreps.
_Avoid_: ω-regular (as the term itself), spinor class, good class

**Frobenius–Schur indicator**:
The value +1, 0 or −1 attached to an irrep stating whether it is real,
complex, or quaternionic (pseudo-real).
_Avoid_: reality, realness type

### Conventions

**Canonical label**:
The convention-free, ASCII label of an element (`C6^1`, `C2'(2)`, `S4^-1`,
`sigma_d(1)`, `C3^-1(4)`), derived by a fixed formula from the abstract group
alone and present in every output.
_Avoid_: default label, internal name, structural label

**Reference label**:
The name a reference work (Altmann & Herzig, Bradley & Cracknell, Bilbao)
gives an element, produced from the canonical label by a label convention and
emitted only when that convention is selected.
_Avoid_: display label, alias, literature label

**Label convention**:
A named, selectable set of rewrite rules and per-group overrides that turns
canonical labels into the reference labels of one reference work. A label
convention implies an orientation of the group in space, but the tool itself
never fixes coordinates.
_Avoid_: setting convention, orientation, standard setting, naming scheme

**Master group**:
The dihedral group `D_Mh` (`M` even) that contains an axial point group and
whose secondary directions a reference work numbers once; the reference's
labels of every axial group restrict from its master (Altmann & Herzig's
Rule 2). A label convention's direction table is keyed by the master's
`M`.
_Avoid_: parent group, supergroup, holohedry

**Basis convention**:
A named, selectable data file (rules for the axial families, generator
matrices for the polyhedral families) that fixes the basis in which irrep
matrices are written, as one reference work writes them. Irrep matrices are
only defined up to unitary equivalence, so a table is meaningless without
naming its basis convention. A convention declares its coverage and may
cover only some families; a group outside it is an error, never a silent
fallback to the canonical basis.
_Avoid_: phase convention, gauge

**Spin-harmonic rule**:
The axial rule of a basis convention that writes every two-dimensional
irrep as a reference work does in a basis of two spherical or spin harmonics
`|j ±m⟩` (Altmann & Herzig): a closed form in the axis order and the irrep's
exponent that fixes which partner comes first (by the chain of subgroups
obtained by dividing the axis order by its smallest prime factor at each
step) and the phase of the secondary element, its only free parameter being
the azimuth of the secondary element the label convention anchors. Verified
against every symmetrized basis the reference prints (`n ≤ 10`).
_Avoid_: harmonic basis (ambiguous with the canonical basis), Wigner
convention, book basis

**Coverage**:
Whether a reference convention's output for a group is verified against the
reference's printed tables (`reference`) or produced by the convention's
rules beyond what the reference prints (`extrapolated`). A property of a
group under one convention, reported in the document, never a silent
fallback: an uncovered group is an error, an extrapolated one is labelled.
_Avoid_: supported, validated, in range

**Canonical basis**:
The convention-free basis of every irrep, used by default: the principal
element is diagonal with its eigenvalue exponents descending, the remaining
phases fixed by formula for the axial families and by shipped data for the
polyhedral families.
_Avoid_: default basis, standard basis, natural basis

**Spinor map**:
The one fixed homomorphism from unit quaternions to SU(2) that turns a
lifted rotation into its two-by-two matrix, with the positive exponent first
on the principal axis.
_Avoid_: SU(2) representation (ambiguous), Cayley–Klein parameters

**Mulliken symbol**:
The name of an irrep in the A/B/E/T/G/H/I notation, extended to arbitrary
axis order (`E_k` for ordinary, `E_{k/2}` for projective irreps). Like
element labels it has two layers: the **canonical Mulliken symbol**, an
ASCII name derived by fixed rules from characters on the principal and
secondary elements and always present, and the **reference irrep label**,
the spelling a reference work uses, produced by a label convention.
_Avoid_: species, symmetry type, F (for three-dimensional irreps)

**Bethe symbol**:
The Γ_i name of an irrep in the numbering of Koster, Dimmock, Wheeler and
Statz, kept only as an alias of the Mulliken symbol and present only for the
32 crystallographic point groups and their double groups.
_Avoid_: Koster symbol, BSW symbol, Γ number as an ordering key

**Principal element**:
The element `(1, 0)` of the rotation group `R` of an axial group, with its
parity: `C_n` for most families, `S_2n` for `S2n` and `Dnd` with even `n`,
`S_n^-1` for `Cnh` and `Dnh` with odd `n`. Its order `N` and its characters
decide A versus B, the E subscript `k`, and the projective subscript `k/2`.
_Avoid_: principal axis (a geometric notion), generator (ambiguous)

**Polyhedral principal element**:
The element of a polyhedral group whose character the symbol rules read:
the element of the principal order (`3` for `T`, `4` for `O`, `5` for `I`)
with the smallest canonical index, with its parity: `C3^1(1)` in `T` and
`Th`, `C4^1(1)` in `O` and `Oh`, `S4^-1(1)` in `Td`, `C5^1(1)` in `I` and
`Ih`. Its character decides `1E`/`2E` and `1E3/2`/`2E3/2` of `T`, the
subscripts `1`/`2` of `A` and `T` and `E1/2`/`E5/2` of `O`, and `T1`/`T2`
and `E1/2`/`E7/2` of `I`; the letter follows the dimension alone.
_Avoid_: principal axis, generator (the basis data's generator pair is a
different, though overlapping, choice)

**Secondary element**:
The element `(0, 1)` of a dihedral rotation group `R = D_N`: `C2'(1)` in the
`Dn`, `Dnh`, `Dnd` families and `sigma_v(1)` in `Cnv`. The sign of its
character decides the subscript 1 or 2 of an A or B irrep.
_Avoid_: secondary axis (a geometric notion)

**Partner**:
The other member of a pair of complex-conjugate one-dimensional irreps (the
`1E`/`2E` split). Both members are separate irreps in every table; `1E`
carries the positive exponent on the principal element.
_Avoid_: real E, degenerate pair, physically irreducible representation

### Output

**Document**:
The single output of one invocation of a table subcommand: the tables for
one point group, under one set of options, with one metadata block and one
field. Deterministic: the same arguments and tool version give the same
bytes.
_Avoid_: file (it may go to standard output), result, dump

**Table**:
One named object inside a document: the element list, the Cayley table, the
factor table, the character table or the irrep matrix table. `all` is the
union of the five; the element list is present in every document.
_Avoid_: section, sheet

**Record**:
The structured object emitted for one element, one conjugacy class or one
irrep (index, label or symbol, and the decomposed fields a consumer would
otherwise parse from the label). Tables refer to records by index only.
_Avoid_: entry, row (for a record), header (for an irrep record)

### Verification

**Self-consistency check**:
An exact property of the tables of one group that holds by mathematics
alone (group axioms, cocycle condition, orthogonality, unitarity, trace
equals character). Proves that a table is correct, not that it follows any
reference.
_Avoid_: sanity check, self-test, assertion

**Reference comparison**:
A test of the tables of one group against a transcribed reference work.
Proves that names, class memberships, factor-system signs and bases agree
with the literature. Comes in two modes: **structural comparison**, which
matches classes by their invariants (size, parity, element order) and
irreps by their characters without using any label, and **labelled
comparison**, which matches rows and columns by the labels of a label
convention.
_Avoid_: validation (ambiguous), cross-check

**Reference data**:
A transcription of a reference work (character table, element labels,
multiplication table, irrep matrices, Bethe symbols) kept as test data with
its source cited; never shipped, and distinct from the convention data
derived from it.
_Avoid_: fixtures (for these files), ground truth

**Check ledger**:
The single table in the verification chapter that lists every check with
an identifier (`AG-1`, `CH-5`, `REF-3`, ...), its oracle, the groups it runs
on and whether the self-check flag runs it.
_Avoid_: test plan, test matrix

**Self-check**:
The optional run of the self-consistency checks on the document a single
invocation is about to emit, selected with `--verify`; a failure reports
every failing check by its ledger identifier and emits nothing.
_Avoid_: verification mode, validation flag

**Golden document**:
A checked-in document for a fixed group, subcommand and format that the
tool must reproduce byte for byte; the regression and determinism oracle.
_Avoid_: snapshot, fixture

**Mutation test**:
A test that corrupts one entry of a document model on purpose and requires
the named checks to fail; the proof that the checks detect what they claim.
_Avoid_: negative test, fault injection

### Number representation

**Field**:
The cyclotomic field `Q(ζ_M)` in which every value of a document lives,
declared once per document.
_Avoid_: number field (unqualified), domain

**Conductor**:
The integer `M` naming a document's field, fixed per group and basis
convention as the smallest `M` that contains every irrep matrix entry.
_Avoid_: field order, root order

**Exact value**:
A number written as its rational coefficient vector over the power basis of
the document's field, never rounded.
_Avoid_: symbolic value, analytic value, coefficient list

**Numeric value**:
The correctly rounded double-precision complex approximation `[re, im]` of an
exact value, emitted next to it for consumers that cannot use exact values.
_Avoid_: float, approximate value

**Pretty value**:
A simplified human-readable string form of an exact value (`"-1 - zeta"`),
optional, never used for comparison or reconstruction.
_Avoid_: symbolic string, display value
