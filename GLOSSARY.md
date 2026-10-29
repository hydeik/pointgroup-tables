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

### Group chains and coupling

**Group chain**:
A descending sequence of point groups `G ⊃ H ⊃ …`, each a subgroup of the
previous one, named by Schoenflies symbols and resolved to concrete
subgroups of the abstract group by a fixed rule. A chain is an input that
selects a basis; it never fixes coordinates.
_Avoid_: subgroup chain, chain of subgroups (when the sequence is meant),
coupling scheme

**Canonical chain**:
A group chain in which every restriction from one group to the next is
multiplicity-free for both ordinary and projective irreps, so that the
chain labels of a component are unique. Only canonical chains are
accepted; a chain with a repeated restriction is an error, never resolved
by an extra multiplicity label.
_Avoid_: complete chain, good chain, proper chain

**Step**:
One group of a chain below the top group `G`, resolved to one concrete
subgroup of the previous group by the embedding rule.
_Avoid_: level, link, member

**Terminal**:
The last group of a chain. An accepted chain ends at a group of the cyclic
shape (`Cn`, `Cnh`, `S2n`), whose double group is abelian, so that the
chain labels name every component.
_Avoid_: bottom group, leaf, end of the chain

**Class anchor**:
An optional list of class labels of `G` attached to a step's symbol
(`C2[C2']`, `C1h[sigma_h]`, `C2v[C2,sigma_h]`), requiring the resolved
subgroup to contain an element of every listed class; the only way a user
steers the embedding rule.
_Avoid_: qualifier, generator, selector

**Embedding rule**:
The fixed rule that resolves a step's symbol to one subgroup and one
identification of its elements: among the admissible isomorphisms of the
symbol's group into the previous group (preserving parity, label
invariants and lifts), those meeting the class anchors and keeping the
chain canonical, the one with the lexicographically smallest image in
canonical indices, the images of the elements of highest order compared
first.
_Avoid_: subgroup selection, orientation rule, conjugate choice

**Canonical spelling**:
The spelling of a chain that carries a class anchor exactly where another
accepted chain class with the same symbols has a different class-label
signature at that step, the anchor being the shortest list of class labels
that excludes those classes, so that it reproduces the chain independently
of the unanchored defaults. A canonical spelling names one chain class
unless that class has twins.
_Avoid_: normalised chain, fully qualified chain

**Class-label signature**:
The set of class labels of `G` whose classes a subgroup meets; the
information a class anchor can see about a subgroup. Two subgroups with the
same signature cannot be told apart by anchors.
_Avoid_: class content, class fingerprint

**Chain class**:
A set of chains of `G` conjugate to one another under `G`, that is, equal up
to a simultaneous conjugation of every step. The unit the subgroup
catalogue counts; conjugate chains give bases related by one irrep matrix,
so a chain class is what a chain spelling can at most name.
_Avoid_: chain orbit, equivalence class of chains

**Twin chain classes**:
Chain classes of `G` that are not conjugate but have the same class-label
signature at every step, so that they share one canonical spelling. The
embedding rule reaches exactly one of them (the *reached* twin); the others
are not accepted, since no chain spelling can name them. The twins found
agree on every ordinary irrep and differ on the projective irreps only.
_Avoid_: ambiguous chains, shadow chains, duplicate classes

**Default chain**:
The one chain of every group that `--chain default` selects: the shortest
canonical descent to the largest cyclic-shape subgroup containing the
principal element (`Dn ⊃ Cn`, `Dnh ⊃ Cnh`, `O ⊃ C4`, `I ⊃ D5 ⊃ C5`). It is
never applied unless asked for; the canonical basis stays the default.
_Avoid_: standard chain, natural chain, preferred chain

**Chain-adapted basis**:
The basis of every irrep of `G` in which the irrep matrices of every
subgroup of the chain are block-diagonal, each component carrying its chain
labels, with ordering and phases fixed by the tool's own deterministic rule
(the eigenfunction method computed exactly). A chain-adapted basis is the
tool's own basis selected by a chain, not a basis convention (which follows
a reference work).
_Avoid_: symmetry-adapted basis, subgroup-adapted basis, CSCO basis,
reduced basis

**Chain label**:
The sequence of canonical Mulliken symbols, one per group of the chain below
`G`, that names a component of an irrep in the chain-adapted basis (the
"quantum numbers of the group chain" of Chen and co-workers).
_Avoid_: quantum numbers, row label, component index (which is the position
in the matrix)

**Label space**:
The subspace of an irrep that one chain label selects, the common range of
the isotypic projectors of the chain's steps; one-dimensional for every
label of an accepted chain.
_Avoid_: eigenspace (unqualified), weight space, sector

**Block**:
The components of an irrep that share the first chain label, on which the
restriction to the first step acts as that step's own chain-adapted
matrices.
_Avoid_: subduction block, sub-block, irreducible part

**Exponent order**:
The order of the labels of a cyclic-shape step in a chain-adapted basis:
by the exponent of each label's eigenvalue at that step's own principal
element, descending, the way version 1 orders the canonical components of
an irrep. A step of the dihedral or polyhedral shape keeps its canonical
irrep order; across steps the chain labels are compared step by step.
_Avoid_: lexicographic order (the comparison across steps, not within
one), irrep order (for a cyclic step), eigenvalue order

**Norm-factor phase rule**:
The one primitive that fixes a phase: the quantity a rule designates is set
equal to the conjugate of the canonical norm factor of its squared modulus,
which is the real positive square root whenever that root lies in the
document field.
_Avoid_: phase convention (unqualified), real-positive rule, Condon–Shortley
convention

**Canonical norm factor**:
The element `μ` of the document field with `μ · conj(μ) = N` that the
arithmetic rule assigns to a totally positive `N`; the real positive square
root when the field contains it.
_Avoid_: normalisation constant, square root of the norm

**Designated component**:
The canonical component of the first vector of a block to which the
norm-factor phase rule is applied between blocks: the first non-vanishing
one fixed by the principal element of the group, or the first
non-vanishing one of all when none is fixed.
_Avoid_: reference component, pivot, leading coefficient

**Clebsch–Gordan coefficient**:
An entry of the unitary matrix that reduces the direct product of two irreps
of `G` (ordinary or projective, with the product factor system) into irreps
of `G`, written in the basis of the document. The coefficients depend on the
basis, so a table is meaningless without the basis the document names.
_Avoid_: coupling coefficient, vector-coupling coefficient, 3-j symbol,
Wigner coefficient

**Irrep pair**:
Two irreps of `G` whose direct product is coupled, written `(α, β)` with
`α ≤ β` in the canonical irrep order. A `cg` document prints every such
unordered pair once; the reversed pair follows from the exchange relation.
_Avoid_: ordered pair, product pair, coupling (which names the operation)

**Product basis**:
The basis `|α i⟩|β j⟩` of the direct product of two irreps, ordered with
the first factor's component outer and the second's inner (row-major); the
row order of every coupling matrix.
_Avoid_: uncoupled basis, tensor basis, Kronecker basis

**Coupling matrix**:
The unitary matrix of one irrep pair whose rows are the product basis and
whose columns are the coupled vectors, ordered by constituent, coupled copy
and component; the object a `cg` document prints per pair.
_Avoid_: CG matrix, transformation matrix, reduction matrix

**Constituent**:
An irrep that occurs, with a positive coupling multiplicity, in the direct
product of two irreps. A `cg` document lists the constituents of every
pair and no irrep that does not occur.
_Avoid_: component (which is a basis vector of an irrep), term, product
irrep

**Isotypic projector**:
The orthogonal projector onto all coupled copies of one constituent in a
direct product, written in the product basis. It does not depend on the
phases of the coupled vectors or on how the copies are mixed, so it is the
quantity two Clebsch–Gordan tables of the same pair must agree on exactly,
whatever their conventions.
_Avoid_: projection operator (unqualified), multiplicity space, isotypic
component (which is the subspace, not the projector)

**Monomial intertwiner**:
A unitary matrix relating two realisations of one irrep that is a
permutation matrix times a diagonal matrix of phases. Two chain-adapted
bases of the same chain are related by one, the permutation being the map
between their chain labels; it is the form in which a reference table of
chain-adapted matrices is compared with the tool's.
_Avoid_: phase matrix (which omits the permutation), relabelling,
basis change (unqualified)

**Coupling multiplicity**:
The number of times an irrep occurs in the direct product of two irreps.
When it exceeds one, the copies are distinguished by a fixed rule of the
tool, never left to the consumer.
_Avoid_: outer multiplicity, repetition index, branching multiplicity (which
is the restriction to a subgroup)

**Coupled copy**:
One of the `n` orthonormal sets of coupled vectors of an irrep that occurs
`n` times in a direct product, numbered by a fixed rule of the tool
(symmetric copies before antisymmetric ones when the two factors are the
same irrep, then the order of their trial product states).
_Avoid_: multiplicity index, repetition, occurrence

**Trial product state**:
The product basis vector from which a coupled copy is projected: the first
one, in the product order of the document's basis, that the projector of
the copy does not annihilate. It is also the first product state whose
coefficient in the designated component of the copy is non-zero.
_Avoid_: seed vector, starting vector, intrinsic quantum number (Chen's
label for the same choice)

**Designated coefficient**:
The coefficient of the trial product state in the first component of a
coupled copy; the one coefficient the phase rule fixes, every other
coefficient of the copy following from the irrep matrices. It is real and
positive whenever the field of the document contains the square root of
its squared modulus.
_Avoid_: leading coefficient, normalisation constant, phase convention

**Exchange sign**:
The sign `+1` or `−1` by which a coupled copy of the product of an irrep
with itself changes when the two factors are exchanged; the symmetric
copies carry `+1`. Defined only when the two factors are the same irrep.
_Avoid_: permutation symmetry, symmetry type, parity (which names
improper elements)

**Squared norm**:
The totally positive number `N = x · conj(x)` of the quantity a phase rule
designates: the diagonal entry of the projector that produces a
chain-adapted or coupled vector, and the only source of new numbers beyond
the irrep matrices.
_Avoid_: norm (unqualified, which may mean the square root), normalisation
constant

**Hermitian norm**:
A number of the form `μ · conj(μ)` for some `μ` of the document field; a
squared norm is one exactly when a phase puts the vector into the field.
_Avoid_: norm from the field, relative norm (unqualified)

**Chain field**:
The smallest cyclotomic field holding the irrep matrices of the group and
the characters of every step of a chain, the field in which the projectors
of a chain live.
_Avoid_: subgroup field, projector field

**Chain conductor**:
The conductor a document produced with `--chain` declares: the smallest
enlargement of the chain field in which every chain-adapted matrix of the
group can be written, fixed by the group and the chain class alone.
_Avoid_: chain field (which may be smaller), enlarged conductor

**Coset decomposition**:
The partition of one group of a chain into the left cosets of the next
step, each coset named by the member with the smallest canonical index; the
data from which symmetry-adapted functions are built by the consumer.
_Avoid_: right cosets, transversal, coset representatives (which are only
the first members)

**Double coset**:
A union of left cosets of a step that is closed under multiplication by the
step on both sides; the double cosets of one step are the generator set of
the double-induced technique of the eigenfunction method.
_Avoid_: double-coset class, orbit

**Branching rule**:
For one step of a chain, which irreps of the step occur in the restriction
of each irrep of the previous group, ordinary irreps to ordinary and
projective to projective; along an accepted chain no irrep occurs twice.
Also called subduction.
_Avoid_: correlation table, compatibility relation, decomposition
(unqualified)

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

**Subgroup record**:
The record emitted for one step of a chain: the step's group, its embedding
into the top group, its coset and double-coset decompositions, its irreps
and the branching rule of the step. Positional: a step is never referred to
by an index.
_Avoid_: step record, chain entry

**Component record**:
The record emitted for one component of an irrep in a chain-adapted basis:
its position and its chain label, as irrep indices of the steps and as
symbols. Absent without a chain.
_Avoid_: label record, basis-vector record

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

**Test oracle**:
A closed form from the literature implemented inside the test suite to
produce the expected values of a check for an arbitrary axis order; a
source of reference values that is code rather than a transcription, and
never part of the shipped tool.
_Avoid_: reference implementation, model answer, golden formula

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
