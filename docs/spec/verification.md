# Verification

Decided in wayfinder ticket
[Verification strategy: self-consistency tests and reference-table comparison](https://github.com/hydeik/pointgroup-tables/issues/12).
Every other chapter ends with a "Verification" section listing the
properties its tables must have. This chapter is the index of all of them:
it fixes where each check runs, on which groups, against which oracle, and
what it takes for the specification to count as implementable. A check is
defined once, in the chapter that owns the table; this chapter assigns it an
identifier and a place in the test suite. Rationale for redistributing the
reference data in ADR-0006.

## Two kinds of proof

An implementation proves its tables in two independent ways:

1. **Self-consistency checks** prove that the tables are *mathematically
   correct*: the Cayley table is a group, the factor system is a normalised
   cocycle with values `±1`, the characters are orthogonal, the matrices
   are unitary homomorphisms with the right traces, and so on. They need no
   external data and hold for every group the tool can construct, including
   an axis order no book tabulates. Every self-consistency check is exact
   (ADR-0002); there is no tolerance anywhere.
2. **Reference comparison** proves that the tables *agree with the
   literature*: that the names, the class memberships, the sign conventions
   of the factor system and the bases of the matrices are the ones the
   reference works print. Self-consistency cannot detect a wrongly named
   irrep or a barred element assigned to the wrong class; only a
   transcribed reference can.

Self-consistency is the proof of correctness; reference comparison is the
proof of convention.

## Where checks run

The **test suite** (`pytest`, under `tests/`) is the home of every check.
It is the proof that ships with a version: the implementation is correct
for every group the suite covers.

In addition, every table subcommand accepts the flag **`--verify`**
([Command line](./cli.md#options)). It runs the self-consistency checks
that apply to the document about to be emitted, for the group and options
actually requested, and only emits the document when all of them pass. It
exists because the axis order `n` is unbounded: the test suite covers a
finite set of groups, and a user who generates `D37h` has `--verify` to
prove that document. `--verify` is off by default; its cost is roughly
`O(|G†|² · Σ d³)` and reaches minutes for large `n`.

`--verify` runs:

- the **base checks** for every table subcommand: the group checks
  (`AG-*`), the class checks (`CL-*`), and, whenever the double group is
  involved (`--double-group`, or a document that contains projective
  irreps), the factor-system checks (`DG-*`);
- the **table checks** for each table in the document: `CH-*` for
  `characters`, `MX-*` for `matrices`, `VL-*` for every table that carries
  values. `IL-*` runs with `characters` and `matrices`.

It never runs a reference comparison (reference data is test data, not
shipped), nor the output-schema, CLI, golden or mutation checks.

On failure `--verify` does not stop at the first check: it runs all of them
and reports every failure, one line each, then exits with code `1` and
writes no document ([Command line](./cli.md#exit-codes-and-errors)):

```
pointgroup-tables: error: verification of 'D3' failed (2 checks)
verify: CH-1: rows 1 and 2 are not orthogonal: <D(1), D(2)> = 6
verify: MX-3: D(C3^1) of irrep 2 is not unitary
```

The identifier after `verify:` is the check identifier of the ledger below.
Every checker is a pure function of the computed document model (the
tables as they will be serialised), never of internal state, so the same
function serves `--verify`, the test suite and the mutation tests.

## Group sets

The test suite is parametrised over named sets of groups. Aliases
(`Cs = C1h`, `D1h = C2v`, `S2 = Ci`, ...) are resolved by the group
argument grammar and counted once.

| Set | Contents | Marker |
| --- | -------- | ------ |
| `seed` | the 43 groups of the seed table in [Double groups](./double-groups.md#verification-data): the 32 crystallographic groups, `I`, `Ih`, `C5`, `C5v`, `D5h`, `D5d`, `C7v`, `D7h`, `D8d`, `C8h` | none |
| `sweep` | every axial family `Cn`, `Cnv`, `Cnh`, `S2n`, `Dn`, `Dnd`, `Dnh` for `n = 1 ... 12` | none |
| `large` | every axial family for `n ∈ {13, 16, 17, 24, 30}` | `slow` |
| `default` | `seed ∪ sweep` | none |

`sweep` exhausts the case distinctions of the labelling and basis rules
(odd and even `n`, `n` divisible by 4 or not, the degenerate `n = 1, 2`
families); `large` proves that the conductor formulas and the `E_k`
subscripts are general and not an artefact of the crystallographic range,
with `13` and `17` prime and `16`, `24`, `30` highly composite.

A check that names no set runs on `default`.

## Exhaustive, never sampled

Every check is exhaustive over its domain: all pairs for the group axioms
and the homomorphism check, all triples for the cocycle condition, every
element of every group in the set. The "sampled beyond `|G| ≤ 48`" clauses
of the abstract-group and double-group chapters are superseded by this
chapter: a sampled check is not a proof, and the cost is bounded (the
largest group of `default` is `Ih†` with 240 elements; `D30h†` in `large`
has 240 as well). Checks that are slow on the largest groups are not
weakened but scheduled: they carry the `slow` marker for groups with
`|G†| > 96`.

## Check ledger

The ledger is the complete list of checks. An implementation session turns
it into the parameter table of the test suite one row at a time; each row
names the chapter that defines the check, the oracle it uses, and the
group set and marker it runs under.

Oracles: **self** (a property of the tables alone), **seed** (the
checked-in tables of known values), **reference** (transcribed reference
data under `tests/reference/`), **schema** (`schema.json`), **golden**
(checked-in documents).

Markers: none (the default battery), `slow`, `reference` (skipped with a
reason when the reference file it needs is absent). The `--verify` column
says whether the flag runs the check.

### Groups and classes

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `AG-1` | [Abstract groups](./abstract-groups.md#verification) | Cayley table of `G` satisfies closure, identity, inverses, associativity (all triples) | self | `default`, `large` | `slow` for `large` | yes |
| `AG-2` | Abstract groups | `(order, class count, histogram of element orders)` equals the seed table | seed | `seed` | | yes (`seed` groups only) |
| `AG-3` | Abstract groups | canonical element order is proper-then-improper with contiguous classes | self | `default`, `large` | | yes |
| `DG-1` | [Double groups](./double-groups.md#factor-system) | `ω(E, g) = ω(g, E) = 1` and every value is `1` or `−1` | self | `default`, `large` | | yes |
| `DG-2` | Double groups | cocycle condition over all triples | self | `default`, `large` | `slow` for `\|G†\| > 96` | yes |
| `DG-3` | Double groups | `ω(g, g⁻¹) = ω(g⁻¹, g)` | self | `default`, `large` | | yes |
| `DG-4` | Double groups | Cayley table of `G†` satisfies the group axioms; `\|G†\| = 2\|G\|`; `~E` is central of order 2 | self | `default`, `large` | `slow` for `large` | yes |
| `DG-5` | Double groups | `ω` depends on rotation parts only; for a type II group it is the pull-back of `ω` of `G⁺`; `ω(C2, C2) = −1` for every binary rotation | self | `default`, `large` | | yes |
| `DG-6` | Double groups | class counts of `G†` and the number of projective irreps equal the seed table | seed | `seed` | | yes (`seed` groups only) |
| `CL-1` | [Conjugacy classes](./conjugacy-classes.md#verification) | classes partition the elements; sizes divide `\|G\|`; `size × \|centraliser\| = \|G\|` | self | `default`, `large` | | yes |
| `CL-2` | Conjugacy classes | class labels match the grammar and are pairwise distinct | self | `default`, `large` | | yes |
| `CL-3` | Conjugacy classes | `regular` classes are exactly the split ones; their number equals the number of projective irreps | self | `default`, `large` | | yes |
| `CL-4` | Conjugacy classes | `inverse` is an involution; `base_class`, `split`, `barred` agree with the element records | self | `default`, `large` | | yes |
| `CL-5` | Conjugacy classes | polyhedral class labels reproduce the structures listed in [Label conventions](./label-conventions.md) | self | polyhedral groups | | yes |

### Irrep labels and characters

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `IL-1` | [Irrep labels](./irrep-labels.md#verification) | exactly one canonical symbol per irrep, grammar-conformant, pairwise distinct | self | `default`, `large` | | yes |
| `IL-2` | Irrep labels | the symbol is a function of the character (recomputed from the character table it is reproduced) | self | `default`, `large` | | yes |
| `IL-3` | Irrep labels | `dimension = χ(E)`; `frobenius_schur = (1/\|G\|) Σ ω(g, g) χ(g²)`; `partner` present iff the indicator is `0` and the conjugate is a different irrep | self | `default`, `large` | | yes |
| `IL-4` | Irrep labels | canonical irrep order is the one the chapter fixes | self | `default`, `large` | | yes |
| `CH-1` | this chapter | row orthogonality: `Σ_g χ_i(g) χ_j(g)* = \|G\| δ_ij` over `G` (ordinary) and over `G†` (all irreps of `G†`) | self | `default`, `large` | | yes |
| `CH-2` | this chapter | column orthogonality: `Σ_i χ_i(C) χ_i(C')* = \|centraliser\| δ_CC'` | self | `default`, `large` | | yes |
| `CH-3` | this chapter | `Σ d² = \|G\|` over the ordinary irreps and again over the projective irreps; the number of ordinary irreps equals the class count of `G` | self | `default`, `large` | | yes |
| `CH-4` | this chapter | characters are class functions of `G†`; `χ(~g) = +χ(g)` for ordinary and `−χ(g)` for projective irreps; projective characters vanish off regular classes | self | `default`, `large` | | yes |
| `CH-5` | this chapter | Frobenius–Schur consistency: the indicator is `±1` iff the character is real, `0` iff it is not; `Σ_i ν_i d_i = #{x ∈ G† : x² = E}` over all irreps of `G†` | self | `default`, `large` | | yes |
| `CH-6` | this chapter | characters equal those computed independently by the Burnside class-matrix method | self | `default`, `large` | `slow` for `\|G†\| > 96` | no |
| `VL-1` | [Values](./values.md) | every exact value lies in the declared field and its coefficient vector has length `φ(M)` | self | `default`, `large` | | yes |
| `VL-2` | Values | every numeric value is the correctly rounded double of its exact value | self | `default` | | yes |
| `VL-3` | Values | the declared conductor is minimal: no proper divisor `M'` of `M` contains every value of the document | self | `default`, `large` | `slow` for `large` | no |

### Irrep matrices

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `MX-1` | [Basis conventions](./basis-conventions.md#verification) | homomorphism by closure: every element reached by two words gets the same matrix (up to `ω` for projective irreps) | self | `default`, `large` | | yes |
| `MX-2` | Basis conventions | `D(g) D(h) = ω(g, h) D(gh)` for all pairs and all irreps | self | `default`, `large` | `slow` for `\|G†\| > 96` | yes |
| `MX-3` | Basis conventions | every matrix is unitary under `ζ ↦ ζ⁻¹` | self | `default`, `large` | | yes |
| `MX-4` | Basis conventions | trace of every matrix equals the character | self | `default`, `large` | | yes |
| `MX-5` | Basis conventions | canonical basis: `D(c)` is diagonal with the prescribed eigenvalue order | self | `default`, `large` | | yes |
| `MX-6` | Basis conventions | reference basis convention: the intertwiner to the canonical realisation exists and is unique up to a phase | self | groups the convention covers | | yes |
| `MX-7` | Basis conventions | every entry lies in the declared field | self | `default`, `large` | | yes |
| `MX-8` | Basis conventions | irreps are pairwise inequivalent and irreducible: `(1/\|G\|) Σ_g \|χ(g)\|² = 1` | self | `default`, `large` | | yes |

### Documents and command line

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `OS-1` | [Output schema](./output-schema.md#verification) | every document of every subcommand validates against `schema.json`, for every combination of `--double-group`, `--irreps`, `--values` | schema | `seed` | | no |
| `OS-2` | Output schema | `tables` and `values` list exactly what is present; all forms of a table have the same shape | self | `seed` | | no |
| `OS-3` | Output schema | every index is in range for its index space; `items[i].index = i`, `classes[c].index = c`; canonical irrep order in `characters` and `matrices` | self | `seed` | | no |
| `OS-4` | Output schema | irrep records of `characters` and `matrices` agree in `all`; `characters.classes` equals the classes of the matching `elements` document | self | `seed` | | no |
| `OS-5` | Output schema, [Command line](./cli.md#verification) | `all` is the union of the five single-table documents; every single-table document embeds the same `elements` | self | `seed` | | no |
| `OS-6` | Output schema, Command line | determinism: two runs with the same arguments are byte-identical, every subcommand and format | self | `seed` | | no |
| `CLI-1` | Command line | the group-argument table of `(exit code, symbol)` | self | as listed | | no |
| `CLI-2` | Command line | unknown convention gives `4`; option on the wrong subcommand gives `2`; on failure the document stream is empty and no file is created | self | | | no |
| `CLI-3` | Command line | `--verify` on a passing group emits the same bytes as without it; on `list`/`schema` gives `2`; on a failing check gives `1`, no document, one `verify:` line per failing check | self, mutation | `C3v`, `Oh` | | no |
| `GOLD-1` | this chapter | golden documents: every subcommand × `json`, `yaml`, `toml` for `C3v`, `D6h`, `Oh`, `Ih` with default options equals the checked-in file byte for byte | golden | those four | | no |

### Reference comparison

| ID | Chapter | Check | Oracle | Groups | Marker | `--verify` |
| -- | ------- | ----- | ------ | ------ | ------ | ---------- |
| `REF-1` | this chapter | structural comparison of the character tables of `G` and `G†` with Altmann & Herzig | reference | groups transcribed | `reference` | no |
| `REF-2` | this chapter | structural comparison of the character tables of `G` and `G†` with Bilbao | reference | 32 crystallographic groups | `reference` | no |
| `REF-3` | this chapter | labelled comparison of the character tables with Altmann & Herzig under the `altmann-herzig` label convention | reference | groups transcribed | `reference` | no |
| `REF-4` | [Double groups](./double-groups.md#verification-data) | `ω(g, h) = −1` exactly where the Altmann & Herzig multiplication table (the `[multiplication]` table of the labels file) marks the product with a tilde, after the convention's `barred` swaps | reference | `C2v`, `D3`, `D4`, `D3h`, `T`, `O`, `Ih` | `reference` | no |
| `REF-5` | Irrep labels | canonical symbols after the `altmann-herzig` convention equal the book's symbols row by row | reference | groups transcribed | `reference` | no |
| `REF-6` | Label conventions | the `altmann-herzig` label convention reproduces the transcribed element labels of every group, element by element | reference | groups transcribed | `reference` | no |
| `REF-7` | Basis conventions | under `--basis-convention altmann-herzig`, the matrices of the transcribed generators equal the book's entry by entry (polyhedral groups from the polyhedral transcription, axial samples from the character-table transcription) | reference | groups transcribed | `reference` | no |
| `REF-8` | Irrep labels | the shipped Koster table: each `Γ_i` entry's identifying characters match the computed irrep it is attached to | reference | 32 crystallographic groups | `reference` | no |

### Mutation tests

| ID | Check | Oracle | Marker |
| -- | ----- | ------ | ------ |
| `MUT-1` | one Cayley entry replaced by another element is caught by `AG-1` | mutation | |
| `MUT-2` | one `ω(g, h)` sign flipped is caught by `DG-2` and by `MX-2` | mutation | |
| `MUT-3` | one matrix multiplied by a phase `ζ` is caught by `MX-1` or `MX-2` and `MX-4` | mutation | |
| `MUT-4` | one character replaced by another value of the field is caught by `CH-1` and `CH-2` | mutation | |
| `MUT-5` | two class labels swapped is caught by `CL-4` or `CL-2`; two irrep symbols swapped by `IL-2` | mutation | |
| `MUT-6` | one element moved to another class is caught by `CL-1` and `CH-4` | mutation | |
| `MUT-7` | one row of a reference character table altered is caught by `REF-1` | mutation | `reference` |

A mutation test builds the document model of `C3v` and `Oh` (one axial, one
polyhedral group, both with projective irreps), applies the mutation to the
model, and runs the checkers on the corrupted model. Every check with
`--verify = yes` in the ledger must be the target of at least one mutation;
the suite fails otherwise. This is the proof that the proof works.

## Reference comparison

### Structural comparison

A **structural comparison** matches a computed character table with a
reference table without using any label. Every reference class carries the
invariants the reference lets one read off: its `size`, its `parity`
(proper or improper), and, when the reference lists its elements, the
`order` of its elements in the group of the table. The comparison succeeds
when there is a bijection of classes that preserves the invariants and a
bijection of irreps such that the two tables agree entry by entry, exactly.
The search is a backtracking over the small candidate sets the invariants
leave; a match must exist, and the test reports the matched bijection so a
human can inspect it.

The structural comparison is available for every reference and needs no
label convention, so it runs as soon as the reference file exists. It
proves the mathematical content of the transcription (and thus also catches
transcription errors, which is why a second, independent reference is
kept), but not that the tool's names agree with the reference's.

### Labelled comparison

A **labelled comparison** applies the selected label convention to the
computed tables, then matches rows to the reference by irrep label and
columns by class label, and requires the entries to agree exactly. It runs
only for a reference that has a label convention (`altmann-herzig` once the
convention data exists) and is what proves the convention: element labels,
irrep symbols, class membership and the barred assignment. Comparing a
double-group table requires the convention's `barred` swaps to be applied
first ([Conjugacy classes](./conjugacy-classes.md)).

### References and what they check

| Reference | Tables | Groups | Mode |
| --------- | ------ | ------ | ---- |
| Altmann & Herzig, *Point-Group Theory Tables*, 2nd ed. | character tables of `G` and `G†`; element labels; irrep symbols; tilde-marked products; irrep matrices of the polyhedral groups plus axial samples | the 32 crystallographic groups, `I`, `Ih`, `C5v`, `D5h`, `D5d`, `S8`, `D4d` | structural, then labelled |
| Bilbao Crystallographic Server, point-group and double-point-group tables | character tables of `G` and `G†` | the 32 crystallographic groups | structural only (no Bilbao label convention in version 1) |
| Koster, Dimmock, Wheeler & Statz | Bethe symbols | the 32 crystallographic groups | by identifying characters (`REF-8`) |

Bradley & Cracknell is not transcribed; it appears only in the notes of the
Koster table where its `Γ_i` numbering disagrees.

## Reference data

Reference data is **test data** and lives under `tests/reference/`, never
in the package. Data the tool ships (label and basis conventions, the
Koster table) lives under `src/pointgroup_tables/data/conventions/` and is
derived from the transcriptions, not identical to them: a convention is
rules plus overrides, a transcription is the book as printed.

```
tests/reference/
  README.md                          # sources, citation, licensing note (ADR-0006)
  altmann-herzig/
    labels/<group>.toml              # element labels, the book's identification of each operation,
                                     # and the book's multiplication table with its tilde marks
    characters/<group>.toml          # character tables of G and G†
    matrices/<group>.toml            # generator matrices (polyhedral groups, axial samples)
  bilbao/
    fetch.py                         # the script that produced the files; run once, not by the tests
    characters/<group>.toml
  koster/
    bethe/<group>.toml               # Γ_i per canonical symbol with identifying characters
```

Every file is TOML, one per group, named by the normalised Schoenflies
symbol, and starts with the header the label transcription established
(`tests/reference/altmann-herzig/labels/README.md`):

```toml
group = "D3"
table = "T 20.1, p. 176"     # where in the source; a URL for Bilbao
status = "done"              # todo | partial | done
notes = ""                   # anything that could not be read unambiguously
```

The tilde-marked products stay inside the labels file (its
`[multiplication]` table), because the isomorphism between the book's
labels and the canonical ones is derived from that table; `REF-4` reads it
from there. Each kind documents its own body in a `README.md` next to the
files.

Values use the value grammar of [Basis conventions](./basis-conventions.md#value-grammar)
(`i`, `zeta`, `sqrt2`, `sqrt3`, `sqrt5`, `phi`), so characters and matrices
are parsed by the same exact parser as the convention files. A character
table file lists the reference's classes with their invariants and then
one row per irrep:

```toml
double = false

[[classes]]
label = "2C3"           # the reference's own spelling
size = 2
order = 3               # omitted when the reference does not identify its elements
parity = 1

[[irreps]]
label = "E"
characters = ["2", "-1", "0"]
bethe = "Gamma3"        # optional
```

A Bilbao file is produced by `fetch.py` from the server's HTML and is then
authoritative: the tests read the TOML, never the network, and a change in
the server's pages is not a change in the data. The script stays in the
repository as the record of provenance.

### Redistribution

The transcribed tables are redistributed with the repository under its
licence (ADR-0006). They are mathematical facts (characters, products,
labels) re-expressed in the project's own format, restricted to what the
tests need, and each file cites its source. `tests/reference/README.md`
carries the citations, the Bilbao terms of use, and a statement that the
files are not reproductions of the books' pages.

## Test organisation

- **Framework**: `pytest`, parametrised over the group sets above. No
  property-based testing library; every domain is a finite set that is
  enumerated.
- **Markers**: `slow` (large groups and the exhaustive pair checks on
  `|G†| > 96`) and `reference` (needs a file under `tests/reference/`;
  skipped with a reason naming the file when it is absent). A plain
  `pytest` runs the default battery without either; continuous integration
  runs `slow` too.
- **Golden documents**: `tests/golden/<group>/<subcommand>.<format>` for
  `C3v`, `D6h`, `Oh`, `Ih`, every table subcommand, all three formats,
  default options. They are the regression and determinism oracle
  (`GOLD-1`). A script under `tools/` regenerates them; a regeneration is
  reviewed as a diff, never accepted blindly. The samples in
  `prototypes/samples/` are their seed.
- **Mutation tests**: as in the ledger; a fixture applies one mutation to
  the document model and asserts the expected check identifiers fail and
  no other check crashes.
- **`--verify` tests**: `CLI-3` runs the flag end to end on a passing group
  and, through an injected mutation, on a failing one.

## Acceptance criterion

The specification is **implementable** when all of the following hold:

1. **The ledger is complete.** Every item of every chapter's "Verification"
   section has a row in the check ledger with an identifier, an oracle, a
   group set and a marker. An implementation session can transcribe the
   ledger into the test suite without inventing a check.
2. **The schema is covered.** Every field of `schema.json` is defined by
   a chapter of `docs/spec/`.
3. **Nothing is deferred.** No chapter contains a `TODO`, `TBD`, or a
   decision deferred to another ticket. References to the open
   transcription tasks (element labels, character tables and axial matrices,
   polyhedral matrices, Koster table, Bilbao fetch) are allowed: the tests
   that need their data carry the `reference` marker and skip until the
   data exists, so the transcriptions are not a precondition of
   implementation.
4. **The seed documents validate.** For the 43 groups of `seed`, the
   documents of every subcommand validate against `schema.json`.

The transcriptions themselves are not on this list on purpose: they are
human work with the books at hand and may finish after the first
implementation; a `reference` test turns from skipped to green or red the
day its file lands.
