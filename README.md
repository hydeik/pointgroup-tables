# pointgroup-tables

Exact character tables, irrep matrices and double-group tables of the finite
point groups, as machine-readable documents.

`pointgroup-tables` computes the tables of any finite point group of the
Schoenflies families `Cn`, `Cnv`, `Cnh`, `S2n`, `Dn`, `Dnd`, `Dnh` (any axis
order `n`) and the polyhedral groups `T`, `Td`, `Th`, `O`, `Oh`, `I`, `Ih`
from an abstract-group construction, with no fixed coordinate system, and
emits one self-describing, deterministic document per invocation in JSON,
YAML or TOML. Every value is exact (a coefficient vector in one cyclotomic
field per document) with a correctly rounded double next to it. The tool has
no runtime dependencies and needs Python 3.13 or later.

```
$ pointgroup-tables characters C3v
{
  "schema_version": 1,
  "generator": {
    "name": "pointgroup-tables",
    "version": "0.1.0"
  },
  "group": {
    "symbol": "C3v",
    "input_symbol": "C3v",
    "family": "Cnv",
    "n": 3,
    "order": 6,
    "type": "twisted",
    "rotation_group": "D3"
  },
  "tables": ["elements", "characters"],
  "double": false,
  "irreps": "all",
  "values": ["exact", "numeric"],
  ...
}
```

## Installation

```
uv tool install pointgroup-tables      # or: pipx install pointgroup-tables
pointgroup-tables --version
```

The package is pure Python with no dependencies; any installer that reads
PyPI works in a Python 3.13 environment.

## What the tables are

| Subcommand | Tables in the document |
| ---------- | ---------------------- |
| `elements` | the element list of `G` (or of the spin double group `G†` with `--double-group`): canonical label, kind, order, inverse and conjugacy class of every element, and the class records |
| `cayley` | the Cayley (multiplication) table of `G` or `G†` |
| `factor` | the factor system `ω(g, h) ∈ {±1}` of the spin double group |
| `characters` | the character table of the ordinary and the projective (double-valued) irreps, with class headers and the Frobenius–Schur indicator of every irrep |
| `matrices` | explicit irrep matrices of every irrep, in the selected basis |
| `all` | the union of the five, under one metadata block |

Every document embeds its `elements` table, so indices in every other table
can be resolved without a second invocation. Elements, classes and irreps
carry convention-free canonical labels (`C3^1`, `sigma_v(2)`, `~C2'(1)`;
`2C3`, `3sigma_v`; canonical Mulliken symbols such as `A1g`, `1E1`,
`E1/2`, `B3/2`, `G3/2`, `I5/2`), and the 32 crystallographic groups also
carry Koster's Bethe symbols (`Gamma4+`).

## Usage

```
pointgroup-tables <subcommand> <group> [options]
pointgroup-tables list groups | list conventions | schema
pointgroup-tables --version | --help
```

The group is named by its Schoenflies symbol, case-insensitively, with the
usual aliases (`Cs`, `Ci`, `V`, `Vh`, `Vd`, `D1h`, `S2`, ...). `n` is
unbounded: `C7v`, `D12h`, `S16` and `D37h` are ordinary arguments. Run
`pointgroup-tables list groups` for the accepted families and aliases.

### Options

Common to every table subcommand:

| Option | Values | Default |
| ------ | ------ | ------- |
| `-f`, `--format` | `json`, `yaml`, `toml` | `json` |
| `-o`, `--output <path>` | a file path; `-` is standard output | standard output |
| `--values <set>` | comma-separated subset of `exact`, `numeric`, `both`, `pretty` | `both` |
| `--label-convention <name>` | a shipped label convention | none (canonical labels only) |
| `--basis-convention <name>` | `canonical` or a shipped basis convention | `canonical` |
| `--verify` | run every self-consistency check on the document before emitting it | off |

Specific to some subcommands (a usage error elsewhere):

| Option | Subcommands | Values | Default |
| ------ | ----------- | ------ | ------- |
| `--double-group` | `elements`, `cayley`, `all` | flag | off |
| `--irreps` | `characters`, `matrices`, `all` | `ordinary`, `projective`, `all` | `all` |

The output format is never inferred from the file name. Every option is
recorded in the document metadata, so a document is self-describing; the same
arguments and the same tool version always produce byte-identical output.

### Examples

```
pointgroup-tables characters C3v
pointgroup-tables characters c3v --format yaml --values exact
pointgroup-tables matrices D6h --basis-convention altmann-herzig -o d6h.json
pointgroup-tables elements Oh --double-group --label-convention altmann-herzig
pointgroup-tables all Ih --irreps projective --values numeric,pretty
pointgroup-tables characters D37h --verify
pointgroup-tables schema -o schema.json
```

### Values

Every table value is emitted in the requested forms, as parallel arrays of
the same shape:

- `exact`: the coefficient vector of the value over the power basis of the
  document's cyclotomic field `Q(ζ_M)`, each coefficient an integer pair
  `[numerator, denominator]`. The field (`conductor`, `degree`, `generator`,
  `basis`) is declared once in the metadata; it is fixed by the group and
  the basis convention.
- `numeric`: `[re, im]` doubles, correctly rounded from a 40-digit
  evaluation, printed as the shortest round-tripping representation.
- `pretty`: a human-readable spelling in the value grammar, such as
  `zeta^7`, `(1+sqrt5)/4` or `-1/2*i*sqrt3`.

`--values both` is `exact` and `numeric`.

### Double groups and projective irreps

`--double-group` lists the elements of `G†` as the elements of `G` followed
by their barred partners (`~C3^1`), in the canonical section (rotation angle
in `(−π, π]`, pole in the positive hemisphere), which fixes
`ω(C2, C2) = −1`. The character table lists ordinary and projective irreps
together by default; its columns are the classes of `G†` unless
`--irreps ordinary` is given, in which case they are the classes of `G`.
Indices in the `G†` index space follow one rule: an index `i ≥ |G|` names
the barred element of element `i − |G|`.

### Conventions and `coverage`

Canonical labels and the canonical basis are always emitted. A **label
convention** adds a reference work's element labels, class labels and irrep
symbols next to the canonical ones; a **basis convention** replaces the
irrep matrices by the reference's (the characters do not change). Both are
shipped as data files, versioned with the tool, and verified against
transcriptions of the reference in the test suite. Run
`pointgroup-tables list conventions` for the shipped names.

| Convention | Kind | Source |
| ---------- | ---- | ------ |
| `altmann-herzig` | label | S. L. Altmann and P. Herzig, *Point-Group Theory Tables*, 2nd ed. (Wien: Clarendon, 2011) |
| `altmann-herzig` | basis | the same book: its polyhedral generator images and its spin-harmonic axial bases |
| `canonical` | basis | the tool's own basis (default) |

The metadata of a document under a reference convention carries
`coverage`: `"reference"` when the convention's output for the group is
verified against the reference's printed tables (for `altmann-herzig`: every
polyhedral group and every axial group with `n ≤ 10`), `"extrapolated"`
when the convention's rules are applied beyond what the reference prints
(`D14h`, `C24v`, ...). A convention that does not cover a group at all is an
error (exit code `4`), never a silent fallback. Labels are output only and
are never accepted as input.

### `--verify`

`--verify` runs the self-consistency checks of the verification ledger
(group axioms, the cocycle condition of the factor system, class structure,
orthogonality of the characters, unitarity and the homomorphism property of
the matrices, the declared field, ...) on the document about to be emitted.
A passing document is emitted byte for byte as without the flag; a failure
reports every failing check by its ledger identifier and emits nothing:

```
pointgroup-tables: error: verification of 'D3' failed (2 checks)
verify: CH-1: rows 1 and 2 are not orthogonal: <D(1), D(2)> = 6
verify: MX-3: D(C3^1) of irrep 2 is not unitary
```

The checks are exhaustive, so a group no book tabulates is proven on the
spot. Their cost grows as `|G†|²` and reaches minutes for large `n`.

### Exit codes

| Code | Meaning |
| ---- | ------- |
| `0` | success |
| `1` | any other failure: output file not writable, internal error, a `--verify` check failing |
| `2` | usage error: unknown subcommand or option, an option on a subcommand that does not take it, a group argument that does not match the grammar (`3v`, `C03`), a missing group |
| `3` | unknown group: the argument matches the grammar but names no group (`Xh`, `Cnd`, `T5`) |
| `4` | unknown convention, or a convention that does not cover the requested group |
| `5` | unsupported `n`: `C0`, `D0h`, `S3` (`S2n` takes even subscripts only) |

Errors go to standard error as one line, `pointgroup-tables: error:
<message>`, optionally followed by `hint: <hint>`. On every failure the
document stream is empty and no output file is created.

### Documents

`pointgroup-tables schema` prints the JSON Schema (draft 2020-12) that every
table document validates against. Metadata keys come first, in a fixed order
(`schema_version`, `generator`, `group`, `tables`, `double`, `irreps`,
`values`, `field`, `label_convention`, `basis_convention`), then the tables
in the order `elements`, `cayley`, `factor`, `characters`, `matrices`. YAML
and TOML documents are plain serialisations of the same object. Indices in
positional arrays are the only identifiers: element index in the canonical
element order, class index in the canonical class order, and global irrep
index in the canonical irrep order (ordinary then projective), which
`--irreps` filters without renumbering.

## Specification and design

The normative text is the specification under [`docs/spec/`](docs/spec/)
(one chapter per module: abstract groups, double groups, values, conjugacy
classes, irrep labels, label conventions, basis conventions, output schema,
verification, command line) and the decision records under
[`docs/adr/`](docs/adr/). The domain vocabulary is fixed in
[`CONTEXT.md`](CONTEXT.md).

## Development

```
uv sync
uv run pytest                 # the default battery
uv run pytest -m ""           # everything, including the large groups (slow)
uv run ruff check . && uv run ruff format --check . && uv run ty check
```

- **Markers.** `slow` (the `large` group set, `n ∈ {13, 16, 17, 24, 30}`,
  and the exhaustive pair checks on `|G†| > 96`) is deselected by default;
  continuous integration runs it. `reference` tests need a transcription
  under `tests/reference/` and skip with a reason naming the file when it
  is absent.
- **Golden documents.** `tests/golden/<group>/<subcommand>.<format>` for
  `C3v`, `D6h`, `Oh`, `Ih`, every table subcommand and format with default
  options, byte for byte (`GOLD-1`). `uv run python tools/golden.py check`
  compares, `regenerate` rewrites; a regeneration is reviewed as a diff.
- **Reference data.** The transcriptions of Altmann & Herzig, the Bilbao
  Crystallographic Server and Koster et al. under `tests/reference/` are
  test data and are never shipped; the shipped conventions are generated
  and proven from them by the scripts under `tools/`. A misprint in a
  reference is recorded as a known discrepancy in its transcription and
  corrected before comparison, never reproduced.
- **Prototypes.** `prototypes/` holds the reference code of the arithmetic,
  the constructions and the document layout; nothing there is shipped.
- **Releases.** See [`docs/release.md`](docs/release.md).

## References

- S. L. Altmann and P. Herzig, *Point-Group Theory Tables*, 2nd ed. (Wien:
  Clarendon, 2011).
- G. F. Koster, J. O. Dimmock, R. G. Wheeler and H. Statz, *Properties of
  the Thirty-Two Point Groups* (MIT Press, 1963).
- Bilbao Crystallographic Server, programs *Point Group Tables* and
  *REPRESENTATIONS DPG*: M. I. Aroyo et al., *Bulg. Chem. Commun.* 43 (2011)
  183; L. Elcoro et al., *J. Appl. Cryst.* 50 (2017) 1457.
