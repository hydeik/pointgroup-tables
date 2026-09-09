# Command line

Decided in wayfinder ticket
[CLI specification: subcommands, flags, group names and errors](https://github.com/hydeik/pointgroup-tables/issues/8).
Builds on [Abstract groups](./abstract-groups.md) (symbols, aliases, exposed
structure), [Double groups](./double-groups.md) (`--double-group`),
[Label conventions](./label-conventions.md) (`--label-convention`),
[Irrep labels](./irrep-labels.md) (ordinary versus projective irreps) and
[Values](./values.md) (`--values`, the document field). The document layout
itself is the business of the output-schema chapter; this chapter fixes what
the command line accepts, what one invocation produces, and how it fails.

## Shape

```
pointgroup-tables <subcommand> [<group>] [options]
pointgroup-tables --version | --help
```

Every option is written after the subcommand. There are two kinds of
subcommand:

| Kind | Subcommands | Positional argument |
| ---- | ----------- | ------------------- |
| table | `elements`, `cayley`, `factor`, `characters`, `matrices`, `all` | exactly one group |
| utility | `list groups`, `list conventions`, `schema` | none |

One invocation of a table subcommand produces exactly one **document** for
exactly one group. Several groups are handled by several invocations; the
tool never emits an array of documents.

## The group argument

The group is named by its Schoenflies symbol and nothing else. Matching is
**case-insensitive**: `C3v`, `c3v` and `C3V` name the same group. The
argument must match, case-insensitively,

```
^([A-Za-z]+)(0|[1-9][0-9]*)?([A-Za-z]*)$
```

`n` is a decimal integer without sign and without leading zeros (`C03` does
not match the grammar). The letters before and after `n` form the **family
key**, which must be one of:

| Family key | `n` | Names |
| ---------- | --- | ----- |
| `C`, `Cv`, `Ch`, `S`, `D`, `Dh`, `Dd` | required, `n ≥ 1`; `S` requires even `n` | the axial families `Cn`, `Cnv`, `Cnh`, `S2n`, `Dn`, `Dnh`, `Dnd` |
| `T`, `Td`, `Th`, `O`, `Oh`, `I`, `Ih` | forbidden | the polyhedral groups |
| `Cs`, `Ci`, `V`, `Vh`, `Vd` | forbidden | aliases, normalised as in [Abstract groups](./abstract-groups.md#aliases) |

`n` has no upper bound in this specification: a large `n` is slow, not
invalid. (Whether an implementation limit is ever needed is an open point of
the map, not of this chapter.)

The output records the argument as typed in `input_symbol` and the
normalised symbol (canonical spelling: family letter upper case, suffix
lower case, aliases and `D1`, `D1h`, `D1d`, `C1v` resolved) in `symbol`, as
fixed by [Abstract groups](./abstract-groups.md#exposed-structure).

## Table subcommands

| Subcommand | Tables in the document | Notes |
| ---------- | ---------------------- | ----- |
| `elements` | element list | of `G`, or of `G†` with `--double-group` |
| `cayley` | Cayley table | of `G`, or of `G†` with `--double-group` |
| `factor` | factor table `ω(g, h)` | always on `G`; no `--double-group` |
| `characters` | character table | ordinary and projective irreps, filtered by `--irreps` |
| `matrices` | irrep matrix table | in the selected basis convention, filtered by `--irreps` |
| `all` | the union of the five | `--double-group` acts on the element list and Cayley table only |

`all` is exactly the union of the other five documents' tables under the
same options, with one shared metadata block; it adds nothing and reorders
nothing. Whether a single-table document also embeds the element list its
indices refer to is decided by the output-schema chapter.

## Options

Common to every table subcommand (accepted everywhere, always recorded in the
metadata, even where they do not change the table):

| Option | Values | Default |
| ------ | ------ | ------- |
| `-f`, `--format` | `json`, `yaml`, `toml` | `json` |
| `-o`, `--output <path>` | a file path; `-` is standard output | standard output |
| `--values <set>` | comma-separated subset of `exact`, `numeric`, `both`, `pretty` ([Values](./values.md)) | `both` |
| `--label-convention <name>` | a shipped label convention | none (canonical labels only) |
| `--basis-convention <name>` | `canonical` or a shipped basis convention | `canonical` |

The output format is never inferred from the file name: `-o table.yaml`
without `--format yaml` writes JSON. An existing output file is overwritten.

Specific to some subcommands (a usage error elsewhere):

| Option | Subcommands | Values | Default |
| ------ | ----------- | ------ | ------- |
| `--double-group` | `elements`, `cayley`, `all` | flag | off |
| `--irreps` | `characters`, `matrices`, `all` | `ordinary`, `projective`, `all` | `all` |

`--irreps` is the only irrep filter of version 1. Irrep names are output
only ([Irrep labels](./irrep-labels.md)), so no option selects an irrep by
symbol; consumers filter the document.

Global: `--version` prints `pointgroup-tables <version>` on one line, the
version taken from the installed package metadata; `--help` prints usage.
There are no verbosity, quiet, colour or configuration-file options.

## Documents

A document is self-describing and **deterministic**: the same arguments and
the same tool version produce byte-identical output. No timestamp and no
echo of the command line is ever embedded. The metadata of every table
document carries:

| Field | Content |
| ----- | ------- |
| `schema_version` | integer, `1` for this specification |
| `generator` | `{ name = "pointgroup-tables", version = <tool version> }` |
| `group` | `symbol`, `input_symbol`, `family`, `n` (absent for polyhedral groups), `order`, `type`, `rotation_group` ([Abstract groups](./abstract-groups.md#exposed-structure)) |
| `tables` | the names of the tables the document holds, in the order above (`["elements"]`, `["elements", "cayley", "factor", "characters", "matrices"]`, ...) |
| `double` | `true` when `--double-group` applied, else `false` |
| `irreps` | the `--irreps` value; present only when the document holds a character or matrix table |
| `values` | the set of value forms present, as requested by `--values` |
| `field` | the document field ([Values](./values.md#the-field-of-a-document)), fixed by group and basis convention |
| `label_convention` | `{ name, source }`; present only when one is selected |
| `basis_convention` | `{ name, source }`; `{ name = "canonical" }` by default |

The exact placement of these keys, and the layout of every table, belong to
the output-schema chapter.

## Utility subcommands

`list groups` describes what the group argument accepts, in `--format`:
three arrays, `families` (one entry per axial family: family key, canonical
family name, the constraint on `n`, the group type rule and rotation group
from the [family table](./abstract-groups.md#family-table), and one
example symbol), `polyhedral` (the seven polyhedral groups with type and
rotation group) and `aliases` (input symbol and normalised symbol). It does
not enumerate concrete axial groups.

`list conventions` lists what `--label-convention` and `--basis-convention`
accept, in `--format`: `label_conventions` and `basis_conventions`, each an
array of `{ name, source }`. `canonical` appears among the basis
conventions; the canonical labelling is not a label convention and is not
listed.

`schema` writes the JSON Schema of the table documents (the `schema.json`
shipped with the output-schema chapter) to standard output or `-o`. It is
always JSON; `--format` is not accepted.

Both `list` subcommands accept `-f` and `-o`; nothing else.

## Exit codes and errors

| Code | Meaning | Examples |
| ---- | ------- | -------- |
| `0` | success | |
| `1` | any other failure | output file not writable, internal error, a self-check failing |
| `2` | usage error | unknown subcommand, unknown option, an option on a subcommand that does not take it, a group argument that does not match the grammar (`3v`, `C03`, `C3-v`), a missing group |
| `3` | unknown group | the argument matches the grammar but its family key names no group: `Xh`, `Cnd`, `Dnv`, `T5`, `Cs3`, a bare `C` |
| `4` | unknown convention | a `--label-convention` or `--basis-convention` name that is not shipped; also any error the basis-convention chapter defines for a convention that does not cover the requested group |
| `5` | unsupported `n` | the family key is known and takes `n`, but `n` is out of range: `C0`, `D0h`, `S3` (`S2n` takes even subscripts only) |

Errors are written to standard error, never to the document stream, and
nothing is written to `-o` on failure. The first line is

```
pointgroup-tables: error: <message>
```

optionally followed by one line `hint: <hint>`. The message names the
offending value in single quotes. Hints point at the utility subcommand that
lists the accepted values:

```
pointgroup-tables: error: unknown group 'S3': S2n takes even subscripts only
hint: run 'pointgroup-tables list groups' for the accepted symbols

pointgroup-tables: error: unknown label convention 'bilbao'
hint: run 'pointgroup-tables list conventions' for the shipped conventions
```

Usage errors (code `2`) use the argument parser's usual usage line followed
by the same `pointgroup-tables: error: ...` line.

## Examples

```
pointgroup-tables characters C3v
pointgroup-tables characters c3v --format yaml --values exact
pointgroup-tables matrices D6h --basis-convention altmann-herzig -o d6h.json
pointgroup-tables elements Oh --double-group --label-convention altmann-herzig
pointgroup-tables all Ih --irreps projective --values numeric,pretty
pointgroup-tables list groups
pointgroup-tables schema -o schema.json
```

## Verification

An implementation must pass:

- a table of group arguments against `(exit code, symbol)`: at least
  `C3v → (0, C3v)`, `c3V → (0, C3v)`, `Cs → (0, C1h)`, `vd → (0, D2d)`,
  `D1h → (0, C2v)`, `S2 → (0, S2)`, `Ih → (0, Ih)`, `3v → 2`, `C03 → 2`,
  `Cnd → 3`, `T5 → 3`, `C → 3`, `Xh → 3`, `C0 → 5`, `S3 → 5`;
- an unknown convention name on each convention option gives `4`, a
  subcommand-specific option on the wrong subcommand gives `2`, and on
  every failure the document stream is empty and no output file is created;
- determinism: two invocations with the same arguments produce
  byte-identical output, for every subcommand and format;
- `all` equals the union of the five single-table documents under the same
  options, table by table.
