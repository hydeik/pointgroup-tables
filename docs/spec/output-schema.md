# Output schema

Decided in wayfinder ticket
[Output schema: shared element and irrep identifiers across all tables](https://github.com/hydeik/pointgroup-tables/issues/11).
Builds on [Command line](./cli.md) (what one invocation produces and the
metadata content), [Values](./values.md) (the document field and the
`exact` / `numeric` / `pretty` forms), [Label conventions](./label-conventions.md)
(the structured element record), [Conjugacy classes](./conjugacy-classes.md)
(class records and the columns of the character table),
[Irrep labels](./irrep-labels.md) (the structured irrep record and the
canonical irrep order), [Double groups](./double-groups.md) (the element
order of `G†` and the factor table) and [Basis conventions](./basis-conventions.md)
(the matrices of the `matrices` table). Rationale in ADR-0005. The machine
form of this chapter is [`schema.json`](./schema.json), a JSON Schema
(draft 2020-12) that the `schema` subcommand prints. Version 2 adds, under
`--chain` only, the `chain` metadata key, the `subgroup_chain` table and
the `components` of the `matrices` entries
([Group chains](./group-chains.md#the-subgroup_chain-table-and-the-chain-labels),
decided in
[The `chain` table and the chain labels of the `matrices` components](https://github.com/hydeik/pointgroup-tables/issues/45));
a document produced without `--chain` is the version-1 document, byte for
byte. Version 2 also adds the `cg` document of the `cg` subcommand, which
holds the [`cg` table](#the-cg-table)
([Clebsch–Gordan coefficients](./clebsch-gordan.md#document), decided in
[The `cg` document and subcommand](https://github.com/hydeik/pointgroup-tables/issues/47));
no version-1 document changes.

Prototype: `prototypes/output_schema_prototype.py` produces every document
below for `C3v`, `D6h`, `Oh` and `Ih` (irrep matrices for the axial groups
only; the polyhedral matrices are the business of the canonical-basis
prototype). The `C3v` documents are checked in under `prototypes/samples/`.

## Principles

1. **One document, one model.** A document is a single JSON object; YAML and
   TOML are serialisations of the same object ([Values](./values.md#yaml-and-toml)).
   The model holds no `null`: an absent value is an absent key.
2. **Identifiers are indices.** Elements, classes and irreps are referred to
   everywhere by their position in the canonical element order, class order
   and irrep order. Labels and symbols are output only and never keys.
   Positional arrays, not objects keyed by number, carry the tables.
3. **Every document embeds its element registry.** The `elements` table is
   present in every document, so a single-table document can be read on its
   own, and `tables` always begins with `"elements"`. `all` remains the plain
   union of the five single-table documents; the `cg` document embeds it
   too. Under `--chain` the
   `subgroup_chain` table is present in every document in the same way,
   right after `elements`.
4. **A table is self-contained.** A table carries the records its own rows
   and columns are indexed by: the character table carries its class
   headers (which are the classes of `G†` in the default case, not present
   in `elements`), the matrices table carries its irrep records. The same
   irrep record appears in `characters` and `matrices` of an `all` document.
5. **Fixed key order.** Keys are emitted in the order this chapter lists
   them, so that documents are deterministic byte for byte.

## Document skeleton

```json
{
  "schema_version": 1,
  "generator": { "name": "pointgroup-tables", "version": "1.0.0" },
  "group": { "symbol": "C3v", "input_symbol": "c3v", "family": "Cnv", "n": 3,
             "order": 6, "type": "twisted", "rotation_group": "D3" },
  "tables": ["elements", "characters"],
  "double": false,
  "irreps": "all",
  "values": ["exact", "numeric"],
  "field": { "kind": "cyclotomic", "conductor": 12, "degree": 4,
             "generator": "exp(2*pi*i/12)", "basis": "power" },
  "label_convention": { "name": "altmann-herzig", "source": "...", "coverage": "reference" },
  "basis_convention": { "name": "canonical" },
  "chain": { "input": "C3", "canonical": "C3", "symbols": ["C3v", "C3"] },
  "elements":       { ... },
  "subgroup_chain": { ... },
  "cayley":         { ... },
  "factor":     { ... },
  "characters": { ... },
  "matrices":   { ... },
  "cg":         { ... }
}
```

Metadata keys are top-level and come first, in this order; the tables
follow, in the order `elements`, `subgroup_chain`, `cayley`, `factor`,
`characters`, `matrices`, `cg`, and only those named in `tables` (`cg`
only in a `cg` document, never with the other value tables). There is no
`metadata` wrapper: `tables` tells a reader which top-level keys are
tables. `chain` and `subgroup_chain` are present exactly when `--chain`
was given (the example above is `--chain C3` on `C3v`).

| Key | Content | Presence |
| --- | ------- | -------- |
| `schema_version` | `1` | always |
| `generator` | `{ name, version }` ([Command line](./cli.md#documents)) | always |
| `group` | `{ symbol, input_symbol, family, n, order, type, rotation_group }`; `n` absent for polyhedral groups; `family` is `Cn`, `Cnv`, `Cnh`, `S2n`, `Dn`, `Dnh`, `Dnd` or the polyhedral symbol | always |
| `tables` | names of the tables present, in the order above; first entry always `"elements"` | always |
| `double` | `true` when `--double-group` applied | always |
| `irreps` | the `--irreps` value | when `characters`, `matrices` or `cg` is present |
| `values` | the value forms present, in the order `exact`, `numeric`, `pretty` | always |
| `field` | the document field ([Values](./values.md#the-field-of-a-document)) | always |
| `label_convention` | `{ name, source, coverage }` | when `--label-convention` selected |
| `basis_convention` | `{ name, source, coverage }` for a reference convention; `{ name = "canonical" }` by default | always |
| `chain` | `{ input, canonical, symbols }` ([Group chains](./group-chains.md#metadata)); its presence declares that `matrices` is in the chain-adapted basis of that chain, `basis_convention` staying `{ name = "canonical" }` | with `--chain` |

**Coverage.** `coverage` is `"reference"` when every table of the document's
group that the convention touches (reference labels, or the basis) is
verified against the reference's printed tables, and `"extrapolated"` when
the convention's rules are applied beyond what the reference prints. The
convention file declares its range (`[coverage]` in
[Label conventions](./label-conventions.md#label-conventions-reference-labels)
and [Basis conventions](./basis-conventions.md#reference-basis-conventions):
`axial_max_n`, `polyhedral`), and the tool compares the group with it; for
`altmann-herzig` every axial group with `n ≤ 10` and every polyhedral group
is `"reference"`. The canonical labelling and the canonical basis have no
reference to be verified against, so the key is absent for them.

## Index spaces

- **Element index**: position in the canonical element order of `G`
  ([Abstract groups](./abstract-groups.md#canonical-element-order)),
  `0 ≤ i < |G|`. In a document with `double = true` the `elements` and
  `cayley` tables use the index space of `G†`, `0 ≤ i < 2|G|`, in which
  element `i + |G|` is the barred element covering `g_i`
  ([Double groups](./double-groups.md#double-group-tables)).
- **Barred indices in a `G` document.** The class records of `G†` inside a
  `characters` table use `G†` indices even when the document's `elements`
  table is that of `G` (the default: `--irreps all` or `projective` without
  `--double-group`). The rule is fixed: an index `i ≥ |G|` names the barred
  element of the element `i − |G|` of the `elements` table. No table ever
  refers to an element it cannot name this way.
- **Class index**: position in the canonical class order of `G`, or of
  `G†` (ordered by base class, [Conjugacy classes](./conjugacy-classes.md#classes-of-g)).
  Which one is named by the `group` / `class_group` key of the table that
  uses it.
- **Irrep index**: position in the canonical irrep order
  ([Irrep labels](./irrep-labels.md#canonical-irrep-order)) **of the full
  irrep list** of the group, ordinary then projective. Filtering by
  `--irreps` removes rows but does not renumber: `--irreps projective` of
  `C3v` lists the irreps with indices `3`, `4`, `5`. Rows of a table are
  therefore matched to irreps through the `index` field of the irrep
  records, and `partner` always refers to this global index.
- The `matrices` table is always indexed by the elements of `G`, also in a
  `double = true` document (the barred matrix is `−D(g)` for a projective
  irrep and `D(g)` for an ordinary one); map through `base_index` there.
- **Step and subgroup-irrep indices** (`--chain` only): a step is the
  position `j` of its subgroup record in `subgroup_chain.items`; an irrep
  of a step is its index in that record's `irreps`, the full irrep list of
  the step's abstract group in its own canonical irrep order. The chain
  labels of `matrices` components are sequences of these indices. The
  `embedding` of a subgroup record uses `G†` indices under the barred rule
  above; every other index of the table is a `G` element index.
- **Component and product indices** (`cg` only): a component index is a
  position `0 ≤ k < d_γ` in the component order of the basis the document
  names (under `--chain` the chain-adapted order, labelled by
  `cg.irreps[*].components`); a row of a coupling matrix is the product
  index `i d_β + j` of the components `i` of `α` and `j` of `β`; a column
  is a position in the pair record's `columns` array
  ([The `cg` table](#the-cg-table)).

## The `elements` table

```json
"elements": {
  "group": "G",
  "items":   [ <element record>, ... ],
  "classes": [ <class record>, ... ]
}
```

`group` is `"G"` or, with `--double-group`, `"G_double"`. `items` lists
every element in canonical order (the record at position `i` has
`index = i`); `classes` lists every conjugacy class of that group in
canonical class order.

**Element record** ([Label conventions](./label-conventions.md#structured-element-record)),
keys in this order:

| Key | Content | Presence |
| --- | ------- | -------- |
| `index` | element index | always |
| `label` | canonical label | always |
| `reference_label` | label of the selected convention | with a label convention |
| `parity` | `1` or `-1` | always |
| `kind` | `identity`, `rotation`, `inversion`, `reflection`, `rotoreflection` | always |
| `order` | `m` | always |
| `power` | `j` | order `≥ 3` |
| `axis_class` | `principal`, `secondary`, `secondary_double_prime`, `polyhedral` | when the label has an axis |
| `axis_index` | `k` | when the label has `(k)` |
| `inverse` | index of the inverse element | always |
| `class` | index of the class containing the element | always |
| `barred` | `true` / `false` | `G_double` only |
| `base_index` | index in `G` of the element covered | `G_double` only |

**Class record** ([Conjugacy classes](./conjugacy-classes.md#class-records)),
keys in this order: `index`, `label`, `reference_label` (with a label
convention), `size`, `representative`, `members`, `inverse`, then for `G`
`regular`, and for `G†` `base_class`, `split`, `barred` (split classes
only). `members` lists element indices with the representative first and
the rest ascending.

Example (`C3v`, class of `G†` containing `Ē·s(σ_v(1))`):

```json
{ "index": 5, "label": "3~sigma_v", "size": 3, "representative": 9,
  "members": [9, 4, 5], "inverse": 4, "base_class": 2, "split": true, "barred": true }
```

## The `subgroup_chain` table

Present exactly when `--chain` was given, listed in `tables` right after
`"elements"`, in every table subcommand; content fixed in
[Group chains](./group-chains.md#the-subgroup_chain-table-and-the-chain-labels).
It depends on the group and the resolved chain only (not on
`--double-group`, `--irreps` or `--values`), carries no value of the
document field, and is therefore byte-identical across every document of
one group under one chain.

```json
"subgroup_chain": {
  "irreps": [ <irrep record>, ... ],
  "items": [
    { "group": { "symbol": "C3", "family": "Cn", "n": 3, "order": 3, "type": "proper", "rotation_group": "C3" },
      "embedding": [0, 1, 2],
      "cosets": [[0, 1, 2], [3, 4, 5]],
      "double_cosets": [ { "representative": 0, "size": 3, "cosets": [0] },
                         { "representative": 3, "size": 3, "cosets": [1] } ],
      "irreps": [ <irrep record>, ... ],
      "branching": [[0], [0], [1, 2], [4, 5], [3], [3]] }
  ]
}
```

- `irreps`: the full irrep list of `G` (never filtered by `--irreps`), the
  same records as `characters.irreps` would hold unfiltered, with
  `reference_symbol` under a label convention.
- `items`: one subgroup record per step in chain order; empty for the
  empty chain. A subgroup record has no `index`.

**Subgroup record**, keys in this order:

| Key | Content | Presence |
| --- | ------- | -------- |
| `group` | `{ symbol, family, n, order, type, rotation_group }` of the step's abstract group (the `group` metadata record without `input_symbol`; `n` absent for a polyhedral step) | always |
| `embedding` | `G†` index of the image of the lifted preimage of every element of the step's abstract group, in its canonical element order (`order` entries) | always |
| `cosets` | the left cosets of the step in the previous group: arrays of `G` element indices, representative (smallest index) first, members ascending, cosets by representative ascending | always |
| `double_cosets` | records `{ representative, size, cosets }` (`cosets`: ascending indices into the `cosets` array), by representative ascending | always |
| `irreps` | the full irrep list of the step's abstract group, irrep records with canonical symbols only (no `reference_symbol`) | always |
| `branching` | one ascending array of indices into this record's `irreps` per irrep of the previous group (`subgroup_chain.irreps` for the first step, the previous record's `irreps` afterwards), in that group's irrep order | always |

## The `cayley` table

```json
"cayley": { "group": "G", "table": [[0, 1, 2, 3, 4, 5], [1, 2, 0, 5, 3, 4], ...] }
```

`table[i][j]` is the index of the product `g_i · g_j` (row: left factor,
column: right factor), a square array of `|G|` (or `2|G|` with
`--double-group`, `group = "G_double"`) rows of element indices.

## The `factor` table

```json
"factor": { "group": "G", "table": [[1, 1, 1, 1, 1, 1], [1, -1, 1, 1, 1, 1], ...] }
```

`table[i][j]` is `ω(g_i, g_j)` as the integer `1` or `-1`
([Double groups](./double-groups.md#factor-system)), a `|G| × |G|` array in
the canonical order of `G`. `group` is always `"G"`. The entries are plain
integers, not values of the document field.

## The `characters` table

```json
"characters": {
  "class_group": "G_double",
  "irreps":  [ <irrep record>, ... ],
  "classes": [ <class record>, ... ],
  "exact":   [ [ <exact value>, ... ], ... ],
  "numeric": [ [ <numeric value>, ... ], ... ],
  "pretty":  [ [ "...", ... ], ... ]
}
```

- `class_group` is `"G"` for `--irreps ordinary` and `"G_double"` otherwise
  ([Conjugacy classes](./conjugacy-classes.md#columns-of-the-character-table)).
- `irreps` holds one irrep record per row, in canonical irrep order,
  filtered by `--irreps`.
- `classes` holds the class records of `class_group`, one per column, in
  canonical class order: the full records, identical to those the
  `elements` table of that group would carry.
- `exact`, `numeric`, `pretty`: one two-dimensional array per requested
  form, `[row][column]` = `[irrep][class]`, all of the same shape
  ([Values](./values.md#layout-in-tables)).

**Irrep record** ([Irrep labels](./irrep-labels.md#structured-irrep-record)),
keys in this order: `index`, `symbol`, `reference_symbol` (with a label
convention), `bethe` (crystallographic groups), `dimension`, `projective`,
`letter`, `component` and `partner` (members of a complex-conjugate pair),
`subscript` (an integer, or a string such as `"3/2"`), `suffix`,
`frobenius_schur`.

## The `matrices` table

```json
"matrices": {
  "irreps": [
    { "index": 2, "symbol": "E", "dimension": 2, "projective": false, "letter": "E",
      "frobenius_schur": 1,
      "exact":   [ <exact matrix of element 0>, <exact matrix of element 1>, ... ],
      "numeric": [ <numeric matrix of element 0>, ... ] },
    ...
  ]
}
```

One entry per irrep, in canonical order, filtered by `--irreps`: the irrep
record followed by the requested value forms, each an array over the
elements of `G` in canonical order whose entries are `d × d` matrices as
nested rows ([Values](./values.md#layout-in-tables)). Matrices are in the
basis named by `basis_convention`, or, when the `chain` metadata key is
present, in the chain-adapted basis of that chain
([Group chains](./group-chains.md#chain-adapted-basis)); one-dimensional
irreps appear with `1 × 1` matrices. Values are per irrep, not table-wide,
because the dimension varies from irrep to irrep.

**Component records** (`--chain` only). Under `--chain` every entry carries,
between `frobenius_schur` and the value forms, a `components` array with
one record per component in the component order of the chain-adapted
basis, the chain labels of
[Group chains](./group-chains.md#chain-labels-of-the-matrices-components):

```json
{ "index": 2, "symbol": "E", "dimension": 2, "projective": false, "letter": "E",
  "frobenius_schur": 1,
  "components": [ { "index": 0, "labels": [1], "symbols": ["1E"] },
                  { "index": 1, "labels": [2], "symbols": ["2E"] } ],
  "exact": [ ... ], "numeric": [ ... ] }
```

| Key | Content |
| --- | ------- |
| `index` | the component's position, `0 ≤ index < dimension` |
| `labels` | one irrep index per step: `labels[j]` indexes `subgroup_chain.items[j].irreps` |
| `symbols` | the canonical Mulliken symbols of those irreps, display only |

The array has `dimension` entries; `labels` and `symbols` have one entry
per step (none under the empty chain). The irrep record proper is
unchanged, so it still equals the record in `characters`.

## The `cg` table

The one table of the `cg` document (version 2, decided in
[The `cg` document and subcommand](https://github.com/hydeik/pointgroup-tables/issues/47);
content fixed in [Clebsch–Gordan coefficients](./clebsch-gordan.md#document)).
It is never part of `all`, and a `cg` document embeds no `matrices` table:
the basis the coefficients are written in is named by `basis_convention`
and `chain`, exactly as for `matrices`.

```json
"cg": {
  "irreps": [ <irrep record>, ... ],
  "pairs": [
    { "irreps": [2, 2],
      "dimension": 4,
      "constituents": [ { "irrep": 0, "multiplicity": 1, "columns": [0, 1] },
                        { "irrep": 1, "multiplicity": 1, "columns": [1, 1] },
                        { "irrep": 2, "multiplicity": 1, "columns": [2, 2] } ],
      "columns": [ { "index": 0, "irrep": 0, "copy": 0, "component": 0, "exchange": 1 },
                   { "index": 1, "irrep": 1, "copy": 0, "component": 0, "exchange": -1 },
                   { "index": 2, "irrep": 2, "copy": 0, "component": 0, "exchange": 1 },
                   { "index": 3, "irrep": 2, "copy": 0, "component": 1, "exchange": 1 } ],
      "exact":   [ <4 rows of 4 exact values> ],
      "numeric": [ <4 rows of 4 numeric values> ] },
    ...
  ]
}
```

- `irreps`: the **full** irrep list of `G` in canonical irrep order, never
  filtered by `--irreps` (the same records as an unfiltered
  `characters.irreps`, `reference_symbol` under a label convention); under
  `--chain` every record carries, after `frobenius_schur`, the `components`
  array of the [`matrices` table](#the-matrices-table).
- `pairs`: one pair record per unordered pair `(α, β)`, `α ≤ β` in the
  global irrep index, both factors in the `--irreps` selection (`ordinary`:
  both ordinary; `projective`: at least one projective; `all`: every
  pair), in lexicographic order of `(α, β)`. The reversed pair is defined
  by the exchange relation of the chapter and is not printed.

**Pair record**, keys in this order:

| Key | Content |
| --- | ------- |
| `irreps` | `[α, β]`, global irrep indices, `α ≤ β` |
| `dimension` | `d_α d_β`, the order of the coupling matrix |
| `constituents` | constituent records, one per irrep `γ` with `n_γ > 0`, in canonical irrep order |
| `columns` | column records, `dimension` entries, by `γ` (canonical order), copy, component |
| `exact`, `numeric`, `pretty` | the coupling matrix, `dimension × dimension`, `[row][column]`, one per requested form ([Values](./values.md#layout-in-tables)) |

**Constituent record**, keys in this order: `irrep` (global index of `γ`),
`multiplicity` (`n_γ ≥ 1`), `columns` (`[first, count]`, `count = n_γ d_γ`).

**Column record**, keys in this order: `index` (position), `irrep` (`γ`),
`copy` (`0 ≤ copy < n_γ`), `component` (`0 ≤ component < d_γ`), `exchange`
(`1` or `-1`, present exactly when `α = β`).

Rows have no records: row `i d_β + j` is the product state of component `i`
of `α` and component `j` of `β` ([Index spaces](#index-spaces)).

## JSON formatting

The JSON serialisation is fixed so that documents are deterministic:

- UTF-8, ASCII only (`\uXXXX` escapes for anything else), LF line ends, a
  trailing newline;
- objects one key per line, two-space indentation, keys in the order of
  this chapter;
- an array whose entries are scalars or (recursively) arrays of scalars is
  written on one line, with `", "` and `": "` separators, when the line
  fits in 100 columns, and one entry per line otherwise;
- integers as decimal integers, floats as the shortest round-tripping
  representation ([Values](./values.md#numeric-form)).

YAML and TOML documents are dumps of the same object with the same key
order; their whitespace layout is that of the serialiser and is fixed per
tool version, not by this chapter. In TOML the metadata keys come first
(scalars and inline tables), followed by one `[elements]`, `[cayley]`, ...
table each; record lists are arrays of inline tables.

## `schema.json`

[`schema.json`](./schema.json) (`$schema` draft 2020-12, `$id`
`urn:pointgroup-tables:schema:1`) describes every table document:

- every metadata key and every table object, with `additionalProperties`
  closed; each table named in `tables` must be present, `irreps` must
  accompany `characters`, `matrices` or `cg`, and `chain` and `subgroup_chain`
  are present together or not at all (`components` is allowed in a
  `matrices` entry or a `cg.irreps` record only when `chain` is present);
- `$defs` for the element, class, irrep, subgroup, double-coset,
  component, pair, constituent and column records, the exact, numeric and pretty values and matrices, and
  the label grammars (element label, class label, irrep symbol) as regular
  expressions taken from their chapters;
- what a schema cannot express (indices in range, `members` partitioning
  the elements, array shapes matching `dimension` and `degree`, equality of
  the irrep records across tables) belongs to the verification chapter.

The `list groups` and `list conventions` documents are not covered by
`schema.json`; their shape is fixed by [Command line](./cli.md#utility-subcommands).

## Sizes

Measured on the prototype, JSON, `--values both` (exact and numeric):

| Document | Size |
| -------- | ---: |
| `all C3v` | 15 KB |
| `elements D6h --double-group` | 18 KB |
| `characters D6h` | 38 KB |
| `matrices D6h --values exact` | 68 KB |
| `all D6h` | 137 KB |
| `all Oh` without `matrices` | 94 KB |
| `elements Ih` | 34 KB |
| `all Ih` without `matrices` | 419 KB |

The `subgroup_chain` table is small next to these (indices and irrep
records only: a few KB for `Ih --chain default`, whose largest record holds
the twenty irrep records of `S10`); `components` adds a few bytes per
matrix row.

`matrices Ih` (240 matrix entries per element, 120 elements, conductor 20,
both forms) is estimated at 5–6 MB. The `cg` documents are estimated
([Clebsch–Gordan coefficients](./clebsch-gordan.md#sizes)) at 20 KB for
`C3v`, 150 KB for `D6h`, 0.5 MB for `Oh` and 15–20 MB for `Ih`, the
largest document of the tool (31428 coefficients of degree 32, each
exact value broken across lines by the 100-column rule). Whether a cap or a streaming layout is
ever needed for large `n` is an open point of the map, not of this chapter.

## Verification

An implementation must check:

- every document of every table subcommand, for every group in the seed
  list and every combination of `--double-group`, `--irreps` and
  `--values`, validates against `schema.json`;
- `tables` lists exactly the table keys present, beginning with
  `elements`; `values` lists exactly the forms present, and all forms of a
  table have the same shape;
- every index is in range for its index space (the barred rule for `G†`
  classes in a `G` document included); `items[i].index = i`,
  `classes[c].index = c`, and `matrices` entries and `characters` rows
  follow the canonical irrep order;
- the irrep records of `characters` and `matrices` in an `all` document are
  identical, and `characters.classes` equals `elements.classes` of the
  corresponding `elements --double-group` (or plain) document;
- `all` equals the union of the five single-table documents under the same
  options, table by table, and the `elements` table of every single-table
  document equals the `elements` document's;
- two runs with the same arguments produce byte-identical output in every
  format ([Command line](./cli.md#verification));
- under `--chain`: `chain` and `subgroup_chain` are both present,
  `subgroup_chain` is the second entry of `tables`, every `matrices` entry
  has `dimension` component records whose `labels` index the step records'
  `irreps` in range, and the structural checks of
  [Group chains](./group-chains.md#verification-of-the-table) hold; without `--chain`
  none of the three keys is present;
- a `cg` document: `tables` is `["elements", "cg"]` or
  `["elements", "subgroup_chain", "cg"]`, `irreps` is present, `double` is
  `false`, and the structural checks of
  [Clebsch–Gordan coefficients](./clebsch-gordan.md#document-checks) hold
  (full `irreps` list, pair set and order, `dimension`, `constituents`,
  `columns`, the `exchange` rule).
