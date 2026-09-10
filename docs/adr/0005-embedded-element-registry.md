---
status: accepted
---

# Every document embeds its element registry and refers to elements, classes and irreps by index

Tables refer to elements, conjugacy classes and irreps only by their
position in the canonical orders; labels and symbols are output, never
keys. To keep a single-table document readable on its own, every document
carries the `elements` table (element and class records of `G`, or of `G†`
with `--double-group`), and every table carries the records its own rows and
columns are indexed by: the character table its class headers, the matrices
table its irrep records. `all` stays the plain union of the five documents.
See `docs/spec/output-schema.md`.

## Considered Options

- **Bare indices, no embedded registry.** Smallest documents, but a
  `cayley` or `factor` table is then a square of integers that cannot be
  read without a second invocation, and the character table of `G†` has no
  column headers anywhere in the document (the default `characters`
  document indexes `G†` classes while `elements` holds those of `G`).
  Rejected for usability; the registry costs a few percent of a matrices
  document and 34 KB for `Ih`.
- **Labels as keys** (`"C3^1": ...`). Rejected: labels depend on the
  selected convention and would make the JSON, YAML and TOML forms diverge
  (object keys are strings), and the canonical orders already give every
  object a stable integer identity.
- **A separate `classes` table or a document-level irrep registry.**
  Rejected: the conjugacy-class chapter places classes in the `elements`
  table, and `irreps` is already the name of a metadata key. Duplicating a
  small irrep record in `characters` and `matrices` is cheaper than a third
  place to look.
- **A `metadata` wrapper object.** Rejected: `tables` already discriminates
  tables from metadata, and a flat top level maps directly onto TOML.

## Consequences

- A document of any table subcommand is self-describing: every index it
  contains resolves inside it, with one fixed rule for `G†` indices in a
  `G` document (`i ≥ |G|` is the barred element of `i − |G|`).
- Irrep indices are global to the group's full irrep list and survive
  `--irreps` filtering, so `partner` and cross-table references stay valid.
- Any change to the element, class or irrep records changes every document;
  `schema_version` guards consumers.
