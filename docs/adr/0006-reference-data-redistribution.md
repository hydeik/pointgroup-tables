---
status: accepted
---

# Transcribed reference tables are redistributed in the repository as test data

The verification chapter compares the tool's tables against reference works
(Altmann & Herzig, the Bilbao Crystallographic Server, Koster et al.).
The transcriptions of those references, character tables, element labels,
tilde-marked products and generator matrices, are committed under
`tests/reference/` and redistributed with the repository under its licence.
Each file is written in the project's own TOML format, restricted to what
the tests need, and cites its source (book, edition, table, page; or the
server page). `tests/reference/README.md` carries the citations, the Bilbao
terms of use, and the statement that the files are not reproductions of the
printed pages. See `docs/spec/verification.md`.

## Considered Options

- **Redistribute the transcriptions** (chosen). Character tables,
  multiplication tables and symmetry-operation names are mathematical
  facts, and facts are not protected; the risk lies only in copying the
  selection and arrangement of a whole book, which the project's own
  format and test-driven selection avoid. Open-source projects (GAP's
  character table library, spglib, phonopy) ship point-group character
  tables on the same basis. With the data in the repository the reference
  comparison runs in continuous integration for everyone, and the
  transcription is reviewed like code.
- **Keep the Altmann & Herzig transcriptions out of version control**
  (`.gitignore`, every developer transcribes locally; Bilbao data only).
  Rejected: the reference comparison would run on one machine, the
  transcription would be unreviewed and unreproducible, and the proof of
  convention would be optional in practice.
- **Keep every transcription out of version control.** Rejected for the
  same reasons, with no reference comparison left at all.

## Consequences

- The proof of convention is part of the shipped test suite; a
  contributor's change that breaks agreement with the literature fails
  continuous integration.
- Every reference file must carry its source header; a file without a
  citation is rejected in review.
- If a rights holder objects, the affected files are removed and the
  tests that need them fall back to the `reference` marker's skip. The
  decision to publish is hard to reverse in the history, which is why it
  is recorded here.
- Version 2 (ticket #48, 2026-10-28) extends the decision unchanged to the
  journal papers of Chen, Fan, Ping and co-workers
  (`tests/reference/chen-fan-ping/`): multiplication tables, irrep
  matrices and Clebsch–Gordan coefficients are mathematical facts, the
  files hold the entries the tests need in the project's format, and each
  cites paper, table and page.
