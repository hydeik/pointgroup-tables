# Reference data

Transcriptions of reference works used by the `reference`-marked tests
(`docs/spec/verification.md`). They are test data: nothing under this
directory is shipped in the package. Layout is
`<source>/<kind>/<group>.toml`; each kind has a `README.md` describing its
body, and every file carries the `group` / `table` / `status` / `notes`
header.

## Sources

- **Altmann & Herzig**: S. L. Altmann and P. Herzig, *Point-Group Theory
  Tables*, 2nd ed. (Wien: Clarendon/Oxford; reprinted). Element labels,
  multiplication tables with tilde marks, character tables and irrep
  matrices, transcribed by hand, table and page cited per file.
- **Bilbao**: Bilbao Crystallographic Server, point-group and
  double-point-group tables (M. I. Aroyo et al., *Z. Kristallogr.* 221
  (2006) 15–27; *Acta Cryst.* A62 (2006) 115–128). Character tables of the
  32 crystallographic point groups, fetched once by `bilbao/fetch.py`, page
  URL cited per file. The server asks that its use be cited; users of this
  repository who republish the data should keep this citation.
- **Koster**: G. F. Koster, J. O. Dimmock, R. G. Wheeler and H. Statz,
  *Properties of the Thirty-Two Point Groups* (MIT Press, 1963). Bethe
  symbols `Γ_i` with identifying characters.

## Licensing

The files re-express mathematical facts (characters, products, names of
symmetry operations) in this project's own format, limited to what the
tests need. They are not reproductions of the printed pages, tables or
their arrangement. They are distributed under the repository licence; the
decision and its reasoning are recorded in
`docs/adr/0006-reference-data-redistribution.md`.
