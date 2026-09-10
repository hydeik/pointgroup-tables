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
- **Bilbao**: Bilbao Crystallographic Server (`https://www.cryst.ehu.eus`),
  programs *Point Group Tables* (POINT) and *REPRESENTATIONS DPG*. Character
  tables of the 32 crystallographic point groups and of their double groups,
  fetched once by `bilbao/fetch.py` (September 2026), page URL cited per
  file. The server's pages are published under CC BY-NC-SA 4.0 and the
  server asks to be cited: generically M. I. Aroyo et al., *Bulg. Chem.
  Commun.* 43 (2011) 183–197; *Z. Kristallogr.* 221 (2006) 15–27; *Acta
  Cryst.* A62 (2006) 115–128; for POINT the *Acta Cryst.* A62 paper; for
  REPRESENTATIONS DPG L. Elcoro et al., *J. Appl. Cryst.* 50 (2017)
  1457–1477. Users of this repository who republish the data should keep
  these citations.
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
