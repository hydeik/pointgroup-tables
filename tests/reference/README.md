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
  multiplication tables with tilde marks, character tables, irrep
  matrices and symmetrized bases, transcribed by hand, table and page
  cited per file.
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
- **Chen–Fan–Ping** (`chen-fan-ping/`, version 2): the papers of J.-Q.
  Chen, P.-D. Fan, J.-L. Ping and co-workers on the eigenfunction method
  (J. Math. Phys. 39 (1998) 467, 5502, 5519; 41 (2000) 8196; Int. J.
  Quantum Chem. 74 (1999) 7; 75 (1999) 67; 76 (2000) 585; 83 (2001) 259;
  91 (2003) 541; Acta Cryst. A55 (1999) 871, 1049; Comput. Phys. Commun.
  120 (1999) 71), cited by the keys of `docs/research/eigenfunction-method.md`.
  Multiplication tables with tildes and element labels of `T†`, `O†`,
  `D6†`, Euler angles of `I†`, chain-adapted irrep matrices and
  Clebsch–Gordan tables, transcribed from the page images, paper, table
  and page cited per file. Layout in
  `docs/spec/verification.md` ("The Chen–Fan–Ping transcriptions").

## Licensing

The files re-express mathematical facts (characters, products, names of
symmetry operations) in this project's own format, limited to what the
tests need. They are not reproductions of the printed pages, tables or
their arrangement. They are distributed under the repository licence; the
decision and its reasoning are recorded in
`docs/adr/0006-reference-data-redistribution.md`.
