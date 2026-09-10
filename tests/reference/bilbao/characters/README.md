# Bilbao character tables

Character tables of the 32 crystallographic point groups and of their double
groups, fetched from the Bilbao Crystallographic Server for wayfinder ticket
[Fetch Bilbao character tables for the 32 crystallographic groups](https://github.com/hydeik/pointgroup-tables/issues/20).
Two TOML files per group, named by the normalised Schoenflies symbol (`Cs` is
`C1h`, `Ci` is `S2`, `C3i` is `S6`):

- `<group>.toml` (`double = false`): the **Point Group Tables** page (program
  POINT, `https://www.cryst.ehu.eus/cgi-bin/rep/programs/sam/point.py?sg=<sg>&num=<n>`):
  the classes of `G` under the server's Seitz-like headers and the ordinary
  irreps with their Mulliken and `Gamma` labels;
- `<group>-double.toml` (`double = true`): the **REPRESENTATIONS DPG** page
  (`https://www.cryst.ehu.eus/cgi-bin/cryst/programs/representations_out.pl?tipogrupo=dbg&pointspace=point&num=<sg>&super=<n>`):
  the classes of `G†` listed operation by operation and every irrep, ordinary
  and double-valued, with three labels and the server's reality.

Produced and checked by

```
uv run python tests/reference/bilbao/fetch.py --raw-dir RAW --turnstile-cookie VALUE   # fetch (once) and write the files
uv run python tools/bilbao_characters.py check                                          # structural comparison with the computed tables (REF-2)
```

The tests read the TOML only; the script is the record of provenance and is
not run by the tests. See `tests/reference/README.md` for the server's
citation request and licence.

## File layout

```toml
group = "D3h"
table = "https://www.cryst.ehu.eus/..."      # the page the file was written from
status = "done"
source = "Bilbao Crystallographic Server, ..."   # program and citation
fetched = "2026-09-10"
double = true                       # false: POINT page (G); true: REPRESENTATIONS DPG page (G†)
hermann_mauguin = "-6m2"
notes = ""
conductor = 12                      # field of `characters`; zeta = exp(2 pi i / conductor)
server_w = "exp(2 pi i / 3)"        # POINT pages only, when the page defines w
server_notes = ["..."]              # POINT pages only: notes printed under the table

[[classes]]                         # the server's columns, in the server's order
label = "C5"                        # DPG: C1 .. Ck; POINT: the header (1, 3+, 2_100, -4, m_1-10, ...)
elements = ["m_110", "m_100", "m_010", "dm_110", "dm_100", "dm_010"]   # DPG only: the operations, Seitz symbols in ASCII
size = 6
order = 4                           # order of the elements (in G† for the double file)
parity = -1                         # +1 proper, -1 improper
barred = "mixed"                    # DPG only: "none", "all" or "mixed" d-operations in the class

[[irreps]]                          # the server's rows, in the server's order
label = "2~E1"                      # Mulliken column, ASCII
bethe = "Gamma7"                    # Bradley & Cracknell Gamma column (Koster numbering), "Gamma5+" in centrosymmetric groups
cdml = "~Gamma7"                    # DPG only: the CDML column
projective = true                   # double-valued (the character on the class of d1 is minus the dimension)
reality = -1                        # DPG only: the server's reality, 1 real, 0 complex, -1 pseudoreal
server_values = ["2", "-1", ...]    # the printed expressions (i, sqrt2, sqrt3, w)
characters = ["2", "-1", ...]       # the same values in the value grammar for `conductor`
```

ASCII spellings. A Seitz symbol keeps its sense and direction (`3+_001`,
`-4-_001`, `m_1-10`), the `d` superscript of a spin-rotated operation becomes
a `d` prefix (`d2_010`, `d1` for the full turn), and overlined minus signs
become `-`. An irrep label keeps the superscript prefix and the subscript
(`1E2g`, `A1''`), and the overline of a double-valued label becomes `~`
before the barred token, as the `~` of a barred element (`1~E1g`, `~F`,
`~Gamma7`).

Values. The server prints `w` (= `exp(2 pi i / 3)`, stated in `server_w`)
and `1j` (= `i`) on the POINT pages and closed expressions in `i`, `sqrt2`,
`sqrt3` on the DPG pages; `fetch.py` evaluates them exactly in
`Q(zeta_conductor)` and writes `characters` with `to_grammar` of the basis
prototype, checking every value numerically against the printed expression.
Element orders in the double file come from the server's own 3 x 3 and SU(2)
matrices (the order of the pair) and are checked against the Seitz symbol.

## What the checker found

Every one of the 64 files matches the computed table structurally: a
bijection of classes (by size, element order and parity) and of irreps
reproduces the server's characters exactly, and the server's reality equals
the Frobenius-Schur indicator. The checker also reports how the server's
labels sit on the canonical ones, with the server's `n+` about `001` taken
as the canonical `Cn` and `-n+` as `i Cn`; the outcome is recorded in the
resolution of the ticket and summarised here.

- **Exponent sign of complex pairs.** The server's `1E` carries the
  *negative* exponent on the principal rotation (`exp(-2 pi i / n)` on
  `n+`), the opposite of the canonical rule, in every cyclic-type group
  (`C3`, `S6`, `C4`, `S4`, `C4h`, `C3h`, `C6`, `C6h`) and for the
  double-valued pairs of `E3/2` type (`C3h`, `C4`, `S4`, `C4h`, `C6`, `C6h`,
  `T`, `Th`) and of `E1/2` type in `C3` and `S6`; it agrees with the
  canonical rule for the `E2`-type pairs of `C6` and `C6h`, for the `E`
  pairs of `T` and `Th`, and for the double-valued pairs of `C1h`, `C2`,
  `C2h`, `C3v`, `D3`, `D3d` and the `E1/2` type of `C4`, `S4`, `C4h`, `C3h`,
  `C6`, `C6h`.
- **The two pages disagree** on the ordinary irreps of `C6` and `C6h`: the
  POINT page's `1E2`/`2E2` (`Gamma3`/`Gamma2`) are the DPG page's
  `1E1`/`2E1` (`Gamma6`/`Gamma5`) and vice versa (the POINT page swaps the
  `E1`/`E2` subscripts relative to Mulliken's rule, which the DPG page
  follows, and prints a note about differences between Bradley & Cracknell
  and Altmann & Herzig); on `D3h`: the POINT page's `A2'` (`Gamma2`) is the
  DPG page's `A1''` (`Gamma3`) and vice versa (the DPG page's `A1''` is
  symmetric under `sigma_h`); and on the spelling of the trivial irrep of
  `C3` (`A` on POINT, `A1` on DPG). Every other group carries the same
  Mulliken and `Gamma` labels on both pages.
- **Spinor section.** The server's SU(2) matrices follow the conjugate sign
  convention (`3+_001` is `diag(exp(i pi/3), exp(-i pi/3))`, the canonical
  lift of `C3` is `diag(exp(-i pi/3), exp(i pi/3))`), so in the server's
  double group `(6+)^3` is the barred `d2_001` while canonically `C6^3 = C2`.
  The structural comparison therefore pairs the server's lifted `6+`, `6-`,
  `-6+`, `-6-` of `C3h`, `C6` and `C6h` with canonical barred elements; the
  match is exact under that pairing.
- **Class counts** agree with the computed tables for every group; the
  server merges a `d` operation into the class of its partner where the two
  are conjugate (`C2` axes perpendicular to a `C4` or `C6` axis and the
  corresponding planes), as the canonical tables do.
