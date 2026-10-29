# Chain prototype data: summary

Generated from `chain_prototype_conductors.json` and `chain_prototype_norm_factors.json` by
`uv run python3 prototypes/chain_cg_prototype.py chains` (or `summary`);
never edited by hand. `docs/spec/group-chains.md` points here for every
count.

## Chain classes and conductors per group

`Spellings` counts the entries of `list chains` (one per canonical
spelling, the reached chain); `Classes` the accepted chain classes, which
exceed the spellings by the unreached twin chain classes (ticket #59).

| Group | Spellings | Classes | Default | `M_chain` × spellings |
| ----- | --------: | ------: | ------- | --------------------- |
| `T` | 2 | 2 | `C3` | 12 × 1, 24 × 1 |
| `Td` | 11 | 12 | `S4` | 8 × 7, 24 × 4 |
| `Th` | 13 | 15 | `S6` | 12 × 3, 24 × 10 |
| `O` | 11 | 12 | `C4` | 8 × 7, 24 × 4 |
| `Oh` | 143 | 160 | `C4h` | 8 × 107, 24 × 36 |
| `I` | 4 | 4 | `D5>C5` | 20 × 2, 60 × 1, 120 × 1 |
| `D6h` | 108 | 108 | `C6h` | 12 × 108 |
| `Ih` | 32 | 34 | `D5d>S10` | 20 × 15, 60 × 5, 120 × 12 |

324 spellings (347 chain classes) in 8 groups.

## Twin chain classes

Accepted chain classes that share one canonical spelling (ticket #59): the
embedding rule reaches the first, `--chain` cannot name the other. `Step`
is the step at which the unreached class differs, with the subgroup each
class holds there; `same lines` counts the irreps whose label spaces
coincide under the two chains.

| Group | Spelling | Step | Reached | Unreached | Same lines |
| ----- | -------- | ---- | ------- | --------- | ---------- |
| `Td` | `D2d>D2>C2` | 3 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 5/5, projective 0/3 |
| `Th` | `C2v>C1h` | 2 `C1h` | `{sigma_h(2)}` | `{sigma_h(3)}` | ordinary 8/8, projective 0/6 |
| `Th` | `D2h>C2v>C1h` | 3 `C1h` | `{sigma_h(2)}` | `{sigma_h(3)}` | ordinary 8/8, projective 0/6 |
| `O` | `D4>D2>C2` | 3 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 5/5, projective 0/3 |
| `Oh` | `D2d[sigma_d]>D2>C2` | 3 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4>D2>C2` | 3 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>C2v[C2,sigma_h]>C1h` | 2 `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}` | `{C2(2), sigma_h(1), sigma_h(3)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>C2v[C2,sigma_h]>C1h` | 2 `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}` | `{C2(2), sigma_h(1), sigma_h(3)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>C2v[C2,sigma_h]>C2` | 2 `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}` | `{C2(2), sigma_h(1), sigma_h(3)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2>C2` | 3 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2h>C1h` | 3 `C1h` | `{sigma_h(1)}` | `{sigma_h(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2h>C2` | 3 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2h>C2h` | 3 `C2h` | `{C2(1), i, sigma_h(1)}` | `{C2(2), i, sigma_h(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2d[sigma_d]>D2>C2` | 4 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2h>C2v>C1h` | 3 `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}` | `{C2(2), sigma_h(1), sigma_h(3)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2h>C2v>C1h` | 3 `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}` | `{C2(2), sigma_h(1), sigma_h(3)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2h>C2v>C2` | 3 `C2v` | `{C2(1), sigma_h(2), sigma_h(3)}` | `{C2(2), sigma_h(1), sigma_h(3)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D2h>D2>C2` | 4 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `D4h>D4>D2>C2` | 4 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `O>D4>D2>C2` | 4 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Oh` | `Td>D2d>D2>C2` | 4 `C2` | `{C2(1)}` | `{C2(2)}` | ordinary 10/10, projective 0/6 |
| `Ih` | `Th>C2v>C1h` | 3 `C1h` | `{sigma(11)}` | `{sigma(14)}` | ordinary 10/10, projective 0/8 |
| `Ih` | `Th>D2h>C2v>C1h` | 4 `C1h` | `{sigma(11)}` | `{sigma(14)}` | ordinary 10/10, projective 0/8 |

## Norm-factor entries

37 entries `(M, N) ↦ μ`, 8 of them printed in
the chapter; `uses` counts the chain classes and
irreps (sub-chains included) whose designated squared norm is `N`.

| `M` | `N` | printed | uses | first use |
| --: | --- | ------- | ---: | --------- |
| 24 | `1/2 + 1/4*sqrt2` | no | 36 | O D3>C2 E1/2 |
| 24 | `1/2 + 1/6*sqrt3` | yes | 133 | I T>D2>C2 (sub-chain T D2>C2) 2E3/2 |
| 24 | `1/2 - 1/4*sqrt2` | no | 8 | Oh C3v>C1h E1/2u |
| 24 | `1/2 - 1/6*sqrt3` | yes | 155 | I T>D2>C2 (sub-chain T D2>C2) 1E3/2 |
| 24 | `1/4 + 1/12*sqrt3` | yes | 16 | O T>C3 G3/2 |
| 24 | `1/4 + 1/8*sqrt2` | no | 4 | Oh C3v>C1h G3/2u |
| 24 | `1/4 - 1/12*sqrt3` | yes | 32 | O D3>C3 G3/2 |
| 24 | `1/4 - 1/8*sqrt2` | no | 18 | O D3>C2 G3/2 |
| 60 | `1/3 + 1/15*sqrt5` | no | 12 | I T>C3 T1 |
| 60 | `1/3 - 1/15*sqrt5` | no | 12 | I T>C3 T2 |
| 60 | `1/9 + 2/45*zeta^2 + 1/45*zeta^4 + 2/45*zeta^6 + 4/45*zeta^8 + 1/45*zeta^10 - 1/15*zeta^12 - 1/15*zeta^14` | no | 10 | I T>C3 I5/2 |
| 60 | `11/15 + 2/15*zeta^2 + 1/15*zeta^4 - 2/15*zeta^6 - 4/15*zeta^8 - 1/5*zeta^10 + 1/15*zeta^12 + 1/3*zeta^14` | no | 12 | I T>C3 E1/2 |
| 60 | `11/30 + 1/15*zeta^2 - 1/15*zeta^4 - 2/15*zeta^6 - 1/15*zeta^8 - 1/15*zeta^10 + 1/5*zeta^14` | no | 10 | I T>C3 G3/2 |
| 60 | `2/9 - 2/45*zeta^2 - 1/45*zeta^4 - 2/45*zeta^6 - 4/45*zeta^8 - 1/45*zeta^10 + 1/15*zeta^12 + 1/15*zeta^14` | no | 10 | I T>C3 I5/2 |
| 60 | `7/15 + 4/15*zeta^2 + 2/15*zeta^4 + 1/15*zeta^6 + 2/15*zeta^8 - 1/15*zeta^10 - 1/5*zeta^12` | no | 12 | I T>C3 E7/2 |
| 60 | `7/30 - 1/15*zeta^2 - 2/15*zeta^4 - 1/15*zeta^6 + 1/15*zeta^8 + 1/15*zeta^10` | no | 12 | I T>C3 G3/2 |
| 60 | `8/45 - 4/45*zeta^2 - 2/45*zeta^4 - 1/45*zeta^6 - 2/45*zeta^8 + 1/45*zeta^10 + 1/15*zeta^12` | no | 12 | I T>C3 I5/2 |
| 120 | `1/2 + 1/10*zeta^2 + 2/5*zeta^6 - 3/10*zeta^18 - 1/10*zeta^22 + 1/5*zeta^30` | no | 23 | I T>D2>C2 E1/2 |
| 120 | `1/2 + 3/10*zeta^2 + 1/5*zeta^6 + 1/10*zeta^18 - 3/10*zeta^22 + 1/10*zeta^30` | no | 3 | Ih Th>C2v>C1h E1/2g |
| 120 | `1/2 - 1/10*zeta^2 - 2/5*zeta^6 + 3/10*zeta^18 + 1/10*zeta^22 - 1/5*zeta^30` | no | 3 | Ih Th>C2v>C1h E7/2u |
| 120 | `1/2 - 3/10*zeta^2 - 1/5*zeta^6 - 1/10*zeta^18 + 3/10*zeta^22 - 1/10*zeta^30` | no | 23 | I T>D2>C2 E7/2 |
| 120 | `1/4 + 1/10*zeta^2 + 1/10*zeta^6 - 1/10*zeta^22 + 1/20*zeta^30` | no | 2 | Ih Th>C2v>C1h G3/2u |
| 120 | `1/4 + 1/10*zeta^6 - 1/10*zeta^18 + 1/20*zeta^30` | no | 1 | Ih Th>D2h>C1h G3/2u |
| 120 | `1/4 - 1/10*zeta^2 - 1/10*zeta^6 + 1/10*zeta^22 - 1/20*zeta^30` | no | 2 | Ih Th>C2v>C1h G3/2g |
| 120 | `1/4 - 1/10*zeta^6 + 1/10*zeta^18 - 1/20*zeta^30` | no | 21 | I T>D2>C2 G3/2 |
| 120 | `1/6 + 1/10*zeta^2 + 1/15*zeta^6 + 1/30*zeta^18 - 1/10*zeta^22 + 1/30*zeta^30` | no | 21 | I T>D2>C2 I5/2 |
| 120 | `1/6 + 1/10*zeta^2 + 2/15*zeta^6 - 1/30*zeta^18 - 1/10*zeta^22 + 1/15*zeta^30` | no | 1 | Ih Th>D2h>C1h I5/2u |
| 120 | `1/6 + 1/30*zeta^2 + 2/15*zeta^6 - 1/10*zeta^18 - 1/30*zeta^22 + 1/15*zeta^30` | no | 2 | Ih Th>C2v>C1h I5/2u |
| 120 | `1/6 + 1/30*zeta^2 - 1/15*zeta^6 + 1/10*zeta^18 - 1/30*zeta^22 - 1/30*zeta^30` | no | 2 | Ih Th>C2v>C1h I5/2u |
| 120 | `1/6 - 1/10*zeta^2 - 1/15*zeta^6 - 1/30*zeta^18 + 1/10*zeta^22 - 1/30*zeta^30` | no | 1 | Ih Th>D2h>C1h I5/2u |
| 120 | `1/6 - 1/10*zeta^2 - 2/15*zeta^6 + 1/30*zeta^18 + 1/10*zeta^22 - 1/15*zeta^30` | no | 19 | I T>D2>C2 I5/2 |
| 120 | `1/6 - 1/30*zeta^2 + 1/15*zeta^6 - 1/10*zeta^18 + 1/30*zeta^22 + 1/30*zeta^30` | no | 2 | Ih Th>C2v>C1h I5/2g |
| 120 | `1/6 - 1/30*zeta^2 - 2/15*zeta^6 + 1/10*zeta^18 + 1/30*zeta^22 - 1/15*zeta^30` | no | 2 | Ih Th>C2v>C1h I5/2g |
| 120 | `1/7` | yes | 0 |  |
| 120 | `2/7` | yes | 0 |  |
| 120 | `7/10` | yes | 0 |  |
| 120 | `7/20` | yes | 0 |  |
