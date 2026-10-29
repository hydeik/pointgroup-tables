"""Prototype: chain-adapted bases, the ``subgroup_chain`` table and the ``cg``
document, exactly over the cyclotomic field (wayfinder ticket #50).

Throwaway reference code, never shipped, like the version-1 prototypes. It
drives the modules ``chain_prototype_exact`` (exact real square root, norm
factors), ``chain_prototype_basis`` (chain classes, embedding, label spaces,
chain-adapted matrices, ``subgroup_chain``), ``chain_prototype_cg`` (coupling
matrices) and ``chain_prototype_papers`` (the printed values of the papers):

    uv run python3 prototypes/chain_cg_prototype.py chains [GROUP ...] [--out SUFFIX]
        every accepted chain class of the groups (default: the seven polyhedral
        groups and D6h): conductor, squared norms, the checks of the chapter;
        writes prototypes/chain_prototype_conductors.json,
        prototypes/chain_prototype_norm_factors.json (the shipped data of #44, the
        frozen data the chapter group-chains.md points at) and the generated
        prototypes/chain_prototype_summary.md (counts per group and per entry).
        Runs are long (T Td Th O Oh I D6h about 13 minutes, Ih about four hours,
        one process each); --out SUFFIX writes chain_prototype_*_SUFFIX.json instead of
        the frozen files, for a parallel or partial run, and
    uv run python3 prototypes/chain_cg_prototype.py merge SUFFIX
        merges those rows into the frozen files and rewrites the summary.
        --order sequence / --designated first reproduce the rules before ticket
        #58 (suffixes _sequence / _first) for comparison only.
    uv run python3 prototypes/chain_cg_prototype.py summary
        rewrites prototypes/chain_prototype_summary.md from the frozen files
    uv run python3 prototypes/chain_cg_prototype.py twins [GROUP ...] [--rewrite]
        the twin chain classes (ticket #59: accepted classes sharing one canonical
        spelling, of which the embedding rule reaches one); --rewrite brings the
        frozen data to that rule without recomputing a basis (rows respelled, the
        unreached twins dropped, norm-factor uses recomputed for the chains touched)
    uv run python3 prototypes/chain_cg_prototype.py matrices GROUP CHAIN
        a sample ``matrices --chain`` document (JSON) in prototypes/samples/
    uv run python3 prototypes/chain_cg_prototype.py cg GROUP [--chain CHAIN] [--irreps all|ordinary|projective] [--no-document]
        the coupling matrices of every pair, the checks CG-1 ... CG-5, timing and
        a sample ``cg`` document (JSON) in prototypes/samples/
    uv run python3 prototypes/chain_cg_prototype.py compare
        the T > C3 and O > D4 > C4 matrices and the T-dagger Clebsch-Gordan
        coefficients against the papers' printed values
    uv run python3 prototypes/chain_cg_prototype.py all
        chains, the sample documents and the comparison (not the Ih timing)

Nothing here is production code.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))

import chain_prototype_basis as cb
import chain_prototype_cg as cgm
import chain_prototype_papers as papers
from chain_prototype_exact import NormFactorTable, real_sqrt

from pointgroup_tables.document import (
    Request,
    build_document,
    build_model,
)
from pointgroup_tables.field import (
    CyclotomicField,
    cyclotomic_field,
    field_metadata,
)
from pointgroup_tables.formats import to_json
from pointgroup_tables.symbols import parse_group
from pointgroup_tables.values import numeric, parse, pretty

SAMPLES = ROOT / "prototypes" / "samples"  # sample documents (gitignored but C3v)
DATA = (
    ROOT / "prototypes"
)  # the generated conductor rows and norm-factor entries (tracked)
POLYHEDRAL = ["T", "Td", "Th", "O", "Oh", "I", "Ih"]

# the printed entries of docs/spec/group-chains.md, "The canonical norm factor" (normative)
PRINTED_ENTRIES = [
    (24, "(3 - sqrt3)/6", "(3 - sqrt3 + i*sqrt6)/6"),
    (24, "(3 + sqrt3)/6", "(3 + sqrt3 + i*sqrt6)/6"),
    (24, "(3 - sqrt3)/12", "(3*sqrt2 - sqrt6 + 2*i*sqrt3)/12"),
    (24, "(3 + sqrt3)/12", "(3*sqrt2 + sqrt6 + 2*i*sqrt3)/12"),
    (120, "1/7", "(2 + i*sqrt3)/7"),
    (120, "2/7", "(2*sqrt2 + i*sqrt6)/7"),
    (120, "7/10", "(2*sqrt10 + i*sqrt30)/10"),
    (120, "7/20", "(2*sqrt5 + i*sqrt15)/10"),
]


def _value_text(F: CyclotomicField, text: str):
    """The value grammar of the chapter with ``sqrtN`` for composite ``N``."""
    for n, spelled in (
        (30, "sqrt2*sqrt3*sqrt5"),
        (15, "sqrt3*sqrt5"),
        (10, "sqrt2*sqrt5"),
        (6, "sqrt2*sqrt3"),
    ):
        text = text.replace(f"sqrt{n}", spelled)
    return parse(F, text)


def seeded_table(generated: bool = True) -> NormFactorTable:
    """The printed entries of the chapter, plus the entries a previous ``chains``
    run wrote to ``prototypes/samples/chain_prototype_norm_factors.json`` (each
    verified ``mu * conj(mu) = N`` again when loaded)."""
    table: NormFactorTable = {}
    for M, N_text, mu_text in PRINTED_ENTRIES:
        F = cyclotomic_field(M)
        N, mu = _value_text(F, N_text), _value_text(F, mu_text)
        assert F.norm_squared(mu) == N, (M, N_text)
        assert real_sqrt(F, N) is None, (M, N_text)
        table[(M, N)] = mu
    path = DATA / "chain_prototype_norm_factors.json"
    if generated and path.exists():
        for entry in json.loads(path.read_text()):
            F = cyclotomic_field(entry["conductor"])
            N, mu = (
                F.from_exact(entry["N"]["exact"]),
                F.from_exact(entry["mu"]["exact"]),
            )
            assert F.norm_squared(mu) == N
            table.setdefault((entry["conductor"], N), mu)
    return table


# ---- chains ------------------------------------------------------------------------------


def _record_usage(
    symbol: str,
    text: str,
    resolved: cb.ResolvedChain,
    table: NormFactorTable,
    usage: dict[tuple[int, Any], set[str]],
) -> tuple[int, list[Any]]:
    """The conductor of a reachable chain and the norm-factor entries it and its
    sub-chains use (``usage``: ``(M, N) -> {"G text [sub-chain] irrep"}``)."""
    M, norms = cb.chain_conductor(resolved, table)
    FL = cyclotomic_field(resolved.chain_field)
    for irrep, label, N in norms:
        NM = FL.lift(N, cyclotomic_field(M))
        if real_sqrt(cyclotomic_field(M), NM) is None:
            usage.setdefault((M, NM), set()).add(f"{symbol} {text} {irrep}")
    for sub in _sub_chains(resolved):
        Ms, sub_norms = cb.chain_conductor(sub, table)
        FLs = cyclotomic_field(sub.chain_field)
        for irrep, label, N in sub_norms:
            NM = FLs.lift(N, cyclotomic_field(Ms))
            if real_sqrt(cyclotomic_field(Ms), NM) is None:
                usage.setdefault((Ms, NM), set()).add(
                    f"{symbol} {text} (sub-chain {sub.G.spec.symbol} {sub.text}) {irrep}"
                )
    return M, norms


def _subgroup_text(model: Any, S: frozenset[int]) -> str:
    return "{" + ", ".join(model.labels[i].text for i in sorted(S) if i != 0) + "}"


def twin_rows(symbol: str) -> list[dict[str, Any]]:
    """The twin chain classes of a group (ticket #59): per shared spelling the
    reached class, the unreached ones (the step at which each differs and its
    subgroup there) and the label-space agreement with the reached class."""
    model = build_model(parse_group(symbol))
    rows = []
    for twins in cb.twin_sets(symbol):
        reached = twins.classes[twins.reached]
        reached_resolved = cb.resolve(symbol, reached, twins.spelling)
        unreached = []
        for k, c in enumerate(twins.classes):
            if k == twins.reached:
                continue
            step = cb.differing_step(model, reached, c)
            other = cb.resolve(symbol, c, twins.spelling)
            unreached.append(
                {
                    "step": step,
                    "symbol": c.symbols[step],
                    "subgroup": _subgroup_text(model, other.steps[step - 1].members),
                    "reached_subgroup": _subgroup_text(
                        model, reached_resolved.steps[step - 1].members
                    ),
                    "label_spaces": cb.label_space_agreement(reached_resolved, other),
                }
            )
        rows.append(
            {
                "chain": twins.spelling,
                "classes": len(twins.classes),
                "symbols": reached.symbols,
                "unreached": unreached,
            }
        )
    return rows


def run_chains(
    groups: list[str], table: NormFactorTable, out: list[str]
) -> dict[str, Any]:
    conductors: dict[str, Any] = {}
    usage: dict[tuple[int, Any], set[str]] = {}
    for symbol in groups:
        chains = cb.reachable_chains(symbol)
        default = cb.default_chain(symbol)
        rows = []
        out.append(
            f"\n## list chains {symbol}: {len(chains)} spellings of {len(cb.accepted_chains(symbol))} chain classes, default {default!r}"
        )
        for text, chain in chains:
            started = time.perf_counter()
            resolved = cb.resolve(symbol, chain, text)
            M, norms = _record_usage(symbol, text, resolved, table, usage)
            basis = cb.chain_basis(resolved, table)
            failures = cb.check_basis(basis)
            FL = cyclotomic_field(resolved.chain_field)
            distinct = sorted({pretty(FL, N) for _, _, N in norms})
            rows.append(
                {
                    "chain": text,
                    "symbols": chain.symbols,
                    "default": text == default,
                    "conductor": M,
                    "chain_field": resolved.chain_field,
                }
            )
            status = "ok" if not failures else "; ".join(failures[:3])
            out.append(
                f"  {text:22s} L={resolved.chain_field:3d} M_chain={M:3d}  norms={distinct}  checks={status}  [{time.perf_counter() - started:.1f}s]"
            )
        conductors[symbol] = {
            "default": default,
            "chains": rows,
            "twins": twin_rows(symbol),
        }
    entries = _entries(table, usage)
    return {"conductors": conductors, "norm_factors": entries}


def _entries(
    table: NormFactorTable, usage: dict[tuple[int, Any], set[str]]
) -> list[dict[str, Any]]:
    entries = []
    for (M, N), mu in sorted(table.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        F = cyclotomic_field(M)
        assert F.norm_squared(mu) == N
        where = sorted(usage.get((M, N), set()))
        printed = any(
            M == pm and N == _value_text(F, pn) for pm, pn, _ in PRINTED_ENTRIES
        )
        entries.append(
            {
                "conductor": M,
                "N": {"pretty": pretty(F, N), "exact": F.exact(N)},
                "mu": {
                    "pretty": pretty(F, mu),
                    "exact": F.exact(mu),
                    "numeric": numeric(F, mu),
                },
                "printed": printed,
                "where": where,
            }
        )
    return entries


def run_twins(
    groups: list[str], table: NormFactorTable, out: list[str], rewrite: bool
) -> None:
    """List the twin chain classes of the groups; with ``rewrite``, bring the frozen
    data to the rule of ticket #59 without recomputing any basis: rows respelled by
    the chapter's spelling rule, the unreached twins' rows dropped, the twins
    recorded, and the ``where`` lists of the norm-factor entries recomputed for
    every chain whose spelling changed or has twins (the entries only an unreached
    twin used are dropped unless printed)."""
    conductors_path, norm_factors_path, _ = data_paths("")
    conductors = json.loads(conductors_path.read_text())
    entries = json.loads(norm_factors_path.read_text())
    usage: dict[tuple[int, Any], set[str]] = {}
    dirty_prefixes: list[str] = []
    for symbol in groups:
        twins = twin_rows(symbol)
        out.append(
            f"\n## {symbol}: {len(cb.accepted_chains(symbol))} chain classes, {len(cb.reachable_chains(symbol))} spellings, {len(twins)} shared"
        )
        for row in twins:
            for u in row["unreached"]:
                ls = u["label_spaces"]
                out.append(
                    f"  {row['chain']:28s} step {u['step']} {u['symbol']}: reached {u['reached_subgroup']}, unreached {u['subgroup']}; "
                    f"same lines: ordinary {len(ls['same_ordinary'])}/{len(ls['same_ordinary']) + len(ls['different_ordinary'])}, "
                    f"projective {len(ls['same_projective'])}/{len(ls['same_projective']) + len(ls['different_projective'])}"
                )
        if not rewrite or symbol not in conductors:
            continue
        record = conductors[symbol]
        old_rows = record["chains"]
        classes = cb.accepted_chains(symbol)
        if len(old_rows) != len(classes):
            raise ValueError(
                f"{symbol}: {len(old_rows)} rows for {len(classes)} classes"
            )
        reachable = {id(c): text for text, c in cb.reachable_chains(symbol)}
        twin_texts = {row["chain"] for row in twins}
        new_rows = []
        for old, c in zip(old_rows, classes, strict=True):
            if old["symbols"] != c.symbols:
                raise ValueError(f"{symbol}: row order differs from the class order")
            text = reachable.get(id(c))
            if text is None:
                dirty_prefixes.append(f"{symbol} {old['chain']} ")
                continue
            if text != old["chain"] or text in twin_texts:
                dirty_prefixes.append(f"{symbol} {old['chain']} ")
                started = time.perf_counter()
                M, _ = _record_usage(
                    symbol, text, cb.resolve(symbol, c, text), table, usage
                )
                if M != old["conductor"]:
                    raise ValueError(
                        f"{symbol} {text}: conductor {M} != {old['conductor']}"
                    )
                out.append(
                    f"  recomputed uses of {text} [{time.perf_counter() - started:.1f}s]"
                )
            new_rows.append({**old, "chain": text})
        new_rows.sort(key=lambda r: (len(r["symbols"]), r["chain"]))
        record["chains"] = new_rows
        record["twins"] = twins
    if not rewrite:
        return
    index: dict[tuple[int, str], dict[str, Any]] = {}
    for e in entries:
        e["where"] = [
            w for w in e["where"] if not any(w.startswith(p) for p in dirty_prefixes)
        ]
        index[(e["conductor"], json.dumps(e["N"]["exact"]))] = e
    for (M, N), where in usage.items():
        F = cyclotomic_field(M)
        key = (M, json.dumps(F.exact(N)))
        if key not in index:
            continue  # N = 1/2 (factor (1 + i)/2) needs no entry, as in run_chains
        index[key]["where"] = sorted(set(index[key]["where"]) | where)
    kept = [e for e in entries if e["where"] or e["printed"]]
    out.append(
        f"\nrewrote the frozen data: {sum(len(r['chains']) for r in conductors.values())} spellings, {len(kept)} entries ({len(entries) - len(kept)} dropped as used only by unreached twins)"
    )
    write_data({"conductors": conductors, "norm_factors": kept}, "")


def data_paths(suffix: str) -> tuple[Path, Path, Path]:
    """The conductor, norm-factor and summary files of a run (``""``: the frozen data)."""
    return (
        DATA / f"chain_prototype_conductors{suffix}.json",
        DATA / f"chain_prototype_norm_factors{suffix}.json",
        DATA / f"chain_prototype_summary{suffix}.md",
    )


def write_data(data: dict[str, Any], suffix: str) -> Path:
    """Write the two data files of a run and its summary; returns the summary path."""
    conductors_path, norm_factors_path, _ = data_paths(suffix)
    conductors_path.write_text(json.dumps(data["conductors"], indent=1) + "\n")
    norm_factors_path.write_text(json.dumps(data["norm_factors"], indent=1) + "\n")
    return write_summary(suffix)


def write_summary(suffix: str) -> Path:
    """Generate the Markdown summary of a run's data files: the counts the chapter
    does not state (chain classes and conductors per group, norm-factor entries
    and their uses), so that a regeneration cannot leave the chapter stale."""
    conductors_path, norm_factors_path, summary_path = data_paths(suffix)
    conductors = json.loads(conductors_path.read_text())
    entries = json.loads(norm_factors_path.read_text())
    lines = [
        "# Chain prototype data: summary",
        "",
        f"Generated from `{conductors_path.name}` and `{norm_factors_path.name}` by",
        "`uv run python3 prototypes/chain_cg_prototype.py chains` (or `summary`);",
        "never edited by hand. `docs/spec/group-chains.md` points here for every",
        "count.",
        "",
        "## Chain classes and conductors per group",
        "",
        "`Spellings` counts the entries of `list chains` (one per canonical",
        "spelling, the reached chain); `Classes` the accepted chain classes, which",
        "exceed the spellings by the unreached twin chain classes (ticket #59).",
        "",
        "| Group | Spellings | Classes | Default | `M_chain` × spellings |",
        "| ----- | --------: | ------: | ------- | --------------------- |",
    ]
    printed_count = sum(1 for e in entries if e["printed"])
    total_classes = 0
    total_spellings = 0
    twin_lines = []
    for group, record in conductors.items():
        rows = record["chains"]
        twins = record.get("twins", [])
        unreached = sum(len(t["unreached"]) for t in twins)
        total_spellings += len(rows)
        total_classes += len(rows) + unreached
        counts: dict[int, int] = {}
        for row in rows:
            counts[row["conductor"]] = counts.get(row["conductor"], 0) + 1
        distribution = ", ".join(f"{M} × {n}" for M, n in sorted(counts.items()))
        lines.append(
            f"| `{group}` | {len(rows)} | {len(rows) + unreached} | `{record['default']}` | {distribution} |"
        )
        for t in twins:
            for u in t["unreached"]:
                ls = u["label_spaces"]
                agreement = (
                    f"ordinary {len(ls['same_ordinary'])}/{len(ls['same_ordinary']) + len(ls['different_ordinary'])}, "
                    f"projective {len(ls['same_projective'])}/{len(ls['same_projective']) + len(ls['different_projective'])}"
                )
                twin_lines.append(
                    f"| `{group}` | `{t['chain']}` | {u['step']} `{u['symbol']}` | `{u['reached_subgroup']}` | `{u['subgroup']}` | {agreement} |"
                )
    lines += [
        "",
        f"{total_spellings} spellings ({total_classes} chain classes) in {len(conductors)} groups.",
        "",
        "## Twin chain classes",
        "",
        "Accepted chain classes that share one canonical spelling (ticket #59): the",
        "embedding rule reaches the first, `--chain` cannot name the other. `Step`",
        "is the step at which the unreached class differs, with the subgroup each",
        "class holds there; `same lines` counts the irreps whose label spaces",
        "coincide under the two chains.",
        "",
        "| Group | Spelling | Step | Reached | Unreached | Same lines |",
        "| ----- | -------- | ---- | ------- | --------- | ---------- |",
        *twin_lines,
        "",
        "## Norm-factor entries",
        "",
        f"{len(entries)} entries `(M, N) ↦ μ`, {printed_count} of them printed in",
        "the chapter; `uses` counts the chain classes and",
        "irreps (sub-chains included) whose designated squared norm is `N`.",
        "",
        "| `M` | `N` | printed | uses | first use |",
        "| --: | --- | ------- | ---: | --------- |",
    ]
    for e in entries:
        first = e["where"][0] if e["where"] else ""
        lines.append(
            f"| {e['conductor']} | `{e['N']['pretty']}` | {'yes' if e['printed'] else 'no'} | {len(e['where'])} | {first} |"
        )
    summary_path.write_text("\n".join(lines) + "\n")
    return summary_path


def merge_data(suffix: str) -> str:
    """Merge the rows of a suffixed run into the frozen files (groups replaced,
    norm-factor uses united, new entries appended) and rewrite the summary."""
    main_c, main_n, _ = data_paths("")
    add_c, add_n, _ = data_paths(suffix)
    conductors = json.loads(main_c.read_text())
    conductors.update(json.loads(add_c.read_text()))
    entries = json.loads(main_n.read_text())
    index = {(e["conductor"], json.dumps(e["N"]["exact"])): e for e in entries}
    added = 0
    for e in json.loads(add_n.read_text()):
        key = (e["conductor"], json.dumps(e["N"]["exact"]))
        if key in index:
            target = index[key]
            if target["mu"]["exact"] != e["mu"]["exact"]:
                raise ValueError(f"two norm factors for {key}")
            target["where"] = sorted(set(target["where"]) | set(e["where"]))
        else:
            entries.append(e)
            index[key] = e
            added += 1
    entries.sort(key=lambda e: (e["conductor"], e["N"]["pretty"]))
    write_data({"conductors": conductors, "norm_factors": entries}, "")
    return f"merged {suffix}: groups {sorted(conductors)}, {len(entries)} entries ({added} new)"


def _sub_chains(chain: cb.ResolvedChain) -> list[cb.ResolvedChain]:
    out = []
    while len(chain.steps) > 1:
        chain = cb.pulled_back(chain)
        out.append(chain)
    return out


# ---- documents -----------------------------------------------------------------------------


def _reorder_metadata(
    doc: dict[str, Any],
    chain_meta: dict[str, Any] | None,
    field: dict[str, Any],
    tables: list[str],
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in doc.items():
        if key == "tables":
            out[key] = tables
        elif key == "field":
            out[key] = field
        else:
            out[key] = value
        if key == "basis_convention" and chain_meta is not None:
            out["chain"] = chain_meta
    return out


def matrices_document(
    symbol: str,
    chain_text: str,
    forms: tuple[str, ...],
    table: NormFactorTable,
    which: str = "all",
) -> tuple[dict[str, Any], cb.ChainBasis]:
    resolved = cb.resolve_text(symbol, chain_text)
    basis = cb.chain_basis(resolved, table)
    F = basis.field
    base = build_document(Request("matrices", symbol, values=forms, irreps=which))
    doc = _reorder_metadata(
        base,
        cb.chain_metadata(basis, chain_text),
        field_metadata(basis.conductor),
        ["elements", "subgroup_chain", "matrices"],
    )
    entries = []
    for entry in doc["matrices"]["irreps"]:
        ci = basis.irreps[entry["index"]]
        new: dict[str, Any] = {}
        for key, value in entry.items():
            if key in ("exact", "numeric", "pretty"):
                continue
            new[key] = value
            if key == "frobenius_schur":
                new["components"] = cb.components_of(basis, ci)
        if "exact" in forms:
            new["exact"] = [
                [[F.exact(x) for x in row] for row in D] for D in ci.matrices
            ]
        if "numeric" in forms:
            new["numeric"] = [
                [[numeric(F, x) for x in row] for row in D] for D in ci.matrices
            ]
        if "pretty" in forms:
            new["pretty"] = [
                [[pretty(F, x) for x in row] for row in D] for D in ci.matrices
            ]
        entries.append(new)
    rebuilt: dict[str, Any] = {}
    for key, value in doc.items():
        if key == "elements":
            rebuilt[key] = value
            rebuilt["subgroup_chain"] = cb.subgroup_chain_table(resolved)
        elif key == "matrices":
            rebuilt[key] = {"irreps": entries}
        else:
            rebuilt[key] = value
    return rebuilt, basis


def cg_document(
    symbol: str,
    chain_text: str | None,
    which: str,
    forms: tuple[str, ...],
    table: NormFactorTable,
    out: list[str],
    write: bool = True,
) -> tuple[
    dict[str, Any] | None,
    dict[tuple[int, int], cgm.Pair],
    CyclotomicField,
    cb.ChainBasis | None,
]:
    model = build_model(parse_group(symbol))
    basis = None
    if chain_text is not None:
        resolved = cb.resolve_text(symbol, chain_text)
        basis = cb.chain_basis(resolved, table)
    M_cg = cgm.cg_conductor(model, None if basis is None else basis.conductor)
    Fcg = cyclotomic_field(M_cg)
    started = time.perf_counter()
    pairs: dict[tuple[int, int], cgm.Pair] = {}
    failures: list[str] = []
    norms = set()
    copies = 0
    chain_matrices = None
    if basis is not None:
        chain_matrices = [
            [
                [[basis.field.lift(x, Fcg) for x in row] for row in D]
                for D in ci.matrices
            ]
            for ci in basis.irreps
        ]
    for alpha, beta in cgm.select_pairs(model, which):
        pair = cgm.couple(model, alpha, beta, Fcg, table)
        failures += cgm.check_pair(model, pair, Fcg, table)
        if basis is not None:
            pair = cgm.transform(pair, Fcg, basis, model.irreps)
            failures += cgm.check_pair(model, pair, Fcg, table, chain_matrices)
        pairs[(alpha.index, beta.index)] = pair
        copies += len(pair.copies)
        norms |= {pretty(model.field, c.N) for c in pair.copies}
    elapsed = time.perf_counter() - started
    out.append(
        f"cg {symbol}{'' if chain_text is None else ' --chain ' + chain_text} --irreps {which}: {len(pairs)} pairs, "
        f"{copies} coupled copies, M_cg = {M_cg}, squared norms {sorted(norms, key=lambda s: (len(s), s))}, "
        f"checks {'ok' if not failures else failures[:5]}, {elapsed:.1f} s"
    )
    if not write:
        return None, pairs, Fcg, basis
    base = build_document(Request("matrices", symbol, values=forms, irreps=which))
    tables = (
        ["elements", "cg"] if basis is None else ["elements", "subgroup_chain", "cg"]
    )
    doc = _reorder_metadata(
        base,
        None if basis is None else cb.chain_metadata(basis, chain_text or ""),
        field_metadata(M_cg),
        tables,
    )
    rebuilt: dict[str, Any] = {}
    for key, value in doc.items():
        if key == "elements":
            rebuilt[key] = value
            if basis is not None:
                rebuilt["subgroup_chain"] = cb.subgroup_chain_table(basis.chain)
        elif key == "matrices":
            rebuilt["cg"] = cgm.cg_table(
                model, [pairs[k] for k in sorted(pairs)], Fcg, forms, basis
            )
        else:
            rebuilt[key] = value
    return rebuilt, pairs, Fcg, basis


def write_sample(name: str, doc: dict[str, Any], out: list[str]) -> Path:
    SAMPLES.mkdir(exist_ok=True)
    path = SAMPLES / name
    text = to_json(doc)
    path.write_text(text)
    out.append(
        f"wrote {path.relative_to(ROOT)} ({len(text) / 1024:.0f} KB, {text.count(chr(10))} lines)"
    )
    return path


# ---- compare -------------------------------------------------------------------------------


def run_compare(table: NormFactorTable, out: list[str]) -> None:
    # T > C3 matrices
    T = build_model(parse_group("T"))
    resolved = cb.resolve_text("T", "C3")
    basis = cb.chain_basis(resolved, table)
    principal = next(i for i, label in enumerate(T.labels) if label.text == "C3^1(1)")
    half_turns = [i for i, label in enumerate(T.labels) if label.text.startswith("C2")]
    out.append(
        "\n## T > C3 against [CF98a] eqs 37-38 (D(C2z); D(C3^1) diagonal by construction)"
    )
    paper = {
        name: (weights, {"C2z": rows}, True)
        for name, (weights, rows) in papers.T_C3_C2Z.items()
    }
    intertwiners: papers.Intertwiners = {}
    out += [
        "  " + line
        for line in papers.compare_matrices(
            basis,
            paper,
            {"F": "T", "E1/2": "E1/2", "E3/2": "2E3/2", "E-3/2": "1E3/2"},
            principal,
            {"C2z": half_turns[:1]},
            3,
            found=intertwiners,
        )
    ]
    # O > D4 > C4 (and O > C4) matrices
    O = build_model(parse_group("O"))
    principal = next(i for i, label in enumerate(O.labels) if label.text == "C4^1(1)")
    coordinate = [
        i
        for i, label in enumerate(O.labels)
        if label.text in ("C2(1)", "C2(2)", "C2(3)")
    ]
    diagonal = [i for i, label in enumerate(O.labels) if label.text.startswith("C2'")]
    name_map = {
        "E": "E",
        "F0": "T1",
        "F2": "T2",
        "E1/2": "E1/2",
        "E3/2": "E5/2",
        "G3/2": "G3/2",
    }
    for chain_text in ("D4>C4", "C4"):
        resolved = cb.resolve_text("O", chain_text)
        basis = cb.chain_basis(resolved, table)
        out.append(
            f"\n## O > {chain_text} against [FC99] eqs 32-33 (C2y a coordinate half turn, C2c a diagonal one)"
        )
        out += [
            "  " + line
            for line in papers.compare_matrices(
                basis,
                papers.O_C4,
                name_map,
                principal,
                {"C2y": coordinate, "C2c": diagonal},
                4,
                papers.O_C4_UNCERTAIN,
            )
        ]
    # T-dagger Clebsch-Gordan coefficients in the T > C3 basis
    out.append(
        "\n## T-dagger Clebsch-Gordan coefficients in the T > C3 basis against [CF98c] Table III"
    )
    _, pairs, Fcg, basis = cg_document(
        "T", "C3", "all", ("exact",), table, out, write=False
    )
    assert basis is not None
    principal = next(i for i, label in enumerate(T.labels) if label.text == "C3^1(1)")
    out += [
        "  " + line
        for line in papers.compare_t_cg(T, basis, pairs, Fcg, principal, intertwiners)
    ]


# ---- main ----------------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("chains")
    p.add_argument("groups", nargs="*")
    p.add_argument("--order", choices=["principal", "sequence"], default="principal")
    p.add_argument("--designated", choices=["invariant", "first"], default="invariant")
    p.add_argument("--out", default="", help="suffix of the data files (parallel runs)")
    p = sub.add_parser("matrices")
    p.add_argument("group")
    p.add_argument("chain")
    p.add_argument("--values", default="exact,numeric")
    p = sub.add_parser("cg")
    p.add_argument("group")
    p.add_argument("--chain")
    p.add_argument("--irreps", default="all")
    p.add_argument("--values", default="exact,numeric")
    p.add_argument("--no-document", action="store_true")
    sub.add_parser("compare")
    sub.add_parser("all")
    sub.add_parser("summary")
    p = sub.add_parser("twins")
    p.add_argument("groups", nargs="*")
    p.add_argument(
        "--rewrite", action="store_true", help="bring the frozen data to the rule"
    )
    p = sub.add_parser("merge")
    p.add_argument("suffix")
    args = parser.parse_args(argv)
    out: list[str] = []
    table = seeded_table(generated=args.command != "chains")
    started = time.perf_counter()
    if args.command == "chains":
        cb.ORDER_RULE = args.order
        cb.DESIGNATED_RULE = args.designated
        data = run_chains(args.groups or POLYHEDRAL + ["D6h"], table, out)
        suffix = (
            ("" if args.order == "principal" else "_sequence")
            + ("" if args.designated == "invariant" else "_first")
            + ("_" + args.out if args.out else "")
        )
        summary = write_data(data, suffix)
        out.append(
            f"\n{len(data['norm_factors'])} norm-factor entries written to prototypes/chain_prototype_norm_factors{suffix}.json, summary in {summary.relative_to(ROOT)}"
        )
    elif args.command == "summary":
        out.append(f"summary written to {write_summary('').relative_to(ROOT)}")
    elif args.command == "twins":
        run_twins(args.groups or POLYHEDRAL + ["D6h"], table, out, args.rewrite)
    elif args.command == "merge":
        out.append(merge_data("_" + args.suffix))
    elif args.command == "matrices":
        forms = tuple(
            f for f in ("exact", "numeric", "pretty") if f in args.values.split(",")
        )
        doc, basis = matrices_document(args.group, args.chain, forms, table)
        failures = cb.check_basis(basis)
        out.append(
            f"matrices {args.group} --chain {args.chain}: M_chain = {basis.conductor}, checks {'ok' if not failures else failures[:5]}"
        )
        write_sample(
            f"{args.group}_matrices_chain_{args.chain.replace('>', '-').replace('[', '').replace(']', '').replace(chr(39), 'p')}.json",
            doc,
            out,
        )
    elif args.command == "cg":
        forms = tuple(
            f for f in ("exact", "numeric", "pretty") if f in args.values.split(",")
        )
        doc, _, _, _ = cg_document(
            args.group,
            args.chain,
            args.irreps,
            forms,
            table,
            out,
            write=not args.no_document,
        )
        if doc is not None:
            name = (
                f"{args.group}_cg"
                + (
                    ""
                    if args.chain is None
                    else "_chain_" + args.chain.replace(">", "-")
                )
                + ("" if args.irreps == "all" else "_" + args.irreps)
                + ".json"
            )
            write_sample(name, doc, out)
    elif args.command == "compare":
        run_compare(table, out)
    elif args.command == "all":
        data = run_chains(POLYHEDRAL + ["D6h"], table, out)
        (DATA / "chain_prototype_conductors.json").write_text(
            json.dumps(data["conductors"], indent=1) + "\n"
        )
        (DATA / "chain_prototype_norm_factors.json").write_text(
            json.dumps(data["norm_factors"], indent=1) + "\n"
        )
        for group, chain in (
            ("T", "C3"),
            ("T", "D2>C2"),
            ("O", "D4>C4"),
            ("D6h", "default"),
            ("I", "T>C3"),
        ):
            doc, basis = matrices_document(group, chain, ("exact", "numeric"), table)
            write_sample(
                f"{group}_matrices_chain_{chain.replace('>', '-')}.json", doc, out
            )
        for group, chain in (("T", None), ("T", "C3"), ("O", None), ("D6h", None)):
            doc, _, _, _ = cg_document(
                group, chain, "all", ("exact", "numeric"), table, out
            )
            assert doc is not None
            write_sample(
                f"{group}_cg"
                + ("" if chain is None else "_chain_" + chain.replace(">", "-"))
                + ".json",
                doc,
                out,
            )
        run_compare(table, out)
    out.append(f"\ntotal {time.perf_counter() - started:.1f} s")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
