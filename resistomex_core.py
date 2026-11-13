"""
resistomex_core.py - ResistomeX core logic (no Streamlit imports).

What this module does
---------------------
* Reads NATIVE output tables of AMRFinderPlus, ResFinder (acquired genes),
  ResFinder-PointFinder and CARD/RGI, each with its own parser.
* Validates the header against the tool's schema, so that a valid zero-hit
  result (right header, no rows) is told apart from an empty, mismatched or
  unparseable file.
* Keeps AMR, STRESS and VIRULENCE annotations apart, and keeps acquired genes
  apart from point mutations / variant models.
* Preserves isolate, organism, genome accession, tool, tool version and
  database version on every record.
* Compares tools ONLY within one isolate, and only when organism and genome
  accession are identical for all runs being compared.

What it does NOT do
-------------------
* It does not detect resistance and does not predict phenotype.
* It does not invent antibiotic-class mappings: classes are either
  reported by the tool itself or come from a mapping the user supplies.
* It does not compare point mutations across tools (nomenclature differs).
"""
from __future__ import annotations

import html
import io
import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import pandas as pd

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
TOOL_AMRFINDER = "AMRFinderPlus"
TOOL_RESFINDER = "ResFinder"
TOOL_POINTFINDER = "ResFinder-PointFinder"
TOOL_RGI = "CARD/RGI"
TOOLS = [TOOL_AMRFINDER, TOOL_RESFINDER, TOOL_POINTFINDER, TOOL_RGI]
AUTO = "Auto-detect from header"

STATUS_HITS = "VALID_HITS"
STATUS_ZERO = "VALID_ZERO_HIT"
STATUS_EMPTY = "INVALID_EMPTY_FILE"
STATUS_SCHEMA = "INVALID_UNRECOGNISED_HEADER"
STATUS_MISMATCH = "INVALID_TOOL_MISMATCH"
STATUS_PARSE = "INVALID_UNPARSEABLE"
STATUS_META = "INVALID_MISSING_ISOLATE_ID"
VALID_STATUSES = {STATUS_HITS, STATUS_ZERO}

STATUS_LABEL = {
    STATUS_HITS: "Valid - hits present",
    STATUS_ZERO: "Valid zero-hit result (correct header, no rows)",
    STATUS_EMPTY: "Invalid - file is completely empty (not a zero-hit result)",
    STATUS_SCHEMA: "Invalid - header does not match any supported tool",
    STATUS_MISMATCH: "Invalid - header belongs to a different tool than selected",
    STATUS_PARSE: "Invalid - file could not be parsed",
    STATUS_META: "Invalid - isolate ID missing",
}

STD_COLS = [
    "isolate_id", "organism", "genome_accession", "tool", "tool_version",
    "database_version", "source_file", "category", "detection_kind",
    "gene_symbol", "mutation", "element_symbol", "element_name",
    "reported_class", "reported_subclass", "reported_phenotype", "mechanism",
    "scope", "method", "evidence_level", "identity", "coverage", "contig",
    "start", "stop", "accession", "review_flag",
]

SCOPE_NOTES = {
    TOOL_AMRFINDER: (
        "AMRFinderPlus reports AMR, and (with the plus scope) STRESS and VIRULENCE "
        "elements, and organism-specific point mutations. Its output depends on the "
        "database version and on the --plus / --organism settings used."
    ),
    TOOL_RESFINDER: (
        "ResFinder acquired-gene output covers acquired resistance genes only. "
        "Chromosomal point mutations are reported by PointFinder (separate file, "
        "supported species only). Zero hits means no acquired gene passed the run's "
        "thresholds."
    ),
    TOOL_POINTFINDER: (
        "PointFinder reports curated chromosomal mutations for supported species only."
    ),
    TOOL_RGI: (
        "CARD/RGI reports CARD ARO terms from homolog and variant models, in "
        "Perfect / Strict / Loose tiers. CARD names usually differ from AMRFinderPlus "
        "and ResFinder gene symbols."
    ),
}

LIMITATIONS = [
    "ResistomeX organises results from other tools; it does not detect resistance.",
    "A detected gene or mutation does not prove phenotypic (clinical) resistance, "
    "expression, or treatment outcome.",
    "A zero-hit result means no hits were reported under that run's settings and "
    "database version; it does not prove that no resistance determinants exist.",
    "Tools differ in database, method and scope, so a gene reported by one tool and "
    "not another is not automatically a disagreement.",
    "Gene matching across tools uses normalised symbols. Different nomenclature "
    "(allele names, CARD ARO names) can hide true matches; family-level matches are "
    "flagged for manual verification, not accepted automatically.",
    "Compatibility of genomes is judged only by the organism and genome accession "
    "strings you supply (they must be identical); the application cannot check the "
    "sequences themselves, apart from a warning when contig IDs share nothing.",
    "Point mutations / variant models are listed per tool but are not compared "
    "across tools.",
    "Antibiotic classes are shown as reported by each tool or by a user-supplied "
    "mapping. Different tools use different class vocabularies and these are not "
    "merged automatically.",
    "The PointFinder and CARD/RGI parsers follow the documented column layout; "
    "check the testing notes for what has been validated on real output files.",
]

# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
_NULLS = {"", "NA", "N/A", "NAN", "NONE", "NULL", "-"}


def _norm(header) -> str:
    return re.sub(r"[\s_]+", " ", str(header).strip().lstrip("#").strip()).lower()


def _clean(value) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return None if text.upper() in _NULLS else text


def _num(value) -> Optional[float]:
    text = _clean(value)
    if text is None:
        return None
    try:
        return float(text.replace("%", ""))
    except ValueError:
        return None


def _pick(colmap: Dict[str, str], aliases: List[str]) -> Optional[str]:
    """Exact (normalised) header match only - no substring guessing."""
    for alias in aliases:
        if alias in colmap:
            return colmap[alias]
    return None


def gene_key(symbol: str) -> str:
    """Normalised symbol used for exact matching: tet(A) and tetA are the same."""
    return re.sub(r"[()\s]", "", str(symbol)).upper()


def _base_key(key: str) -> str:
    """Allele-stripped key, e.g. BLATEM-1B -> BLATEM. Used only for 'family-level'
    hints that must be verified manually."""
    return re.sub(r"[-_]?\d+[A-Z]?$", "", key)


_ORG_ALIASES = {
    "e coli": "Escherichia coli",
    "k pneumoniae": "Klebsiella pneumoniae",
    "s enterica": "Salmonella enterica",
}


def normalise_organism(text: str) -> Tuple[str, str]:
    """Return (display_name, species_key). Only abbreviation expansion and
    capitalisation are done; nothing is guessed."""
    raw = re.sub(r"\s+", " ", str(text or "")).strip()
    if not raw:
        return "", ""
    tokens = raw.split(" ")
    short = " ".join(tokens[:2]).lower().replace(".", "")
    if short in _ORG_ALIASES:
        full = _ORG_ALIASES[short]
        display = " ".join([full] + tokens[2:])
        return display, full.lower()
    if len(tokens) >= 2:
        display = " ".join([tokens[0].capitalize(), tokens[1].lower()] + tokens[2:])
        return display, f"{tokens[0].lower()} {tokens[1].lower().rstrip('.')}"
    return raw.capitalize(), raw.lower()


_ACC_RE = re.compile(
    r"^(GC[AF]_\d{9}\.\d+|(NC|NZ|NG|NW|AC)_[A-Z0-9]+(\.\d+)?|"
    r"[A-Z]{1,2}\d{5,8}(\.\d+)?|[A-Z]{4}\d{8,}(\.\d+)?)$"
)


def accession_looks_standard(acc: str) -> bool:
    return bool(_ACC_RE.match(str(acc or "").strip().upper()))


# --------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------
def _empty_records() -> pd.DataFrame:
    return pd.DataFrame(columns=STD_COLS)


@dataclass
class RunMeta:
    isolate_id: str = ""
    organism: str = ""
    genome_accession: str = ""
    tool_version: str = ""
    database_version: str = ""


@dataclass
class ParsedRun:
    filename: str
    selected_tool: str
    meta: RunMeta
    status: str = STATUS_PARSE
    message: str = ""
    detected_tool: Optional[str] = None
    tool: str = ""
    n_rows_read: int = 0
    n_rows_dropped: int = 0
    records: pd.DataFrame = field(default_factory=_empty_records)

    @property
    def is_valid(self) -> bool:
        return self.status in VALID_STATUSES

    @property
    def organism_display(self) -> str:
        return normalise_organism(self.meta.organism)[0]

    @property
    def organism_key(self) -> str:
        return normalise_organism(self.meta.organism)[1]


# --------------------------------------------------------------------------
# Reading and schema detection
# --------------------------------------------------------------------------
def _decode(content) -> str:
    if isinstance(content, str):
        return content
    for enc in ("utf-8-sig", "latin-1"):
        try:
            return content.decode(enc)
        except UnicodeDecodeError:
            continue
    raise ValueError("file encoding could not be decoded")


def _read_table(text: str) -> pd.DataFrame:
    """Read a delimited table. NA stays text, a leading '#' on the header line
    (older AMRFinder versions) is removed rather than treated as a comment."""
    lines = [ln for ln in text.splitlines() if not ln.startswith("##")]
    start = next((i for i, ln in enumerate(lines) if ln.strip()), None)
    if start is None:
        raise ValueError("no header line found")
    lines = lines[start:]
    lines[0] = lines[0].lstrip("#")
    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else "\t")
    df = pd.read_csv(
        io.StringIO("\n".join(lines)), sep=delim, dtype=str,
        keep_default_na=False, na_filter=False,
    )
    df.columns = [str(c).strip() for c in df.columns]
    return df


def detect_tool(columns) -> Optional[str]:
    cols = {_norm(c) for c in columns}
    if ({"element symbol", "gene symbol"} & cols
            and {"type", "element type"} & cols
            and {"method", "class"} & cols):
        return TOOL_AMRFINDER
    if {"best hit aro", "cut off", "model type"} <= cols:
        return TOOL_RGI
    if {"resistance gene", "identity", "coverage"} <= cols:
        return TOOL_RESFINDER
    if ("mutation" in cols and "resistance" in cols
            and {"nucleotide change", "amino acid change"} & cols):
        return TOOL_POINTFINDER
    return None


# --------------------------------------------------------------------------
# Tool-specific parsers (each returns (records, n_dropped))
# --------------------------------------------------------------------------
_AMRF = {
    "symbol": ["element symbol", "gene symbol"],
    "name": ["element name", "sequence name"],
    "scope": ["scope"],
    "type": ["type", "element type"],
    "subtype": ["subtype", "element subtype"],
    "class": ["class"],
    "subclass": ["subclass"],
    "method": ["method"],
    "coverage": ["% coverage of reference", "% coverage of reference sequence"],
    "identity": ["% identity to reference", "% identity to reference sequence"],
    "contig": ["contig id"],
    "start": ["start"],
    "stop": ["stop"],
    "accession": ["closest reference accession", "accession of closest sequence"],
}


def _parse_amrfinder(df: pd.DataFrame):
    cm = {_norm(c): c for c in df.columns}
    f = {k: _pick(cm, v) for k, v in _AMRF.items()}
    out, dropped = [], 0
    for r in df.to_dict("records"):
        def g(k):
            return _clean(r[f[k]]) if f.get(k) else None
        sym = g("symbol")
        if not sym:
            dropped += 1
            continue
        raw_type = (g("type") or "").upper()
        category = raw_type if raw_type in {"AMR", "STRESS", "VIRULENCE"} else "OTHER"
        subtype = (g("subtype") or "").upper()
        method = (g("method") or "").upper()
        is_point = subtype.startswith("POINT") or method.startswith("POINT")
        gene, mutation = sym, None
        if is_point and "_" in sym:
            gene, mutation = sym.split("_", 1)
        flags = []
        if "PARTIAL" in method or "INTERNAL_STOP" in method:
            flags.append(f"Method {method}: partial/truncated match - review")
        if category == "OTHER":
            flags.append(f"Element type '{raw_type or 'missing'}' not recognised")
        out.append({
            "category": category,
            "detection_kind": "point_mutation" if is_point else "gene",
            "gene_symbol": gene, "mutation": mutation, "element_symbol": sym,
            "element_name": g("name"), "reported_class": g("class"),
            "reported_subclass": g("subclass"), "scope": g("scope"),
            "method": g("method"), "identity": _num(r[f["identity"]]) if f["identity"] else None,
            "coverage": _num(r[f["coverage"]]) if f["coverage"] else None,
            "contig": g("contig"), "start": _num(r[f["start"]]) if f["start"] else None,
            "stop": _num(r[f["stop"]]) if f["stop"] else None,
            "accession": g("accession"), "review_flag": "; ".join(flags) or None,
        })
    return out, dropped


def _parse_resfinder(df: pd.DataFrame):
    cm = {_norm(c): c for c in df.columns}
    c_gene = _pick(cm, ["resistance gene"])
    c_id, c_cov = _pick(cm, ["identity"]), _pick(cm, ["coverage"])
    c_contig = _pick(cm, ["contig"])
    c_pos = _pick(cm, ["position in contig"])
    c_pheno = _pick(cm, ["phenotype"])
    c_acc = _pick(cm, ["accession no", "accession no."])
    out, dropped = [], 0
    for r in df.to_dict("records"):
        sym = _clean(r[c_gene])
        if not sym:
            dropped += 1
            continue
        nums = re.findall(r"\d+", r[c_pos]) if c_pos and r[c_pos] else []
        out.append({
            "category": "AMR", "detection_kind": "gene", "gene_symbol": sym,
            "element_symbol": sym, "reported_phenotype": _clean(r[c_pheno]) if c_pheno else None,
            "identity": _num(r[c_id]), "coverage": _num(r[c_cov]),
            "contig": _clean(r[c_contig]) if c_contig else None,
            "start": float(nums[0]) if len(nums) >= 1 else None,
            "stop": float(nums[1]) if len(nums) >= 2 else None,
            "accession": _clean(r[c_acc]) if c_acc else None,
        })
    return out, dropped


def _parse_pointfinder(df: pd.DataFrame):
    cm = {_norm(c): c for c in df.columns}
    c_mut = _pick(cm, ["mutation"])
    c_res = _pick(cm, ["resistance"])
    out, dropped = [], 0
    for r in df.to_dict("records"):
        mut = _clean(r[c_mut])
        if not mut:
            dropped += 1
            continue
        out.append({
            "category": "AMR", "detection_kind": "point_mutation",
            "gene_symbol": mut.split()[0], "mutation": mut, "element_symbol": mut,
            "reported_phenotype": _clean(r[c_res]) if c_res else None,
            "review_flag": "PointFinder parser follows documented columns",
        })
    return out, dropped


def _parse_rgi(df: pd.DataFrame):
    cm = {_norm(c): c for c in df.columns}
    f = {
        "gene": _pick(cm, ["best hit aro"]), "cut": _pick(cm, ["cut off"]),
        "ident": _pick(cm, ["best identities"]), "model": _pick(cm, ["model type"]),
        "snps": _pick(cm, ["snps in best hit aro"]), "drug": _pick(cm, ["drug class"]),
        "mech": _pick(cm, ["resistance mechanism"]), "fam": _pick(cm, ["amr gene family"]),
        "cov": _pick(cm, ["percentage length of reference sequence"]),
        "contig": _pick(cm, ["contig"]), "start": _pick(cm, ["start"]),
        "stop": _pick(cm, ["stop"]), "nudged": _pick(cm, ["nudged"]),
    }
    def g(r, k):
        return _clean(r[f[k]]) if f.get(k) else None
    out, dropped = [], 0
    for r in df.to_dict("records"):
        sym = g(r, "gene")
        if not sym:
            dropped += 1
            continue
        model = (g(r, "model") or "").lower()
        kind = "gene" if "homolog" in model else (
            "point_mutation" if "variant" in model else "other_variant")
        cut = g(r, "cut")
        flags = []
        if cut and cut.lower() == "loose":
            flags.append("CARD Cut_Off=Loose: lower-confidence hit")
        if (g(r, "nudged") or "").lower() == "true":
            flags.append("Loose hit nudged to Strict by RGI")
        out.append({
            "category": "AMR", "detection_kind": kind, "gene_symbol": sym,
            "mutation": g(r, "snps") if kind != "gene" else None,
            "element_symbol": sym, "element_name": g(r, "fam"),
            "reported_class": g(r, "drug"), "mechanism": g(r, "mech"),
            "method": g(r, "model"), "evidence_level": cut,
            "identity": _num(r[f["ident"]]) if f["ident"] else None,
            "coverage": _num(r[f["cov"]]) if f["cov"] else None,
            "contig": g(r, "contig"),
            "start": _num(r[f["start"]]) if f["start"] else None,
            "stop": _num(r[f["stop"]]) if f["stop"] else None,
            "review_flag": "; ".join(flags) or None,
        })
    return out, dropped


_PARSERS = {
    TOOL_AMRFINDER: _parse_amrfinder, TOOL_RESFINDER: _parse_resfinder,
    TOOL_POINTFINDER: _parse_pointfinder, TOOL_RGI: _parse_rgi,
}


# --------------------------------------------------------------------------
# Public parsing entry point
# --------------------------------------------------------------------------
def parse_run(content, filename: str, selected_tool: str, meta: RunMeta) -> ParsedRun:
    """Parse one result file. Never raises for bad user files; it returns a
    ParsedRun whose `status` explains what happened."""
    run = ParsedRun(filename=filename, selected_tool=selected_tool, meta=meta)
    try:
        text = _decode(content)
    except Exception as exc:  # noqa: BLE001
        run.status, run.message = STATUS_PARSE, f"Could not decode file: {exc}"
        return run
    if not text.strip():
        run.status = STATUS_EMPTY
        run.message = ("File has no header and no rows. A real zero-hit run still "
                       "writes a header; this looks like a failed or truncated run.")
        return run
    try:
        df = _read_table(text)
    except Exception as exc:  # noqa: BLE001
        run.status, run.message = STATUS_PARSE, f"Could not parse table: {exc}"
        return run

    run.n_rows_read = len(df)
    run.detected_tool = detect_tool(df.columns)
    if run.detected_tool is None:
        run.status = STATUS_SCHEMA
        run.message = ("Header does not match AMRFinderPlus, ResFinder, PointFinder "
                       "or CARD/RGI native output. Columns found: "
                       + ", ".join(df.columns[:12]) + (" ..." if len(df.columns) > 12 else ""))
        return run
    if selected_tool not in (None, "", AUTO) and selected_tool != run.detected_tool:
        run.status = STATUS_MISMATCH
        run.message = (f"Selected tool is {selected_tool}, but the header looks like "
                       f"{run.detected_tool} output.")
        return run
    run.tool = run.detected_tool
    if not meta.isolate_id.strip():
        run.status, run.message = STATUS_META, "Enter an isolate / sample ID for this file."
        return run

    rows, dropped = _PARSERS[run.tool](df)
    run.n_rows_dropped = dropped
    rec = pd.DataFrame(rows)
    for col in STD_COLS:
        if col not in rec.columns:
            rec[col] = None
    rec["isolate_id"] = meta.isolate_id.strip()
    rec["organism"] = run.organism_display
    rec["genome_accession"] = meta.genome_accession.strip()
    rec["tool"] = run.tool
    rec["tool_version"] = meta.tool_version.strip()
    rec["database_version"] = meta.database_version.strip()
    rec["source_file"] = filename
    rec = rec[STD_COLS].copy()
    for col in ("identity", "coverage", "start", "stop"):
        rec[col] = pd.to_numeric(rec[col], errors="coerce")
    run.records = rec if len(rec) else _empty_records()
    if len(rec):
        run.status = STATUS_HITS
        run.message = f"{len(rec)} record(s) parsed."
    else:
        run.status = STATUS_ZERO
        run.message = ("Header matches the tool and there are no result rows: a valid "
                       "zero-hit result."
                       if not dropped else "Rows present but none had a gene symbol.")
        if dropped:
            run.status = STATUS_PARSE
    return run


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
def validate_runs(runs: List[ParsedRun]) -> List[dict]:
    issues: List[dict] = []

    def add(sev, run, msg):
        issues.append({"severity": sev, "file": run.filename if run else "", "message": msg})

    for r in runs:
        if not r.is_valid:
            add("error", r, f"{STATUS_LABEL.get(r.status, r.status)}. {r.message}")
        if r.n_rows_dropped:
            add("warning", r, f"{r.n_rows_dropped} row(s) dropped: no gene/element symbol.")
        if not r.is_valid:
            continue
        if not r.meta.organism.strip():
            add("warning", r, "Organism missing: this run cannot be used in tool comparison.")
        if not r.meta.genome_accession.strip():
            add("warning", r, "Genome accession missing: this run cannot be used in tool comparison.")
        elif not accession_looks_standard(r.meta.genome_accession):
            add("info", r, "Genome ID is not a recognised accession pattern; accepted as free-text assembly ID.")
        if not r.meta.tool_version.strip() or not r.meta.database_version.strip():
            add("warning", r, "Tool version and/or database version not recorded; "
                              "reports will show 'not recorded'.")
    valid = [r for r in runs if r.is_valid]
    names = [r.filename for r in runs]
    for n in sorted({n for n in names if names.count(n) > 1}):
        add("warning", None, f"File name '{n}' was uploaded more than once.")
    by_iso: Dict[str, List[ParsedRun]] = {}
    for r in valid:
        by_iso.setdefault(r.meta.isolate_id.strip(), []).append(r)
    for iso, rs in by_iso.items():
        orgs = {r.organism_key for r in rs if r.organism_key}
        if len(orgs) > 1:
            add("error", None, f"Isolate ID '{iso}' is used with different organisms "
                               f"({', '.join(sorted(orgs))}). Use unique isolate IDs.")
        tools = [r.tool for r in rs]
        for t in sorted({t for t in tools if tools.count(t) > 1}):
            add("warning", None, f"Isolate '{iso}' has more than one {t} file; "
                                 "it will not be tool-compared until resolved.")
    return issues


# --------------------------------------------------------------------------
# Combining, classes, summaries
# --------------------------------------------------------------------------
def combine_runs(runs: List[ParsedRun]) -> pd.DataFrame:
    frames = [r.records for r in runs if r.is_valid and len(r.records)]
    if not frames:
        data = _empty_records()
        data["gene_key"] = pd.Series(dtype=str)
        return data
    data = pd.concat(frames, ignore_index=True)
    data["gene_key"] = data["gene_symbol"].map(gene_key)
    return data


def load_mapping(df: pd.DataFrame) -> Tuple[Dict[str, str], List[str]]:
    """User-supplied gene -> class mapping. Conflicting duplicates are excluded
    (and reported) instead of silently multiplying rows."""
    warnings: List[str] = []
    cols = {_norm(c): c for c in df.columns}
    if "gene symbol" not in cols or "antibiotic class" not in cols:
        raise ValueError("Mapping needs columns gene_symbol and antibiotic_class.")
    pairs: Dict[str, set] = {}
    for sym, cls in zip(df[cols["gene symbol"]], df[cols["antibiotic class"]]):
        s, c = _clean(sym), _clean(cls)
        if s and c:
            pairs.setdefault(gene_key(s), set()).add(c)
    mapping = {}
    for key, classes in pairs.items():
        if len(classes) == 1:
            mapping[key] = next(iter(classes))
        else:
            warnings.append(f"Mapping for '{key}' lists conflicting classes "
                            f"({', '.join(sorted(classes))}); entry ignored.")
    return mapping, warnings


def harmonise_classes(data: pd.DataFrame, mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Add class_label / class_source. Tool-reported classes are kept verbatim;
    a user mapping is used only where the tool gave none. No built-in mapping."""
    mapping = mapping or {}
    d = data.copy()
    labels, sources, notes = [], [], []
    for _, r in d.iterrows():
        reported = _clean(r["reported_class"])
        tokens = [t.strip() for t in re.split(r"[;,]", reported)] if reported else []
        tokens = [t for t in tokens if t]
        note = None
        if r["category"] != "AMR":
            label, src = f"({r['category']} annotation - not an antibiotic class)", "n/a"
        elif len(tokens) == 1:
            label, src = tokens[0], "tool-reported"
        elif len(tokens) > 1:
            label, src = "Multiple classes (tool-reported)", "tool-reported"
        elif r["gene_key"] in mapping:
            label, src = mapping[r["gene_key"]], "user mapping"
        else:
            label, src = "Not reported by tool", "none"
        if src == "tool-reported" and label.casefold() == "efflux":
            note = "EFFLUX is a mechanism-type term used by the tool, not a single antibiotic class."
        labels.append(label)
        sources.append(src)
        notes.append(note)
    d["class_label"], d["class_source"], d["class_note"] = labels, sources, notes
    return d


def run_table(runs: List[ParsedRun]) -> pd.DataFrame:
    rows = []
    for r in runs:
        rec = r.records
        def n(cat, kind=None):
            if not len(rec):
                return 0
            m = rec["category"] == cat
            if kind == "gene":
                m &= rec["detection_kind"] == "gene"
            elif kind == "mut":
                m &= rec["detection_kind"] != "gene"
            return int(m.sum())
        rows.append({
            "file": r.filename, "isolate_id": r.meta.isolate_id, "organism": r.organism_display,
            "genome_accession": r.meta.genome_accession, "tool": r.tool or r.selected_tool,
            "tool_version": r.meta.tool_version or "not recorded",
            "database_version": r.meta.database_version or "not recorded",
            "status": STATUS_LABEL.get(r.status, r.status), "rows_read": r.n_rows_read,
            "AMR genes": n("AMR", "gene"), "AMR mutations/variants": n("AMR", "mut"),
            "STRESS": n("STRESS"), "VIRULENCE": n("VIRULENCE"), "OTHER": n("OTHER"),
        })
    return pd.DataFrame(rows)


def summarize(runs: List[ParsedRun], data: pd.DataFrame) -> dict:
    valid = [r for r in runs if r.is_valid]
    amr = data[data["category"] == "AMR"] if len(data) else data
    genes = amr[amr["detection_kind"] == "gene"] if len(amr) else amr
    return {
        "files": len(runs),
        "valid_runs": len(valid),
        "zero_hit_runs": sum(r.status == STATUS_ZERO for r in runs),
        "invalid_runs": sum(not r.is_valid for r in runs),
        "isolates": len({r.meta.isolate_id.strip() for r in valid}),
        "organisms": len({r.organism_key for r in valid if r.organism_key}),
        "genomes": len({r.meta.genome_accession.strip() for r in valid if r.meta.genome_accession.strip()}),
        "amr_gene_detections": len(genes),
        "unique_amr_genes": int(genes["gene_key"].nunique()) if len(genes) else 0,
        "amr_mutations_variants": len(amr) - len(genes),
        "stress": int((data["category"] == "STRESS").sum()) if len(data) else 0,
        "virulence": int((data["category"] == "VIRULENCE").sum()) if len(data) else 0,
        "other": int((data["category"] == "OTHER").sum()) if len(data) else 0,
    }


# --------------------------------------------------------------------------
# Tool-to-tool comparison (same isolate, same organism, same genome only)
# --------------------------------------------------------------------------
@dataclass
class ComparisonGroup:
    isolate_id: str
    organism: str = ""
    genome_accession: str = ""
    tools: List[str] = field(default_factory=list)
    status: str = "NOT_COMPARED"
    reasons: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    presence: pd.DataFrame = field(default_factory=pd.DataFrame)
    pairwise: List[dict] = field(default_factory=list)
    run_summary: pd.DataFrame = field(default_factory=pd.DataFrame)


def _contigs_overlap(a: set, b: set) -> bool:
    """True if the two tools share a contig. CARD/RGI appends an ORF number to
    the contig name (contig_12), so 'x' and 'x_12' count as the same contig."""
    for x in a:
        for y in b:
            if x == y or x.startswith(y + "_") or y.startswith(x + "_"):
                return True
    return False


def compare_tools(runs: List[ParsedRun]) -> List[ComparisonGroup]:
    groups: Dict[str, List[ParsedRun]] = {}
    for r in runs:
        if r.meta.isolate_id.strip():
            groups.setdefault(r.meta.isolate_id.strip(), []).append(r)
    out: List[ComparisonGroup] = []
    for iso, rs in sorted(groups.items()):
        g = ComparisonGroup(isolate_id=iso)
        invalid = [r for r in rs if not r.is_valid]
        valid = [r for r in rs if r.is_valid]
        if invalid:
            g.warnings.append("Excluded invalid run(s): " + ", ".join(r.filename for r in invalid))
        pf = [r for r in valid if r.tool == TOOL_POINTFINDER]
        if pf:
            g.notes.append("PointFinder results are not part of the gene comparison "
                           "(mutations are not compared across tools).")
        valid = [r for r in valid if r.tool != TOOL_POINTFINDER]
        g.tools = sorted({r.tool for r in valid})
        g.organism = ", ".join(sorted({r.organism_display for r in valid if r.organism_display}))
        g.genome_accession = ", ".join(sorted({r.meta.genome_accession.strip() for r in valid if r.meta.genome_accession.strip()}))

        for t in g.tools:
            if sum(r.tool == t for r in valid) > 1:
                g.reasons.append(f"More than one {t} file for this isolate; resolve which one to use.")
        orgs = {r.organism_key for r in valid}
        if "" in orgs:
            g.reasons.append("Organism not supplied for at least one run.")
        elif len(orgs) > 1:
            g.reasons.append("Runs list different organisms; they cannot be the same isolate.")
        gens = {r.meta.genome_accession.strip().upper() for r in valid}
        if "" in gens:
            g.reasons.append("Genome accession not supplied for at least one run.")
        elif len(gens) > 1:
            g.reasons.append("Runs list different genome accessions; not comparable as the same genome.")
        if len(g.tools) < 2:
            g.reasons.append("At least two different tools with valid runs are needed.")
        if g.reasons:
            g.status = "NOT_COMPARED"
            out.append(g)
            continue

        sets: Dict[str, Dict[str, str]] = {}
        contigs: Dict[str, set] = {}
        for t in g.tools:
            rec = pd.concat([r.records for r in valid if r.tool == t], ignore_index=True)
            sub = rec[(rec["category"] == "AMR") & (rec["detection_kind"] == "gene")] if len(rec) else rec
            sets[t] = {gene_key(s): s for s in sub["gene_symbol"]} if len(sub) else {}
            contigs[t] = set(sub["contig"].dropna()) if len(sub) else set()
        with_contigs = [t for t in g.tools if contigs[t]]
        for i, a in enumerate(with_contigs):
            for b in with_contigs[i + 1:]:
                if not _contigs_overlap(contigs[a], contigs[b]):
                    g.warnings.append(f"{a} and {b} share no contig IDs; check that both ran on the same assembly.")

        rows = []
        for key in sorted(set().union(*[set(s) for s in sets.values()])):
            present = [t for t in g.tools if key in sets[t]]
            row = {"gene": next(sets[t][key] for t in present)}
            for t in g.tools:
                row[t] = key in sets[t]
            if len(present) == len(g.tools):
                level, hint = "Exact match (all tools)", ""
            elif len(present) >= 2:
                level, hint = "Exact match (some tools)", ""
            else:
                fam = [(t, sym) for t in g.tools if t not in present
                       for k2, sym in sets[t].items() if k2 != key and _base_key(k2) == _base_key(key)]
                if fam:
                    level = "Family-level only - verify manually"
                    hint = "; ".join(f"{t} reports {sym}" for t, sym in fam)
                else:
                    level, hint = "Reported by one tool only", ""
            row["match_level"], row["note"] = level, hint
            rows.append(row)
        g.presence = pd.DataFrame(rows, columns=["gene", *g.tools, "match_level", "note"])
        for i, a in enumerate(g.tools):
            for b in g.tools[i + 1:]:
                g.pairwise.append({
                    "tool_a": a, "tool_b": b,
                    "shared (exact)": len(set(sets[a]) & set(sets[b])),
                    f"only {a}": len(set(sets[a]) - set(sets[b])),
                    f"only {b}": len(set(sets[b]) - set(sets[a])),
                })
        g.run_summary = pd.DataFrame([{"tool": t, "unique AMR genes": len(sets[t]),
                                       "run status": STATUS_LABEL[next(r.status for r in valid if r.tool == t)]}
                                      for t in g.tools])
        for t in g.tools:
            g.notes.append(f"{t}: {SCOPE_NOTES[t]}")
        zero = [t for t in g.tools if not sets[t]]
        if zero:
            g.notes.append(f"{', '.join(zero)} reported zero AMR genes. Differences against a tool "
                           "with broader scope are scope differences, not proof of disagreement.")
        g.status = "COMPARED"
        out.append(g)
    return out


# --------------------------------------------------------------------------
# HTML report
# --------------------------------------------------------------------------
def _tbl(df: pd.DataFrame) -> str:
    if df is None or not len(df):
        return "<p><em>None.</em></p>"
    return df.fillna("").to_html(index=False, escape=True, border=0)


def build_report_html(runs: List[ParsedRun], data: pd.DataFrame,
                      comparisons: List[ComparisonGroup], issues: List[dict]) -> str:
    s = summarize(runs, data)
    e = html.escape
    show = ["isolate_id", "organism", "genome_accession", "tool", "tool_version",
            "database_version", "gene_symbol", "mutation", "detection_kind",
            "class_label" if "class_label" in data.columns else "reported_class",
            "reported_phenotype", "identity", "coverage", "evidence_level",
            "review_flag", "source_file"]
    amr = data[data["category"] == "AMR"] if len(data) else data
    amr_genes = amr[amr["detection_kind"] == "gene"] if len(amr) else amr
    amr_mut = amr[amr["detection_kind"] != "gene"] if len(amr) else amr
    other = data[data["category"] != "AMR"] if len(data) else data
    stat_rows = pd.DataFrame([
        ("Files uploaded", s["files"]), ("Valid runs", s["valid_runs"]),
        ("Valid zero-hit runs", s["zero_hit_runs"]), ("Invalid runs", s["invalid_runs"]),
        ("Isolates", s["isolates"]), ("Organisms", s["organisms"]), ("Genomes", s["genomes"]),
        ("AMR gene detections", s["amr_gene_detections"]),
        ("Unique AMR genes (normalised symbols)", s["unique_amr_genes"]),
        ("AMR mutation / variant records", s["amr_mutations_variants"]),
        ("STRESS annotations (not AMR)", s["stress"]),
        ("VIRULENCE annotations (not AMR)", s["virulence"]),
    ], columns=["Measure", "Value"])
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><title>ResistomeX report</title>",
        "<style>body{font-family:Arial,sans-serif;margin:32px;color:#17212b}"
        "h1{color:#146b78}table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13px}"
        "td,th{border:1px solid #ddd;padding:6px;text-align:left}th{background:#eef5f6}"
        ".warn{background:#fff6e5;padding:8px;border-left:4px solid #e0a100}</style></head><body>",
        "<h1>ResistomeX – AMR Gene Profiling, Harmonisation and Comparative Analysis</h1><p><em>Prototype v0.1 - report</em></p>",
        f"<p>Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>",
        "<p class='warn'>ResistomeX organises output from other AMR tools. Gene detection "
        "does not prove phenotypic resistance.</p>",
        "<h2>Summary</h2>", _tbl(stat_rows),
        "<h2>Input runs and metadata</h2>", _tbl(run_table(runs)),
        "<h2>Validation messages</h2>",
        _tbl(pd.DataFrame(issues)) if issues else "<p>No validation issues.</p>",
        "<h2>AMR genes</h2>", _tbl(amr_genes[[c for c in show if c in amr_genes.columns]] if len(amr_genes) else amr_genes),
        "<h2>AMR point mutations / variant models</h2>",
        _tbl(amr_mut[[c for c in show if c in amr_mut.columns]] if len(amr_mut) else amr_mut),
        "<h2>Non-AMR annotations (STRESS / VIRULENCE / other)</h2>",
        _tbl(other[[c for c in show if c in other.columns]] if len(other) else other),
        "<h2>Tool-to-tool comparison (same isolate and genome only)</h2>",
    ]
    if not comparisons:
        parts.append("<p>No isolates.</p>")
    for g in comparisons:
        parts.append(f"<h3>Isolate {e(g.isolate_id)} - {e(g.organism or 'organism not supplied')} - "
                     f"{e(g.genome_accession or 'genome not supplied')}</h3>")
        if g.status != "COMPARED":
            parts.append("<p><strong>Not compared:</strong></p><ul>"
                         + "".join(f"<li>{e(x)}</li>" for x in g.reasons) + "</ul>")
        else:
            parts += [_tbl(g.run_summary), _tbl(g.presence), _tbl(pd.DataFrame(g.pairwise))]
        for w in g.warnings:
            parts.append(f"<p class='warn'>{e(w)}</p>")
        parts.append("<ul>" + "".join(f"<li>{e(n)}</li>" for n in g.notes) + "</ul>")
    parts.append("<h2>Limitations</h2><ul>" + "".join(f"<li>{e(x)}</li>" for x in LIMITATIONS) + "</ul>")
    parts.append("</body></html>")
    return "\n".join(parts)