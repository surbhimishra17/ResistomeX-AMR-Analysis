import hashlib
import os

import pandas as pd
import plotly.express as px
import streamlit as st

import resistomex_core as rc

st.set_page_config(page_title="ResistomeX | AMR Gene Profiling, Harmonisation and Comparative Analysis", page_icon="🧬", layout="wide")

st.title("🧬 ResistomeX")
st.caption("AMR Gene Profiling, Harmonisation and Comparative Analysis • working prototype (v0.1)")
st.info(
    "ResistomeX organises results produced by AMRFinderPlus, ResFinder and CARD/RGI. It does not "
    "detect resistance, predict phenotype or recommend treatment. A detected gene is not proof of "
    "clinical resistance."
)

CATEGORY_COLORS = {"AMR": "#D55E00", "STRESS": "#E69F00", "VIRULENCE": "#CC79A7", "OTHER": "#8C8C8C"}
_PALETTE = ["#0072B2", "#009E73", "#56B4E9", "#F0E442", "#CC79A7", "#E69F00", "#7F3C8D", "#11A579",
            "#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695", "#CF1C90", "#4B4B8F"]
_FIXED = {
    "Not reported by tool": "#9AA0A6",
    "Multiple classes (tool-reported)": "#4B3B7A",
    "Reported differently by tools": "#6D4C41",
}
SHOW_COLS = ["isolate_id", "organism", "genome_accession", "tool", "tool_version", "database_version",
             "category", "detection_kind", "gene_symbol", "mutation", "class_label", "class_source",
             "reported_class", "reported_phenotype", "identity", "coverage", "evidence_level",
             "review_flag", "source_file"]


def class_colors(labels):
    cmap, i = dict(_FIXED), 0
    for lab in sorted(set(labels)):
        if lab not in cmap:
            cmap[lab] = _PALETTE[i % len(_PALETTE)]
            i += 1
    return cmap


def filename_id(name):
    base = os.path.basename(name)
    for suf in (".tsv", ".csv", ".txt"):
        if base.lower().endswith(suf):
            base = base[: -len(suf)]
    return base or "isolate"


def demo_files():
    """FICTIONAL demo data in native-like layouts (fake gene names). Goes through the real parsers."""
    af = "Contig id\tElement symbol\tType\tSubtype\tClass\tMethod\t% Coverage of reference\t% Identity to reference\n"
    rf = "Resistance gene\tIdentity\tCoverage\tContig\tPhenotype\n"
    rg = "Contig\tCut_Off\tBest_Hit_ARO\tBest_Identities\tModel_type\tDrug Class\tPercentage Length of Reference Sequence\n"
    ec = {"id": "DEMO-EC-01", "org": "Escherichia coli", "acc": "DEMO_ASSEMBLY_EC01"}
    kp = {"id": "DEMO-KP-01", "org": "Klebsiella pneumoniae", "acc": "DEMO_ASSEMBLY_KP01"}
    return [
        ("demo_EC01_amrfinder.tsv", rc.TOOL_AMRFINDER, ec, af + "ctg1\tblaDEMO-1\tAMR\tAMR\tBETA-LACTAM\tEXACTX\t100\t100\n"
         "ctg1\tsulDEMO\tAMR\tAMR\tSULFONAMIDE\tEXACTX\t100\t99.5\nctg2\tstressDEMO\tSTRESS\tBIOCIDE\tNA\tEXACTX\t100\t100\n"),
        ("demo_EC01_resfinder.txt", rc.TOOL_RESFINDER, ec, rf + "blaDEMO-1\t100\t100\tctg1\tDemo-drug-A\nsulDEMO\t99.5\t100\tctg1\tDemo-drug-B\n"),
        ("demo_EC01_rgi.txt", rc.TOOL_RGI, ec, rg + "ctg1_5\tStrict\tblaDEMO-2\t98.0\tprotein homolog model\tclass A; class B\t100\n"),
        ("demo_KP01_amrfinder.tsv", rc.TOOL_AMRFINDER, kp, af + "ctgA\ttetDEMO(A)\tAMR\tAMR\tTETRACYCLINE\tEXACTX\t100\t100\n"),
        ("demo_KP01_resfinder.txt", rc.TOOL_RESFINDER, kp, rf),
    ]


# ------------------------------------------------------------------ sidebar
with st.sidebar:
    st.header("1. Input data")
    demo_mode = st.toggle("Load FICTIONAL demo data", value=False)
    st.caption("Demo records use made-up gene names. Turn off for real results.")
    uploads = st.file_uploader("Upload native tool output files", type=["csv", "tsv", "txt"],
                               accept_multiple_files=True)
    st.caption("Supported: AMRFinderPlus TSV, ResFinder acquired-gene tab file, ResFinder "
               "PointFinder table, CARD/RGI main .txt. Each tool is read with its own parser.")
    st.divider()
    st.header("2. Gene-class mapping (optional)")
    map_upload = st.file_uploader("Mapping CSV", type=["csv"])
    st.caption("Columns: gene_symbol, antibiotic_class. Used only where a tool reported no class. "
               "No mapping is applied by default.")

mapping, map_warnings = {}, []
if map_upload is not None:
    try:
        mapping, map_warnings = rc.load_mapping(pd.read_csv(map_upload, dtype=str))
    except Exception as exc:
        st.sidebar.error(f"Mapping not used: {exc}")
for w in map_warnings:
    st.sidebar.warning(w)

# ------------------------------------------------------------------ build runs
runs = []
if uploads:
    st.subheader("Input files and sample metadata")
    st.caption("Use the SAME Isolate ID for files from the same isolate. Organism and genome accession are "
               "required for tool-to-tool comparison. Tool/database versions are not stored in the output "
               "files, so enter them if you want them recorded.")
    blobs = [(u.name, u.getvalue()) for u in uploads]
    key = "meta_" + hashlib.md5("|".join(f"{n}:{len(b)}" for n, b in blobs).encode()).hexdigest()
    base = pd.DataFrame({
        "File": [n for n, _ in blobs], "Tool": [rc.AUTO] * len(blobs),
        "Isolate ID": [filename_id(n) for n, _ in blobs], "Organism": [""] * len(blobs),
        "Genome accession": [""] * len(blobs), "Tool version": [""] * len(blobs),
        "Database version": [""] * len(blobs),
    })
    edited = st.data_editor(
        base, key=key, hide_index=True, disabled=["File"], num_rows="fixed",
        column_config={
            "Tool": st.column_config.SelectboxColumn("Tool", options=[rc.AUTO] + rc.TOOLS, required=True),
            "Organism": st.column_config.TextColumn("Organism", help="e.g. Escherichia coli, Klebsiella pneumoniae"),
            "Genome accession": st.column_config.TextColumn("Genome accession", help="e.g. GCF_000005845.2"),
        },
    ).fillna("")
    for (name, blob), (_, row) in zip(blobs, edited.iterrows()):
        meta = rc.RunMeta(row["Isolate ID"], row["Organism"], row["Genome accession"],
                          row["Tool version"], row["Database version"])
        runs.append(rc.parse_run(blob, name, row["Tool"], meta))
elif demo_mode:
    st.warning("DEMO MODE - all data below is FICTIONAL (made-up gene names) and only shows how the app behaves.")
    for name, tool, m, text in demo_files():
        runs.append(rc.parse_run(text.encode(), name, tool, rc.RunMeta(m["id"], m["org"], m["acc"],
                                                                       "demo", "demo")))
else:
    st.warning("Upload one or more result files, or enable the fictional demo data.")
    st.stop()

issues = rc.validate_runs(runs)
data = rc.harmonise_classes(rc.combine_runs(runs), mapping)
summary = rc.summarize(runs, data)
comparisons = rc.compare_tools(runs)
data["sample"] = data["isolate_id"] + " | " + data["organism"].replace("", "organism not supplied")

# ------------------------------------------------------------------ run status
st.subheader("Input run status")
st.dataframe(rc.run_table(runs), hide_index=True)
st.caption("A zero-hit result is accepted only when the header matches the tool and there are no rows. "
           "Empty, mismatched or unreadable files are marked invalid and excluded.")
errs = [i for i in issues if i["severity"] == "error"]
for i in errs:
    st.error(f"{i['file']}: {i['message']}" if i["file"] else i["message"])
others = [i for i in issues if i["severity"] != "error"]
if others:
    with st.expander(f"Validation notes ({len(others)})"):
        st.dataframe(pd.DataFrame(others), hide_index=True)

# ------------------------------------------------------------------ overview
st.subheader("Overview")
r1 = st.columns(4)
r1[0].metric("Files", summary["files"])
r1[1].metric("Valid runs", summary["valid_runs"])
r1[2].metric("Valid zero-hit runs", summary["zero_hit_runs"])
r1[3].metric("Invalid runs", summary["invalid_runs"])
r2 = st.columns(3)
r2[0].metric("Isolates", summary["isolates"])
r2[1].metric("Organisms", summary["organisms"])
r2[2].metric("Genome accessions", summary["genomes"])
r3 = st.columns(5)
r3[0].metric("AMR gene detections", summary["amr_gene_detections"])
r3[1].metric("Unique AMR genes", summary["unique_amr_genes"])
r3[2].metric("AMR mutations/variants", summary["amr_mutations_variants"])
r3[3].metric("STRESS annotations", summary["stress"])
r3[4].metric("VIRULENCE annotations", summary["virulence"])
st.caption("STRESS and VIRULENCE annotations are NOT counted as AMR. Unique genes use normalised symbols.")

amr = data[data["category"] == "AMR"]
amr_genes = amr[amr["detection_kind"] == "gene"]
amr_mut = amr[amr["detection_kind"] != "gene"]
non_amr = data[data["category"] != "AMR"]

# ------------------------------------------------------------------ charts
st.subheader("AMR gene profile")
if amr_genes.empty:
    st.info("No AMR gene detections in the valid runs. Zero-hit runs are listed in the status table above; "
            "isolates with only zero-hit runs have no bars or dots.")
else:
    g = amr_genes.groupby(["sample", "gene_key"], as_index=False).agg(
        gene=("gene_symbol", "first"), tools=("tool", lambda s: ", ".join(sorted(set(s)))),
        n_tools=("tool", "nunique"), labels=("class_label", lambda s: "\x1f".join(sorted(set(s)))))
    g["class_label"] = g["labels"].map(lambda ls: ls if "\x1f" not in ls else "Reported differently by tools")
    cmap = class_colors(g["class_label"])
    fig = px.scatter(g, x="sample", y="gene", color="class_label", color_discrete_map=cmap,
                     hover_data={"tools": True, "n_tools": True, "gene_key": False, "labels": False},
                     title="AMR genes per isolate (colour = class as reported by the tool)")
    fig.update_traces(marker=dict(size=16, symbol="square", line=dict(width=1, color="white")))
    fig.update_layout(height=max(380, 34 * g["gene"].nunique() + 180), legend_title_text="Class (as reported)",
                      xaxis_title="Isolate | organism", yaxis_title="AMR gene", xaxis=dict(type="category"))
    st.plotly_chart(fig)
    cA, cB = st.columns(2)
    by_class = g.groupby(["sample", "class_label"], as_index=False).size().rename(columns={"size": "genes"})
    figb = px.bar(by_class, x="sample", y="genes", color="class_label", color_discrete_map=cmap,
                  title="Distinct AMR genes per isolate, by class")
    figb.update_layout(legend_title_text="Class (as reported)", xaxis_title="Isolate | organism",
                       yaxis_title="Distinct AMR genes", xaxis=dict(type="category"))
    cA.plotly_chart(figb)
    if not data.empty:
        cat = data.groupby(["sample", "category"], as_index=False).size().rename(columns={"size": "records"})
        figc = px.bar(cat, x="sample", y="records", color="category", color_discrete_map=CATEGORY_COLORS,
                      title="Annotation categories per isolate (AMR vs STRESS vs VIRULENCE)")
        figc.update_layout(legend_title_text="Category", xaxis_title="Isolate | organism",
                           yaxis_title="Records", xaxis=dict(type="category"))
        cB.plotly_chart(figc)
    st.caption("Classes are shown exactly as each tool reports them (or from your mapping). Different tools use "
               "different class vocabularies and EFFLUX is a mechanism-type term, so they are not merged.")

# ------------------------------------------------------------------ tables
st.subheader("Results by annotation type")
t1, t2, t3, t4 = st.tabs([f"AMR genes ({len(amr_genes)})", f"AMR mutations / variants ({len(amr_mut)})",
                          f"STRESS / VIRULENCE / other ({len(non_amr)})", "Gene-class summary"])
cols = [c for c in SHOW_COLS if c in data.columns]
with t1:
    if len(amr_genes):
        st.dataframe(amr_genes[cols], hide_index=True)
    else:
        st.info("None.")
with t2:
    st.caption("Point mutations and variant models are listed per tool and are not compared across tools.")
    if len(amr_mut):
        st.dataframe(amr_mut[cols], hide_index=True)
    else:
        st.info("None.")
with t3:
    st.caption("Not counted as AMR genes.")
    if len(non_amr):
        st.dataframe(non_amr[cols], hide_index=True)
    else:
        st.info("None.")
with t4:
    if amr_genes.empty:
        st.info("No AMR genes to summarise.")
    else:
        cs = amr_genes.drop_duplicates(["sample", "gene_key", "tool"]).pivot_table(
            index="class_label", columns="sample", values="gene_key", aggfunc="nunique", fill_value=0)
        st.dataframe(cs)
        unmapped = amr_genes[amr_genes["class_source"] == "none"]
        if len(unmapped):
            st.warning(f"{len(unmapped)} AMR record(s) have no class from the tool or your mapping. "
                       "They are kept as 'Not reported by tool' (nothing is guessed).")

# ------------------------------------------------------------------ comparison
st.subheader("Tool-to-tool comparison (same isolate, same organism, same genome only)")
for gcmp in comparisons:
    st.markdown(f"**Isolate {gcmp.isolate_id}** • {gcmp.organism or 'organism not supplied'} • "
                f"{gcmp.genome_accession or 'genome accession not supplied'}")
    if gcmp.status != "COMPARED":
        st.warning("Not compared: " + " ".join(gcmp.reasons))
    else:
        st.dataframe(gcmp.run_summary, hide_index=True)
        pres = gcmp.presence.copy()
        for t in gcmp.tools:
            pres[t] = pres[t].map({True: "✓", False: "–"})
        if len(pres):
            st.dataframe(pres, hide_index=True)
        else:
            st.info("No AMR genes reported by any compared tool.")
        st.dataframe(pd.DataFrame(gcmp.pairwise), hide_index=True)
    for w in gcmp.warnings:
        st.warning(w)
    with st.expander("Scope notes"):
        for n in gcmp.notes:
            st.write("• " + n)
    st.divider()
st.caption("Genome compatibility is judged only from the organism and genome accession you entered; the app cannot verify the sequences. "
           "Exact matches use normalised symbols. 'Family-level only' means the allele/name differs and must be "
           "checked manually. A gene missing from a narrower-scope tool is a scope difference, not proof of disagreement.")

# ------------------------------------------------------------------ downloads
st.subheader("Download outputs")
d1, d2, d3, d4 = st.columns(4)
d1.download_button("Standardised results CSV", data.drop(columns=["sample"]).to_csv(index=False).encode("utf-8"),
                   "resistomex_standardised_results.csv", "text/csv")
d2.download_button("AMR-only CSV", amr.drop(columns=["sample"]).to_csv(index=False).encode("utf-8"),
                   "resistomex_amr_only.csv", "text/csv")
d3.download_button("Run status CSV", rc.run_table(runs).to_csv(index=False).encode("utf-8"),
                   "resistomex_run_status.csv", "text/csv")
d4.download_button("HTML report", rc.build_report_html(runs, data, comparisons, issues).encode("utf-8"),
                   "resistomex_report.html", "text/html")

with st.expander("What ResistomeX does, does not do, and its limitations"):
    st.markdown("**Does:** reads native tool tables with tool-specific parsers; validates headers; keeps zero-hit "
                "results; separates AMR / STRESS / VIRULENCE and genes / mutations; keeps organism, genome "
                "accession, tool and version metadata; compares tools within one isolate and genome.")
    st.markdown("**Does not:** detect resistance, predict phenotype, compare different isolates as if they were one, "
                "or invent gene-class mappings.")
    for item in rc.LIMITATIONS:
        st.write("• " + item)