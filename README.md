# ResistomeX – AMR Gene Profiling, Harmonisation and Comparative Analysis

**Working prototype (v0.1).** A tool for organising and comparing the output files of AMR (antimicrobial resistance) detection tools.

> **Status: working prototype (v0.1).** Tested on one real genome and on small synthetic files. Not tested at
> large scale. For research and learning use only; **not for clinical use**. A detected gene is not proof of
> clinical resistance.

## What it is

ResistomeX is a small Streamlit app. You run AMRFinderPlus, ResFinder or CARD/RGI yourself and upload their output
files. The app puts the results into one standard table together with sample information (isolate, organism, genome
accession, tool, tool/database version), shows summaries and charts, and compares tools for the same isolate when the
organism and genome accession match. It is not a new detection method.

## What it currently does

- Reads AMRFinderPlus, ResFinder (acquired genes), PointFinder and CARD/RGI output, each with its own parser.
- Checks that the file header matches the tool. A file with the correct header and no rows is kept as a **valid
  zero-hit result**; empty, mismatched or unreadable files are marked **invalid** and left out.
- Keeps AMR, STRESS and VIRULENCE annotations separate (only AMR is counted as AMR), and keeps genes separate from
  point mutations / variant models.
- Keeps isolate ID, organism, genome accession, tool, and tool/database version (entered by the user) on each record.
- Compares tools only within the same isolate, and only if organism and genome accession are identical.
- Shows a dashboard, charts coloured by class, tables, CSV downloads and an HTML report.
- Is designed so that several isolates and organisms can be loaded together (see testing status below).

## What it does not do

- It does not run the detection tools and does not detect or predict resistance.
- It does not compare different isolates or organisms as if they were one sample.
- It does not create antibiotic-class mappings. Classes are shown as each tool reports them; an optional
  user-supplied mapping CSV is used only where a tool gave no class.
- It does not compare point mutations across tools.
- It cannot check that two files came from the same sequence; it relies on the organism and genome accession
  you enter.

## How each tool's output is read

| Tool / file | How it is recognised | How it is treated |
|---|---|---|
| AMRFinderPlus TSV | `Element symbol` (or `Gene symbol`), `Type`, `Method`/`Class` | `Type` gives AMR / STRESS / VIRULENCE. POINT subtype or method gives a point mutation. Class kept as reported. |
| ResFinder `ResFinder_results_tab.txt` | `Resistance gene`, `Identity`, `Coverage` | Acquired AMR genes. `Phenotype` kept as text, not as a class. |
| PointFinder results | `Mutation`, `Resistance` | Point mutations. |
| CARD/RGI main `.txt` | `Best_Hit_ARO`, `Cut_Off`, `Model_type` | Gene = `Best_Hit_ARO`. Homolog model = gene, variant model = mutation. Loose hits are flagged. |

## How comparison works

1. Runs are grouped by isolate ID.
2. All runs in a group need the same organism (species) and the same genome accession, at least two different tools,
   and one file per tool. Otherwise the group is shown as "not compared" with the reason.
3. AMR genes are matched by normalised symbol (`tet(A)` = `tetA`). Name or allele differences are shown as
   "family-level only, verify manually" and are not counted as matches.
4. Tools have different scopes (for example ResFinder acquired genes vs AMRFinderPlus), so a gene found by one tool
   and not another is not automatically a disagreement.
5. If one tool reports zero hits and another reports hits, the app shows a caution to confirm that all files came from
   the same assembly (zero-hit files have no contig IDs, so the app cannot check this itself).

## Testing status

| What | Status |
|---|---|
| AMRFinderPlus and ResFinder output for *E. coli* K-12 MG1655 (GCF_000005845.2): 5 AMRFinderPlus rows (3 AMR, 2 STRESS) and a header-only ResFinder file | Tested on real files |
| 38 automated tests (parsers, validation, comparison rules, report, app smoke test) | Pass. Most use small synthetic files with made-up gene names. |
| App interface | Run manually on a small set of files. The automated app test uses a stub, not real Streamlit. |
| Several organisms (*Klebsiella*, *Salmonella*) | Synthetic data only; no real output tested yet |
| CARD/RGI and PointFinder parsers | Written from documented columns; not tested on real output |
| Large files, many isolates, speed | Not tested |
| Other tool versions or output layouts | Not tested (developed with Python 3.14, Streamlit 1.67) |

## Run it

```bash
pip install -r requirements.txt
python -m streamlit run app.py
```

Tests (need only pandas):

```bash
python -m unittest tests.test_core tests.test_app_smoke
```

A fictional demo mode (made-up gene names) is available in the sidebar.

## Known limitations

- Matching across tools depends on gene naming; CARD names often differ from AMRFinderPlus/ResFinder symbols.
- Tool and database versions are not in the output files, so they must be typed in.
- Results depend on each tool's database, settings and thresholds. Zero hits means no hits under those settings,
  not proof that no resistance determinants exist.
- Files are re-read on every interaction, so large inputs may be slow.

## Possible next steps

- Test with real CARD/RGI and PointFinder output.
- Test with real *Klebsiella pneumoniae* and *Salmonella enterica* genomes.
- Test with larger sets of isolates and files.

## Background

ResistomeX started as a team prototype at the Bioverse 2026 hackathon (3rd prize), where it was presented as a working prototype using AMRFinderPlus and ResFinder results. After the event the output files were examined more closely: the ResFinder acquired-gene file had a header but no result rows (zero hits), so the AMRFinderPlus-ResFinder agreement shown in the presentation was not supported by the files. AMRFinderPlus also reports STRESS and VIRULENCE elements, which must not be counted as AMR. This version was rebuilt to handle these points.

## Team Credits

- Surbhi Mishra
- Rashi Jain
- Kajal Agarwal
- Tarun Deb Sharma

## Citing the underlying tools

ResistomeX only organises output from AMRFinderPlus (NCBI), ResFinder/PointFinder (Center for Genomic Epidemiology)
and CARD/RGI. Please cite those tools and record their database versions when using results.