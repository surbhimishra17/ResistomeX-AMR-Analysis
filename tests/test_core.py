"""Unit tests for resistomex_core.py.  Run from the project folder:
    python -m unittest tests.test_core -v
Fixtures amr_results.tsv / ResFinder_results_tab.txt are REAL outputs (E. coli K-12 MG1655).
Everything built in this file (SYNTHETIC_*) is made-up test data with fake gene names;
it only tests parser/logic behaviour and says nothing about real detections."""
import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
import resistomex_core as rc

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
def real(name):
    with open(os.path.join(FIX, name), "rb") as fh:
        return fh.read()

META = dict(isolate_id="Ecoli_K12_MG1655", organism="Escherichia coli", genome_accession="GCF_000005845.2")
def M(**kw):
    d = dict(META); d.update(kw); return rc.RunMeta(**d)

AMRF_HDR = ["Protein id","Contig id","Start","Stop","Strand","Element symbol","Element name","Scope","Type","Subtype","Class","Subclass","Method","Target length","Reference sequence length","% Coverage of reference","% Identity to reference","Alignment length","Closest reference accession","Closest reference name","HMM accession","HMM description","Hierarchy node"]
def amrf_row(sym, typ="AMR", sub="AMR", cls="BETA-LACTAM", method="EXACTX", contig="contigA"):
    r = dict.fromkeys(AMRF_HDR, "NA"); r.update({"Contig id": contig, "Start": "1", "Stop": "100", "Element symbol": sym,
        "Type": typ, "Subtype": sub, "Class": cls, "Method": method, "% Coverage of reference": "100.00", "% Identity to reference": "99.50"})
    return "\t".join(r[h] for h in AMRF_HDR)
def amrf(*rows): return ("\t".join(AMRF_HDR) + "\n" + "\n".join(rows) + "\n").encode()

RGI_HDR = ["ORF_ID","Contig","Start","Stop","Orientation","Cut_Off","Pass_Bitscore","Best_Hit_Bitscore","Best_Hit_ARO","Best_Identities","ARO","Model_type","SNPs_in_Best_Hit_ARO","Other_SNPs","Drug Class","Resistance Mechanism","AMR Gene Family","Predicted_DNA","Predicted_Protein","CARD_Protein_Sequence","Percentage Length of Reference Sequence","ID","Model_ID","Nudged","Note"]
def rgi_row(aro, cut="Strict", model="protein homolog model", snps="n/a", drug="classA; classB", fam="SYNTH family", nudged="FALSE", contig="contigA_1"):
    r = dict.fromkeys(RGI_HDR, "n/a"); r.update({"Contig": contig, "Cut_Off": cut, "Best_Hit_ARO": aro, "Best_Identities": "98.5",
        "Model_type": model, "SNPs_in_Best_Hit_ARO": snps, "Drug Class": drug, "AMR Gene Family": fam,
        "Percentage Length of Reference Sequence": "100.00", "Nudged": nudged, "Start": "5", "Stop": "50"})
    return "\t".join(r[h] for h in RGI_HDR)
def rgi(*rows): return ("\t".join(RGI_HDR) + "\n" + "\n".join(rows) + "\n").encode()

RF_HDR = "Resistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\tPosition in reference\tContig\tPosition in contig\tPhenotype\tAccession no.\n"
def resf(*rows): return (RF_HDR + "".join(r + "\n" for r in rows)).encode()
def resf_row(sym, contig="contigA"): return f"{sym}\t100.00\t861/861\t100.00\t1..861\t{contig}\t10..870\tSynthetic-drug\tXX000001"

class RealFiles(unittest.TestCase):
    def setUp(self):
        self.a = rc.parse_run(real("amr_results.tsv"), "amr_results.tsv", rc.TOOL_AMRFINDER, M(tool_version="x", database_version="y"))
        self.r = rc.parse_run(real("ResFinder_results_tab.txt"), "ResFinder_results_tab.txt", rc.TOOL_RESFINDER, M())
    def test_amrfinder_counts_and_categories(self):
        self.assertEqual(self.a.status, rc.STATUS_HITS)
        rec = self.a.records
        self.assertEqual(len(rec), 5)
        self.assertEqual((rec.category == "AMR").sum(), 3)
        self.assertEqual((rec.category == "STRESS").sum(), 2)
        self.assertEqual(sorted(rec[rec.category == "AMR"].gene_symbol), ["acrF", "blaEC", "mdtM"])
        self.assertTrue((rec.detection_kind == "gene").all())
    def test_metadata_preserved(self):
        rec = self.a.records
        self.assertTrue((rec.organism == "Escherichia coli").all())
        self.assertTrue((rec.genome_accession == "GCF_000005845.2").all())
        self.assertTrue((rec.tool == "AMRFinderPlus").all())
        self.assertEqual(rec.database_version.iloc[0], "y")
        self.assertEqual(rec.identity.iloc[0], 100.0)
    def test_resfinder_zero_hit_valid(self):
        self.assertEqual(self.r.status, rc.STATUS_ZERO)
        self.assertEqual(self.r.tool, rc.TOOL_RESFINDER)
        self.assertEqual(len(self.r.records), 0)
    def test_summary_includes_zero_hit_run(self):
        data = rc.combine_runs([self.a, self.r]); s = rc.summarize([self.a, self.r], data)
        self.assertEqual((s["valid_runs"], s["zero_hit_runs"], s["isolates"]), (2, 1, 1))
        self.assertEqual((s["amr_gene_detections"], s["unique_amr_genes"], s["stress"]), (3, 3, 2))
        rt = rc.run_table([self.a, self.r])
        self.assertEqual(set(rt.tool), {"AMRFinderPlus", "ResFinder"})  # zero-hit tool stays visible
    def test_comparison_real_pair(self):
        g = rc.compare_tools([self.a, self.r])[0]
        self.assertEqual(g.status, "COMPARED")
        self.assertEqual(len(g.presence), 3)
        self.assertTrue((g.presence.match_level == "Reported by one tool only").all())
        self.assertFalse(g.presence["ResFinder"].any())
        self.assertEqual(g.pairwise[0]["shared (exact)"], 0)
        self.assertTrue(any("zero AMR genes" in n for n in g.notes))
    def test_zero_hit_vs_hits_asks_to_confirm_same_assembly(self):
        g = rc.compare_tools([self.a, self.r])[0]
        self.assertTrue(any("same assembly" in w for w in g.warnings))
    def test_classes_not_invented(self):
        d = rc.harmonise_classes(rc.combine_runs([self.a]))
        got = dict(zip(d.gene_symbol, d.class_label))
        self.assertEqual(got["blaEC"], "BETA-LACTAM")
        self.assertEqual(got["acrF"], "EFFLUX")
        self.assertTrue(d[d.gene_symbol == "acrF"].class_note.iloc[0].startswith("EFFLUX is a mechanism"))
        self.assertIn("STRESS annotation", got["ariR"])      # not presented as antibiotic class
        self.assertEqual(len(d), 5)

class Validity(unittest.TestCase):
    def test_empty_file_is_not_zero_hit(self):
        for blob in (b"", b"  \n\n"):
            self.assertEqual(rc.parse_run(blob, "e.tsv", rc.AUTO, M()).status, rc.STATUS_EMPTY)
    def test_junk_header_only_rejected(self):
        self.assertEqual(rc.parse_run(b"foo\tbar\tbaz\n", "j.tsv", rc.AUTO, M()).status, rc.STATUS_SCHEMA)
    def test_tool_mismatch_and_auto(self):
        self.assertEqual(rc.parse_run(real("ResFinder_results_tab.txt"), "r.txt", rc.TOOL_AMRFINDER, M()).status, rc.STATUS_MISMATCH)
        self.assertEqual(rc.parse_run(real("ResFinder_results_tab.txt"), "r.txt", rc.AUTO, M()).status, rc.STATUS_ZERO)
    def test_missing_isolate(self):
        self.assertEqual(rc.parse_run(real("amr_results.tsv"), "a.tsv", rc.AUTO, M(isolate_id=" ")).status, rc.STATUS_META)
    def test_ragged_file_does_not_crash(self):
        bad = b"Resistance gene\tIdentity\tCoverage\nx\t1\t2\t3\t4\t5\n"
        run = rc.parse_run(bad, "b.txt", rc.AUTO, M())
        self.assertIn(run.status, {rc.STATUS_PARSE, rc.STATUS_HITS})
    def test_old_amrfinder_header_with_hash(self):
        hdr = "#Target identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\n"
        row = "x\tc1\t1\t9\t+\tfakeA\tfake protein\tcore\tAMR\tAMR\tFAKE\tFAKE\tBLASTX\n"
        run = rc.parse_run((hdr + row).encode(), "old.tsv", rc.AUTO, M())
        self.assertEqual(run.status, rc.STATUS_HITS); self.assertEqual(run.records.gene_symbol.iloc[0], "fakeA")
    def test_dropped_rows_reported(self):
        run = rc.parse_run(amrf(amrf_row("fakeA"), amrf_row("")), "d.tsv", rc.AUTO, M())
        self.assertEqual((len(run.records), run.n_rows_dropped), (1, 1))
        self.assertTrue(any("dropped" in i["message"] for i in rc.validate_runs([run])))
    def test_validation_messages(self):
        a = rc.parse_run(real("amr_results.tsv"), "a.tsv", rc.AUTO, M())
        msgs = " ".join(i["message"] for i in rc.validate_runs([a]))
        self.assertIn("not recorded", msgs)
        b = rc.parse_run(real("amr_results.tsv"), "b.tsv", rc.AUTO, M(organism="Klebsiella pneumoniae"))
        self.assertTrue(any(i["severity"] == "error" and "different organisms" in i["message"] for i in rc.validate_runs([a, b])))

class Parsers(unittest.TestCase):
    def test_amrfinder_point_mutation_kept_separate(self):
        run = rc.parse_run(amrf(amrf_row("gyrA_S83L", sub="POINT", cls="QUINOLONE", method="POINTX"), amrf_row("fakeA")), "p.tsv", rc.AUTO, M())
        pm = run.records[run.records.detection_kind == "point_mutation"].iloc[0]
        self.assertEqual((pm.gene_symbol, pm.mutation, pm.category), ("gyrA", "S83L", "AMR"))
        self.assertEqual(rc.summarize([run], rc.combine_runs([run]))["amr_gene_detections"], 1)
    def test_amrfinder_virulence_and_partial_flag(self):
        run = rc.parse_run(amrf(amrf_row("vfX", typ="VIRULENCE", sub="VIRULENCE", cls="NA"), amrf_row("fakeA", method="PARTIAL")), "v.tsv", rc.AUTO, M())
        self.assertEqual(rc.summarize([run], rc.combine_runs([run]))["virulence"], 1)
        self.assertIn("partial", run.records[run.records.gene_symbol == "fakeA"].review_flag.iloc[0])
    def test_rgi_native(self):
        run = rc.parse_run(rgi(rgi_row("SynthARO-1"), rgi_row("SynthARO-2", cut="Loose", nudged="TRUE"),
                               rgi_row("Synth gyrA variant", model="protein variant model", snps="S83L")), "r.txt", rc.AUTO, M())
        self.assertEqual(run.tool, rc.TOOL_RGI); rec = run.records
        self.assertEqual(rec.gene_symbol.iloc[0], "SynthARO-1")          # NOT 'AMR Gene Family'
        self.assertEqual(rec.coverage.iloc[0], 100.0); self.assertEqual(rec.identity.iloc[0], 98.5)
        self.assertTrue((rec.category == "AMR").all())
        self.assertEqual(rec.detection_kind.tolist(), ["gene", "gene", "point_mutation"])
        self.assertIn("Loose", rec.review_flag.iloc[1]); self.assertEqual(rec.mutation.iloc[2], "S83L")
        self.assertEqual(rc.harmonise_classes(rc.combine_runs([run])).class_label.iloc[0], "Multiple classes (tool-reported)")
    def test_pointfinder_synthetic(self):
        txt = b"Mutation\tNucleotide change\tAmino acid change\tResistance\tPMID\nfakeP p.A1B\tc.1A>B\tA1B\tSynthetic-drug\t0\n"
        run = rc.parse_run(txt, "pf.txt", rc.AUTO, M())
        self.assertEqual((run.tool, run.records.detection_kind.iloc[0]), (rc.TOOL_POINTFINDER, "point_mutation"))

class Comparison(unittest.TestCase):
    def three(self, **over):
        a = rc.parse_run(amrf(amrf_row("blaTEST-1"), amrf_row("onlyA")), "a.tsv", rc.AUTO, M())
        r = rc.parse_run(resf(resf_row("blaTEST-1"), resf_row("blaFAM-2B")), "r.txt", rc.AUTO, M(**over))
        g = rc.parse_run(rgi(rgi_row("blaTEST-1"), rgi_row("blaFAM-1")), "g.txt", rc.AUTO, M())
        return a, r, g
    def test_three_tools_and_family_hint(self):
        g = rc.compare_tools(list(self.three()))[0]
        self.assertEqual((g.status, len(g.tools)), ("COMPARED", 3))
        lv = dict(zip(g.presence.gene, g.presence.match_level))
        self.assertEqual(lv["blaTEST-1"], "Exact match (all tools)")
        self.assertEqual(lv["onlyA"], "Reported by one tool only")
        self.assertTrue(lv["blaFAM-2B"].startswith("Family-level")); self.assertTrue(lv["blaFAM-1"].startswith("Family-level"))
        self.assertEqual(len(g.pairwise), 3)
    def test_symbol_normalisation(self):
        a = rc.parse_run(amrf(amrf_row("tet(A)")), "a.tsv", rc.AUTO, M()); r = rc.parse_run(resf(resf_row("tetA")), "r.txt", rc.AUTO, M())
        self.assertEqual(rc.compare_tools([a, r])[0].presence.match_level.iloc[0], "Exact match (all tools)")
    def test_blocked_cases(self):
        for over, text in [(dict(genome_accession="GCF_000000001.1"), "different genome"), (dict(genome_accession=""), "Genome accession not supplied"),
                           (dict(organism="Klebsiella pneumoniae"), "different organisms"), (dict(organism=""), "Organism not supplied")]:
            a, r, _ = self.three(**over)
            g = rc.compare_tools([a, r])[0]
            self.assertEqual(g.status, "NOT_COMPARED", over); self.assertTrue(any(text in x for x in g.reasons), (over, g.reasons))
    def test_duplicate_tool_blocked(self):
        a = rc.parse_run(amrf(amrf_row("x1")), "a1.tsv", rc.AUTO, M()); a2 = rc.parse_run(amrf(amrf_row("x2")), "a2.tsv", rc.AUTO, M())
        r = rc.parse_run(resf(resf_row("x1")), "r.txt", rc.AUTO, M())
        self.assertEqual(rc.compare_tools([a, a2, r])[0].status, "NOT_COMPARED")
    def test_invalid_run_excluded_not_compared(self):
        a = rc.parse_run(amrf(amrf_row("x1")), "a.tsv", rc.AUTO, M()); bad = rc.parse_run(b"", "r.txt", rc.AUTO, M())
        g = rc.compare_tools([a, bad])[0]
        self.assertEqual(g.status, "NOT_COMPARED"); self.assertTrue(any("Excluded invalid" in w for w in g.warnings))
    def test_contig_mismatch_warning(self):
        a = rc.parse_run(amrf(amrf_row("x1", contig="c1")), "a.tsv", rc.AUTO, M()); r = rc.parse_run(resf(resf_row("x1", contig="zz9")), "r.txt", rc.AUTO, M())
        self.assertTrue(any("share no contig" in w for w in rc.compare_tools([a, r])[0].warnings))
    def test_rgi_orf_suffix_is_not_a_contig_mismatch(self):
        a = rc.parse_run(amrf(amrf_row("x1", contig="ctg7")), "a.tsv", rc.AUTO, M())
        g = rc.parse_run(rgi(rgi_row("x1", contig="ctg7_12")), "g.txt", rc.AUTO, M())
        self.assertFalse(any("share no contig" in w for w in rc.compare_tools([a, g])[0].warnings))
    def test_mixed_organism_project_kept_separate(self):
        runs = []
        for iso, org, acc in [("EC1", "E. coli", "GCF_000000011.1"), ("KP1", "K. pneumoniae", "GCF_000000022.1"), ("SE1", "Salmonella enterica subsp. enterica serovar Synth", "GCF_000000033.1")]:
            m = dict(isolate_id=iso, organism=org, genome_accession=acc)
            runs.append(rc.parse_run(amrf(amrf_row("g_" + iso)), iso + "_a.tsv", rc.AUTO, rc.RunMeta(**m)))
            runs.append(rc.parse_run(resf(resf_row("g_" + iso)), iso + "_r.txt", rc.AUTO, rc.RunMeta(**m)))
        groups = rc.compare_tools(runs)
        self.assertEqual([g.isolate_id for g in groups], ["EC1", "KP1", "SE1"]); self.assertTrue(all(g.status == "COMPARED" for g in groups))
        self.assertEqual(groups[0].organism, "Escherichia coli"); self.assertEqual(rc.summarize(runs, rc.combine_runs(runs))["organisms"], 3)
        self.assertTrue(all(len(g.presence) == 1 and g.presence.match_level.iloc[0].startswith("Exact") for g in groups))  # no cross-isolate mixing

class Mapping(unittest.TestCase):
    def test_conflicts_dropped_and_no_row_multiplication(self):
        m, w = rc.load_mapping(pd.DataFrame({"gene_symbol": ["sulX", "sulX", "okY"], "antibiotic_class": ["A", "B", "C"]}))
        self.assertEqual(m, {"OKY": "C"}); self.assertTrue(w)
        run = rc.parse_run(resf(resf_row("sulX"), resf_row("okY")), "r.txt", rc.AUTO, M())
        d = rc.harmonise_classes(rc.combine_runs([run]), m)
        self.assertEqual(len(d), 2); self.assertEqual(dict(zip(d.gene_symbol, d.class_source)), {"sulX": "none", "okY": "user mapping"})
    def test_bad_mapping_columns(self):
        with self.assertRaises(ValueError): rc.load_mapping(pd.DataFrame({"a": [1]}))

class Report(unittest.TestCase):
    def test_report_content_and_escaping(self):
        a = rc.parse_run(real("amr_results.tsv"), "<script>x</script>.tsv", rc.AUTO, M(tool_version="v", database_version="d"))
        r = rc.parse_run(real("ResFinder_results_tab.txt"), "r.txt", rc.AUTO, M())
        data = rc.harmonise_classes(rc.combine_runs([a, r]))
        html_ = rc.build_report_html([a, r], data, rc.compare_tools([a, r]), rc.validate_runs([a, r]))
        for needle in ("GCF_000005845.2", "Escherichia coli", "Valid zero-hit result", "STRESS annotations (not AMR)", "Limitations", "not recorded"):
            self.assertIn(needle, html_)
        self.assertNotIn("<script>x", html_); self.assertNotIn("14/14", html_)
    def test_report_with_only_zero_hit_runs(self):
        r = rc.parse_run(real("ResFinder_results_tab.txt"), "r.txt", rc.AUTO, M())
        html_ = rc.build_report_html([r], rc.combine_runs([r]), rc.compare_tools([r]), rc.validate_runs([r]))
        self.assertIn("Valid zero-hit", html_)

if __name__ == "__main__":
    unittest.main()