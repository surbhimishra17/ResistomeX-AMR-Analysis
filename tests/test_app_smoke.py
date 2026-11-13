"""Smoke test of app.py with a STUBBED streamlit/plotly (so it runs without either installed).
It proves app.py executes end-to-end on real/demo inputs without Python errors and that chart
calls reference existing columns. It does NOT test real Streamlit rendering or real Plotly output."""
import os, sys, runpy, types, unittest
from unittest.mock import MagicMock
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import pandas as pd
FIX = os.path.join(ROOT, "tests", "fixtures")

class Up:
    def __init__(self, name, blob): self.name, self._b = name, blob
    def getvalue(self): return self._b
class Stopped(Exception): pass

def run_app(uploads=None, demo=False, edit=None):
    st = MagicMock(); log = {"metrics": {}, "downloads": {}, "errors": [], "charts": [], "dfs": 0}
    st.toggle.return_value = demo
    st.file_uploader.side_effect = [uploads or [], None]
    def mkcol():
        c = MagicMock()
        c.metric.side_effect = lambda label, value, *a, **k: log["metrics"].__setitem__(label, value)
        c.download_button.side_effect = lambda label, data, *a, **k: log["downloads"].__setitem__(label, data)
        return c
    st.columns.side_effect = lambda n: [mkcol() for _ in range(n if isinstance(n, int) else len(n))]
    st.tabs.side_effect = lambda labels: [MagicMock() for _ in labels]
    st.stop.side_effect = Stopped
    def editor(df, **kw):
        df = df.copy()
        if edit: edit(df)
        return df
    st.data_editor.side_effect = editor
    st.error.side_effect = lambda m, *a, **k: log["errors"].append(m)
    st.dataframe.side_effect = lambda *a, **k: log.__setitem__("dfs", log["dfs"] + 1)
    st.download_button.side_effect = lambda label, data, *a, **k: log["downloads"].__setitem__(label, data)
    px = types.ModuleType("plotly.express")
    def mk(kind):
        def f(df, **kw):
            for k in ("x", "y", "color"):
                if k in kw: assert kw[k] in df.columns, (kind, k, kw[k], list(df.columns))
            log["charts"].append(kind); return MagicMock()
        return f
    px.scatter, px.bar = mk("scatter"), mk("bar")
    plotly = types.ModuleType("plotly"); plotly.express = px
    saved = {k: sys.modules.get(k) for k in ("streamlit", "plotly", "plotly.express")}
    sys.modules.update({"streamlit": st, "plotly": plotly, "plotly.express": px})
    try:
        try: runpy.run_path(os.path.join(ROOT, "app.py"), run_name="__main__")
        except Stopped: log["stopped"] = True
    finally:
        for k, v in saved.items():
            sys.modules.pop(k, None)
            if v is not None: sys.modules[k] = v
    return log

def real(n):
    with open(os.path.join(FIX, n), "rb") as fh: return fh.read()

class AppSmoke(unittest.TestCase):
    def test_no_input_stops_cleanly(self):
        self.assertTrue(run_app().get("stopped"))
    def test_demo_mode(self):
        log = run_app(demo=True)
        self.assertEqual(log["metrics"]["Isolates"], 2); self.assertEqual(log["metrics"]["Organisms"], 2)
        self.assertEqual(log["metrics"]["Valid zero-hit runs"], 1); self.assertEqual(log["metrics"]["Invalid runs"], 0)
        self.assertEqual(log["charts"], ["scatter", "bar", "bar"]); self.assertIn("HTML report", log["downloads"])
    def test_real_files_with_metadata(self):
        def edit(df):
            df["Organism"], df["Genome accession"] = "Escherichia coli", "GCF_000005845.2"
            df["Isolate ID"] = "Ecoli_K12_MG1655"
        log = run_app([Up("amr_results.tsv", real("amr_results.tsv")), Up("ResFinder_results_tab.txt", real("ResFinder_results_tab.txt"))], edit=edit)
        m = log["metrics"]
        self.assertEqual((m["Isolates"], m["Valid runs"], m["Valid zero-hit runs"], m["Invalid runs"]), (1, 2, 1, 0))
        self.assertEqual((m["AMR gene detections"], m["Unique AMR genes"], m["STRESS annotations"], m["VIRULENCE annotations"]), (3, 3, 2, 0))
        rep = log["downloads"]["HTML report"].decode()
        for s in ("GCF_000005845.2", "Escherichia coli", "Valid zero-hit result"): self.assertIn(s, rep)
        self.assertEqual(log["errors"], [])
    def test_zero_hit_only_still_gives_report(self):
        log = run_app([Up("r.txt", real("ResFinder_results_tab.txt"))])
        self.assertEqual(log["metrics"]["Valid zero-hit runs"], 1); self.assertEqual(log["charts"], [])
        self.assertIn("HTML report", log["downloads"])
    def test_invalid_files_reported_not_crashing(self):
        log = run_app([Up("e.tsv", b""), Up("j.tsv", b"foo\tbar\n"), Up("amr.tsv", real("amr_results.tsv"))])
        self.assertEqual(log["metrics"]["Invalid runs"], 2); self.assertTrue(len(log["errors"]) >= 2)

if __name__ == "__main__":
    unittest.main()