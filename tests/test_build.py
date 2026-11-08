"""Tests: unit tests for matching/parsing, determinism of an offline rebuild,
and the around() query against the built data.

Run:  python3 -I -m unittest discover -s tests -v
The determinism and data tests need a populated cache/ and data/ (run build.py
first); they are skipped if those are missing.
"""
import filecmp
import json
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from tol import build as B  # noqa: E402
from tol.cache import Cache  # noqa: E402
from tol.query import Dataset  # noqa: E402

DATA = os.path.join(ROOT, "data")
CACHE = os.path.join(ROOT, "cache")
HAVE_DATA = os.path.exists(os.path.join(DATA, "clades.jsonl"))
HAVE_CACHE = os.path.isdir(CACHE) and bool(os.listdir(CACHE))


class TestUnits(unittest.TestCase):
    def test_rank_check(self):
        self.assertEqual(B.rank_check("phylum", "phylum"), "exact")
        self.assertEqual(B.rank_check("superorder", "order"), "near")
        self.assertEqual(B.rank_check("no rank", "unranked clade"), "unranked")
        self.assertIsNone(B.rank_check("phylum", "infraorder"))  # the sloth trap
        self.assertIsNone(B.rank_check("genus", "class"))

    def test_scoped_index_rejects_homonym(self):
        # A scope built from an anchor never contains the sloth 'Tardigrada'.
        idx = B.PbdbIndex()
        idx.add({"taxon_no": "67137", "orig_no": "67137", "accepted_no": "67137", "taxon_name": "Tardigrada",
                 "accepted_name": "Tardigrada", "accepted_rank": "phylum"})
        self.assertEqual(idx.candidates("Tardigrada"), [67137])
        # a synonym row is indexed under its accepted name only
        idx.add({"taxon_no": "182871", "orig_no": "57074", "accepted_no": "337045", "taxon_name": "Tardigrada",
                 "accepted_name": "Folivora", "accepted_rank": "unranked clade", "difference": "replaced by"})
        self.assertEqual(idx.candidates("Tardigrada"), [67137])
        self.assertEqual(idx.resolve(182871), 337045)

    def test_newick(self):
        p = B.parse_newick("((ott1,ott2)mrcaott1ott2,ott3)ott4;")
        self.assertEqual(p["ott1"], "mrcaott1ott2")
        self.assertEqual(p["mrcaott1ott2"], "ott4")
        self.assertEqual(p["ott3"], "ott4")
        self.assertIsNone(p["ott4"])

    def test_env_tags(self):
        self.assertEqual(B.env_tags("terrestrial"), ["terrestrial"])
        self.assertEqual(B.env_tags("marine"), ["marine"])
        self.assertEqual(B.env_tags(None), [])

    def test_dumps_is_canonical(self):
        self.assertEqual(B.dumps({"b": 1, "a": [2, 1]}), '{"a":[2,1],"b":1}')


@unittest.skipUnless(HAVE_CACHE and HAVE_DATA, "needs cache/ and data/ (run build.py)")
class TestDeterminism(unittest.TestCase):
    def test_offline_rebuild_is_byte_identical(self):
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            for out in (a, b):
                B.build(Cache(CACHE, offline=True, verbose=False), ROOT, out, log=lambda *x: None)
            names = sorted(os.listdir(a))
            self.assertEqual(names, sorted(os.listdir(b)))
            match, mismatch, errors = filecmp.cmpfiles(a, b, names, shallow=False)
            self.assertEqual(mismatch + errors, [])
            # and identical to the committed data/
            match, mismatch, errors = filecmp.cmpfiles(a, DATA, names, shallow=False)
            self.assertEqual(mismatch + errors, [], "data/ is stale: rebuild and commit")

    def test_manifest_hashes(self):
        m = json.load(open(os.path.join(DATA, "manifest.json"), encoding="utf-8"))
        for fn, info in m["files"].items():
            self.assertEqual(B.sha256_file(os.path.join(DATA, fn)), info["sha256"], fn)


@unittest.skipUnless(HAVE_DATA, "needs data/ (run build.py)")
class TestData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ds = Dataset(DATA)

    def test_ids_unique_and_parents_present(self):
        ids = [r["id"] for r in self.ds.clades]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(ids, sorted(ids, key=lambda i: int(i[4:])))  # numeric OTT order
        for r in self.ds.clades:
            if r["parent"] is not None:
                self.assertIn(r["parent"], self.ds.by_id, r["id"])
        roots = [r for r in self.ds.clades if r["parent"] is None]
        self.assertEqual([r["name"] for r in roots], ["Eukaryota"])

    def test_tardigrada(self):
        t = self.ds.by_id["ott:111438"]
        self.assertEqual(t["name"], "Tardigrada")
        self.assertEqual(t["pbdb"]["id"], "txn:67137")  # the phylum, not the sloths
        self.assertEqual(t["rank"], "phylum")
        self.assertIn("ott:816442", t["parent"])  # Panarthropoda
        self.assertGreater(t["fossil"]["first_max_ma"], 480)

    def test_no_sloths(self):
        for r in self.ds.clades + self.ds.genera:
            self.assertNotEqual((r.get("higher") or {}).get("class"), "Mammalia", r["id"])
            if r.get("pbdb"):
                self.assertNotIn(r["pbdb"]["name"], ("Folivora", "Eufolivora"))

    def test_divergence_has_citation(self):
        refs = json.load(open(os.path.join(DATA, "references.json"), encoding="utf-8"))
        n = 0
        for r in self.ds.clades:
            for d in r["divergence"]:
                self.assertIn(d["src"], refs)
                n += 1
        self.assertGreater(n, 5)

    def test_query_tardigrade_at_500ma(self):
        res = self.ds.around(500, size_mm=0.5, habitat="marine", kinds=("clade",))
        names = [r["name"] for r in res]
        self.assertIn("Tardigrada", names)
        for r in res:
            self.assertTrue(r["time_range"]["start_ma"] >= 500 >= r["time_range"]["end_ma"])
            self.assertIn("marine", r["habitat"])

    def test_query_time_filter(self):
        # Nothing in the package lived 3000 Ma (eukaryotes are younger).
        self.assertEqual(self.ds.around(3000, include_unknown_size=True, include_unknown_habitat=True), [])

    def test_query_size_filter(self):
        res = self.ds.around(0, size_mm=1000, decades=0.5, kinds=("clade",))
        self.assertNotIn("Tardigrada", [r["name"] for r in res])

    def test_slices_reference_existing_genera(self):
        for s in self.ds.slices:
            for ids in list(s["representatives"].values()) + list(s["by_size_decade"].values()):
                for i in ids:
                    self.assertIn(i, self.ds.by_id)


if __name__ == "__main__":
    unittest.main()
