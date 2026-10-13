import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from export_pea import SCHEMA, build_stage, markdown, taxon_from_leaf, write_report


class ExportPeaTest(unittest.TestCase):
    def test_taxon_parser(self):
        leaf = "Sapindales_Rutaceae_Citrus_sinensis_KEW_ERR"
        self.assertEqual(taxon_from_leaf(leaf, "family"), "Rutaceae")
        self.assertEqual(taxon_from_leaf(leaf, "species"), "Citrus sinensis")
        self.assertIsNone(taxon_from_leaf("NODE_1", "genus"))

    def test_stage_and_readable_report(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            matrix = root / "matrix"
            matrix.mkdir()
            (matrix / "4691.1.cleaned.csv").write_text(
                ",sample_4691_exon_1,Sapindales_Rutaceae_Citrus_sinensis\n"
                "Sapindales_Rutaceae_Citrus_sinensis,0.01,0\n"
                "Sapindales_Sapindaceae_Litchi_chinensis,0.3,0.4\n",
                encoding="utf-8",
            )
            (root / "predictions.csv").write_text(
                "row_name,sum_of_total_value,z_score\nRutaceae,1.5,1.0\nSapindaceae,0.1,-1.0\n",
                encoding="utf-8",
            )
            (root / "candidates.txt").write_text("Rutaceae\n", encoding="utf-8")
            args = SimpleNamespace(
                project="sample", level="family", stage_id=None,
                matrix_dir=matrix, predictions=root / "predictions.csv",
                evidence=None, candidates=root / "candidates.txt", gene_list=None,
                z_threshold=0.0,
            )
            stage = build_stage(args)
            self.assertEqual(stage["genes"][0]["matches"][0]["taxon"], "Rutaceae")
            self.assertEqual(stage["counts"]["query_exons"], 1)
            self.assertTrue(stage["predictions"][0]["pipeline_selected"])
            self.assertEqual(stage["predictions"][0]["locus_votes"], 1)
            self.assertEqual(stage["predictions"][0]["query_votes"], 1)
            stage["predictions"][0]["locus_votes"] = 0  # Simulate an older export.
            stage["predictions"][0]["query_votes"] = 0
            report = {"schema_version": SCHEMA, "sample_id": "sample", "stages": [stage]}
            write_report(report, root / "sample.pea.json.gz", root / "sample.report.md")
            self.assertEqual(stage["predictions"][0]["locus_votes"], 1)
            self.assertLess((root / "sample.pea.json.gz").stat().st_size, 10000)
            self.assertIn("Rutaceae", markdown(report))


if __name__ == "__main__":
    unittest.main()
