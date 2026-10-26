import csv
import sys
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_analysis import FIGURES, RESULTS, fit_standard_curve, read_csv, read_project_data


class EnzymologyDataTests(unittest.TestCase):
    def test_standard_curve_fit_matches_tabulated_points(self):
        standards = read_csv(ROOT / "data" / "maltose_standard_curve.csv")
        fit = fit_standard_curve(standards)
        self.assertAlmostEqual(fit["slope"], 0.81833, places=4)
        self.assertAlmostEqual(fit["intercept"], -0.01578, places=4)
        self.assertGreater(fit["r_squared"], 0.99)

    def test_all_source_readings_are_transcribed_once(self):
        _, assays, _ = read_project_data()
        self.assertEqual(len(assays), 52)
        counts = Counter((row["experiment"], row["condition_value"], row["enzyme"]) for row in assays)
        self.assertEqual(set(counts.values()), {1})

    def test_reported_peak_conditions(self):
        _, assays, _ = read_project_data()
        for enzyme in ("Commercial fungal amylase", "Salivary amylase"):
            ph_rows = [r for r in assays if r["experiment"] == "pH" and r["enzyme"] == enzyme]
            temp_rows = [r for r in assays if r["experiment"] == "temperature" and r["enzyme"] == enzyme]
            ca_rows = [r for r in assays if r["experiment"] == "calcium_chloride" and r["enzyme"] == enzyme]
            self.assertEqual(max(ph_rows, key=lambda r: float(r["absorbance_540nm"]))["condition_value"], "7")
            self.assertEqual(max(temp_rows, key=lambda r: float(r["absorbance_540nm"]))["condition_value"], "37")
            self.assertEqual(max(ca_rows, key=lambda r: float(r["absorbance_540nm"]))["condition_value"], "10")

    def test_output_files_were_generated(self):
        expected = [
            FIGURES / "figure_01_maltose_standard_curve.png",
            FIGURES / "figure_02_ph_response.png",
            FIGURES / "figure_03_temperature_response.png",
            FIGURES / "figure_04_calcium_chloride_response.png",
            FIGURES / "figure_05_starch_response.png",
            FIGURES / "figure_06_ethanol_response.png",
            RESULTS / "recalculated_product_yield.csv",
            ROOT / "reports" / "final_report.pdf",
        ]
        self.assertTrue(all(path.exists() and path.stat().st_size > 0 for path in expected))
        with (RESULTS / "recalculated_product_yield.csv").open(newline="", encoding="utf-8") as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 52)


if __name__ == "__main__":
    unittest.main()
