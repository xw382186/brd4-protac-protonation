"""Verify the archived calculation against current SI Tables S4A–S4C."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.output = Path(cls.temp.name) / "run"
        cls.command = [sys.executable, str(ROOT / "scripts/reproduce.py"),
                       "--output-dir", str(cls.output)]
        cls.process = subprocess.run(cls.command, capture_output=True, text=True)
        if cls.process.returncode:
            raise RuntimeError(cls.process.stderr + cls.process.stdout)
        cls.result = json.loads((cls.output / "pka_model_reproduction.json").read_text())
        cls.reference = json.loads((ROOT / "results/reference/pka_model_reproduction.json").read_text())

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_unchanged_archived_source(self):
        source = ROOT / "source/pka_rdkit_lfer.py"
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(),
                         "305e0cf22629f5e0afbc9c1e14d3cba1d8e3b73a99369ffbf63998c8c622c650")

    def test_table_s4b_calibration_inputs_and_loo(self):
        with (ROOT / "data/calibration_inputs.csv").open() as handle:
            inputs = list(csv.DictReader(handle))
        self.assertEqual(len(inputs), 20)
        for supplied, actual, reference in zip(inputs, self.result["reference_records"],
                                              self.reference["reference_records"]):
            self.assertEqual(supplied["name"], actual["name"])
            self.assertEqual(supplied["smiles"], actual["smi"])
            self.assertEqual(float(supplied["input_pKa"]), actual["pka"])
            for key in ["sigma", "rhet", "predicted_pKa", "loo_predicted_pKa", "loo_residual"]:
                self.assertAlmostEqual(actual[key], reference[key], places=10)

    def test_table_s4c_coefficients_and_validation(self):
        for actual, expected in zip(self.result["beta"], self.reference["beta"]):
            self.assertAlmostEqual(actual, expected, places=10)
        self.assertEqual([round(v, 4) for v in self.result["beta"]],
                         [10.1754, -6.6536, 0.8827, 0.1139, -13.0356])
        for key in ["r2", "fit_rmse", "loo_rmse", "loo_mae"]:
            self.assertAlmostEqual(self.result[key], self.reference[key], places=10)
        self.assertEqual(round(self.result["loo_rmse"], 3), 0.493)

    def test_table_s4a_predictions_sites_and_single_site_fractions(self):
        with (ROOT / "data/si_table_s4a.csv").open() as handle:
            table = {row["compound"]: row for row in csv.DictReader(handle)}
        for compound, record in self.result["compound_records"].items():
            si = table[compound]
            reference = self.reference["compound_records"][compound]
            self.assertEqual(record["basic_N_index_0based"], int(si["nitrogen_index_0based"]))
            self.assertEqual(record["n_basic_N"], 1)
            self.assertAlmostEqual(record["pKa_unrounded"], reference["pKa_unrounded"], places=10)
            self.assertEqual(f"{record['pKa_unrounded']:.2f}", si["estimated_pKa"])
            self.assertEqual(f"{record['protonated_fraction_percent_pH7_4']:.1f}",
                             si["protonated_percent_pH7_4"])
            f = 100.0 / (1.0 + 10.0 ** (7.4 - record["pKa_unrounded"]))
            self.assertAlmostEqual(f, record["protonated_fraction_percent_pH7_4"], places=12)
        for compound in ["1", "2", "3", "4", "9", "10", "11", "12"]:
            self.assertEqual(table[compound]["estimated_pKa"], "")
            self.assertEqual(table[compound]["working_net_charge"], "0")

    def test_original_smiles_order_and_current_table_labels(self):
        self.assertEqual([line.split() for line in (self.output / "molecular_inputs.smi").read_text().splitlines()],
                         [line.split() for line in (ROOT / "data/molecular_inputs.smi").read_text().splitlines()])
        with (self.output / "compound_predictions.csv").open() as handle:
            columns = next(csv.reader(handle))
        self.assertIn("Table_S4A_pKa", columns)
        self.assertNotIn("Table_S5_pKa", columns)

    def test_existing_outputs_are_not_overwritten(self):
        before = (self.output / "pka_model_reproduction.json").read_bytes()
        repeated = subprocess.run(self.command, capture_output=True, text=True)
        self.assertNotEqual(repeated.returncode, 0)
        self.assertIn("already exists", repeated.stderr)
        self.assertEqual(before, (self.output / "pka_model_reproduction.json").read_bytes())

    def test_archived_rule_change_is_one_amine_row(self):
        original = [line for line in (ROOT / "rules/site_substructures.smarts.orig").read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
        modified = [line for line in (ROOT / "rules/site_substructures.smarts").read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
        self.assertEqual(len(original), len(modified))
        changed = [(a, b) for a, b in zip(original, modified) if a != b]
        self.assertEqual(len(changed), 1)
        a, b = [row.split() for row in changed[0]]
        self.assertEqual(a[0], "Amines_primary_secondary_tertiary")
        self.assertEqual(a[0], b[0])
        self.assertEqual(a[2:], b[2:])


if __name__ == "__main__":
    unittest.main()
