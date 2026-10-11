#!/usr/bin/env python3
"""Reproduce the archived project pKa calculation without its original writes.

Dependencies: Python 3.10+, NumPy, RDKit. The verified 2026-10-05 reproduction
used Python 3.11.16, NumPy 2.4.6 and RDKit 2026.03.5. Dimorphite-DL is not a
dependency of this numerical-regression wrapper: state enumeration is a
separate calculation documented by the supplied SMARTS source files.

Usage:
    python reproduce_pka.py --output-dir ./pka_reproduction_run

The output directory must not already exist. If omitted, a unique timestamped
directory is created in the current working directory. The supplied original
script is hash-checked and executed unchanged only up to its final json.dump
section. The wrapper writes JSON, CSV, SMILES and provenance files exclusively
inside the new output directory; no historical project output is overwritten.

The reference pKa inputs are literal values from the archived script, whose
per-entry primary citations have not been verified. LOO RMSE measures error
within that reference set; it is not a PROTAC-specific confidence interval.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
import sys


EXPECTED_SHA256 = "305e0cf22629f5e0afbc9c1e14d3cba1d8e3b73a99369ffbf63998c8c622c650"
ORIGINAL_REMOTE_PATH = "/scratch/wx2237/kaitlyn/BRD4/(1-12)protonation_20260906/pka_rdkit_lfer.py"


def write_json(path: Path, data: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        help="New directory for all reproduction outputs; must not already exist.",
    )
    parser.add_argument(
        "--source-script", type=Path,
        default=Path(__file__).resolve().parents[1] / "source" / "pka_rdkit_lfer.py",
        help="Archived original script (must match the verified SHA256).",
    )
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    try:
        import numpy as np
        from rdkit import rdBase
    except ImportError as error:
        parser.error(f"Missing dependency: {error}. Install NumPy and RDKit in the chosen environment.")

    source_path = args.source_script.expanduser().resolve()
    try:
        source_bytes = source_path.read_bytes()
    except OSError as error:
        parser.error(f"Cannot read the archived source script: {error}")
    source_hash = hashlib.sha256(source_bytes).hexdigest()
    if source_hash != EXPECTED_SHA256:
        parser.error("Source SHA256 differs from the verified original; refusing to execute altered code.")
    source_text = source_bytes.decode("utf-8")
    marker = "\njson.dump("
    if source_text.count(marker) != 1:
        parser.error("Expected exactly one final json.dump output section in the archived source.")
    calculation_prefix = source_text.split(marker, 1)[0]

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    output_dir = (args.output_dir or Path(f"pka_reproduction_{timestamp}")).expanduser().resolve()
    if output_dir.exists():
        parser.error(f"Output directory already exists: {output_dir}. Choose a new directory.")

    namespace = {"__name__": "__pka_reproduction__", "__file__": str(source_path)}
    captured_stdout = io.StringIO()
    # Only the verified original calculation is executed. Its hard-coded final
    # JSON/CSV writes are outside this prefix and cannot run here.
    with contextlib.redirect_stdout(captured_stdout):
        exec(compile(calculation_prefix, str(source_path), "exec"), namespace)

    beta = namespace["beta"]
    reference_records = []
    for row, predicted, residual, loo_residual in zip(
        namespace["rows"], namespace["pred"], namespace["resid"], namespace["loo"]
    ):
        record = dict(row)
        record.update(
            predicted_pKa=float(predicted), residual=float(residual),
            loo_residual=float(loo_residual),
            loo_predicted_pKa=float(row["pka"] - loo_residual),
        )
        reference_records.append(record)

    compound_records = {}
    for compound, smiles in namespace["SM"].items():
        n_subst, sigma, ring_size, charge, index, count, ring_hetero = namespace["descriptors"](smiles)
        design_row = namespace["design"]({"sigma": sigma, "n_subst": n_subst, "rhet": ring_hetero})
        pka = float(np.asarray(design_row) @ beta)
        fraction_percent = float(100.0 / (1.0 + 10.0 ** (7.4 - pka)))
        compound_records[compound] = {
            "smiles": smiles, "n_subst": int(n_subst), "sigma_ind": float(sigma),
            "ring_size": int(ring_size), "gasteiger_N_charge": float(charge),
            "basic_N_index_0based": int(index), "n_basic_N": int(count),
            "ring_hetero": float(ring_hetero), "design_row": list(map(float, design_row)),
            "pKa_unrounded": pka,
            "protonated_fraction_percent_pH7_4": fraction_percent,
            "original_script_rounded_record": namespace["out"][compound],
        }

    results = {
        "provenance": "New numerical reproduction of the archived original script; final file-writing section omitted. Not an original archived output.",
        "script_remote_path": ORIGINAL_REMOTE_PATH,
        "script_repository_path": "source/pka_rdkit_lfer.py", "script_sha256": source_hash,
        "reproduction_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "rdkit_version": rdBase.rdkitVersion, "numpy_version": np.__version__,
        "alpha": float(namespace["ALPHA"]),
        "basis": ["intercept", "sigma_ind", "secondary_indicator", "tertiary_indicator", "ring_hetero"],
        "beta": list(map(float, beta)), "r2": float(namespace["r2"]),
        "fit_rmse": float(namespace["rmse"]), "loo_rmse": float(namespace["loo_rmse"]),
        "loo_mae": float(namespace["loo_mae"]),
        "reference_records": reference_records, "compound_records": compound_records,
        "captured_original_stdout": captured_stdout.getvalue(),
        "interpretation_notes": [
            "Reference pKa values are the literal REF inputs in the original script; per-entry primary bibliographic sources have not been verified.",
            "The captured original stdout retains historical labels. Those labels do not establish verified experimental-source provenance.",
            "LOO RMSE is an internal reference-set error statistic, not a compound-specific confidence interval. Original pka_lo/pka_hi fields are plus/minus one LOO RMSE and are not confidence limits.",
            "OLD_ANCHOR values are historical comparisons and do not enter model fitting or reported model predictions.",
            "Gasteiger nitrogen charge and ring size are diagnostic outputs, not regression features.",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=False)
    write_json(output_dir / "pka_model_reproduction.json", results)
    reference_fields = ["name", "smi", "pka", "n_subst", "sigma", "ring", "q", "rhet",
                        "predicted_pKa", "residual", "loo_residual", "loo_predicted_pKa"]
    with (output_dir / "calibration_inputs_and_loo.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=reference_fields)
        writer.writeheader()
        writer.writerows(reference_records)
    with (output_dir / "compound_predictions.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["compound", "basic_N_index_0based", "n_subst", "ring_size", "sigma_ind",
                         "ring_hetero", "pKa_unrounded", "Table_S4A_pKa",
                         "protonated_fraction_percent_pH7_4", "Table_S4A_fraction_percent"])
        for compound, record in compound_records.items():
            writer.writerow([compound, record["basic_N_index_0based"], record["n_subst"], record["ring_size"],
                             record["sigma_ind"], record["ring_hetero"], record["pKa_unrounded"],
                             f"{record['pKa_unrounded']:.2f}", record["protonated_fraction_percent_pH7_4"],
                             f"{record['protonated_fraction_percent_pH7_4']:.1f}"])
    with (output_dir / "molecular_inputs.smi").open("x", encoding="utf-8") as handle:
        for compound, smiles in namespace["SM"].items():
            handle.write(f"{smiles}\t{compound}\n")
    (output_dir / "reproduction_stdout.txt").write_text(captured_stdout.getvalue(), encoding="utf-8")
    print(f"Reproduction outputs: {output_dir}")
    print("pKa (compounds 5–8): " + ", ".join(f"{record['pKa_unrounded']:.2f}" for record in compound_records.values()))
    print(f"LOO RMSE: {results['loo_rmse']:.12f}; per-entry primary sources remain unverified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
