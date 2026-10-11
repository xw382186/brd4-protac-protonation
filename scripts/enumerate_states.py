#!/usr/bin/env python3
"""Replay the SI's state-enumeration settings without modifying an installation.

This is a new reproducibility wrapper around Dimorphite-DL 2.0.2, not the
historical preparation script and not a numerical pKa prediction method.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = {"ph_min": 7.4, "ph_max": 7.4, "precision": 0.0, "max_variants": 32}
REQUIRED_IDS = {"5", "6", "7", "8"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_inputs(path: Path) -> list[dict[str, str]]:
    records = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f"Input line {line_number}: expected SMILES and compound ID")
        records.append({"compound": parts[1], "input_smiles": parts[0]})
    ids = [record["compound"] for record in records]
    if len(ids) != len(set(ids)) or set(ids) != REQUIRED_IDS:
        raise ValueError("The archived input set must contain compounds 5, 6, 7, and 8 once each")
    return sorted(records, key=lambda record: int(record["compound"]))


def worker() -> None:
    """Run only in a subprocess whose PYTHONPATH starts with the temporary copy."""
    import dimorphite_dl
    from rdkit import Chem, rdBase

    if dimorphite_dl.__version__ != "2.0.2":
        raise RuntimeError(f"Expected Dimorphite-DL 2.0.2; found {dimorphite_dl.__version__}")
    payload = json.load(sys.stdin)
    records = []
    for item in payload["inputs"]:
        input_molecule = Chem.MolFromSmiles(item["input_smiles"])
        if input_molecule is None:
            raise ValueError(f"Invalid input SMILES for compound {item['compound']}")
        variants = dimorphite_dl.protonate_smiles(item["input_smiles"], **PARAMETERS)
        unique_states = {}
        for variant in variants:
            molecule = Chem.MolFromSmiles(variant)
            if molecule is None:
                raise ValueError(f"Invalid enumerated SMILES for compound {item['compound']}")
            canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
            unique_states[canonical] = Chem.GetFormalCharge(molecule)
        if not unique_states:
            raise RuntimeError(f"No states returned for compound {item['compound']}")
        records.append(
            {
                **item,
                "number_of_unique_states": len(unique_states),
                "states": [
                    {"smiles": smiles, "formal_charge": unique_states[smiles]}
                    for smiles in sorted(unique_states)
                ],
            }
        )
    json.dump(
        {
            "software": {"dimorphite_dl": dimorphite_dl.__version__, "rdkit": rdBase.rdkitVersion},
            "compounds": records,
        },
        sys.stdout,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "data/molecular_inputs.smi")
    parser.add_argument("--rules", type=Path, default=ROOT / "rules/site_substructures.smarts")
    parser.add_argument("--output", type=Path, default=ROOT / "results/enumerated_states")
    args = parser.parse_args()
    inputs = read_inputs(args.inputs)
    if importlib.metadata.version("dimorphite-dl") != "2.0.2":
        raise RuntimeError("Install the pinned optional requirements: Dimorphite-DL must be 2.0.2")
    package_spec = importlib.util.find_spec("dimorphite_dl")
    if package_spec is None or not package_spec.submodule_search_locations:
        raise RuntimeError("Cannot locate the Dimorphite-DL package directory")
    source_package = Path(next(iter(package_spec.submodule_search_locations)))
    with tempfile.TemporaryDirectory(prefix="brd4-dimorphite-replay-") as temporary:
        temporary_root = Path(temporary)
        copied_package = temporary_root / "dimorphite_dl"
        shutil.copytree(source_package, copied_package, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        target_rules = copied_package / "smarts/site_substructures.smarts"
        if not target_rules.is_file():
            raise RuntimeError("Unexpected Dimorphite-DL package layout; no SMARTS resource found")
        shutil.copyfile(args.rules, target_rules)
        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(temporary_root) + os.pathsep + environment.get("PYTHONPATH", "")
        # The installed package remains untouched; the worker imports the copy.
        process = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--_worker"],
            input=json.dumps({"inputs": inputs}), text=True, capture_output=True,
            env=environment, check=False,
        )
        if process.returncode:
            raise RuntimeError(f"Enumeration worker failed:\n{process.stderr.strip()}")
        result = json.loads(process.stdout)
    result.update(
        {
            "purpose": "Protonation-state enumeration; not numerical pKa estimation or automatic working-state selection",
            "parameters": PARAMETERS,
            "input_sha256": sha256(args.inputs),
            "rule_table_sha256": sha256(args.rules),
            "ordering": "Canonical isomeric SMILES sorted lexically; order is not a probability ranking",
        }
    )
    args.output.mkdir(parents=True, exist_ok=True)
    json_path = args.output / "enumerated_states.json"
    csv_path = args.output / "enumerated_states.csv"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound", "state_index", "formal_charge", "smiles"])
        writer.writeheader()
        for compound in result["compounds"]:
            for index, state in enumerate(compound["states"], 1):
                writer.writerow({"compound": compound["compound"], "state_index": index, **state})
    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    try:
        if sys.argv[1:] == ["--_worker"]:
            worker()
        else:
            main()
    except (OSError, ValueError, RuntimeError, importlib.metadata.PackageNotFoundError) as error:
        raise SystemExit(str(error)) from error
