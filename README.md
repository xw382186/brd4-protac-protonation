# BRD4 PROTAC protonation calculations

Reproducible numerical pKa calculations and protonation-method records for the Supplementary Information of *Computational Evaluation of Linker-Varied BRD4-Targeting PROTACs Using Deep-Learning Activity Classification and Molecular Dynamics*.

[中文说明](README_zh.md) · [Detailed calculation procedure](docs/METHODS.md) · [Optional state enumeration](docs/ENUMERATION.md)

This repository provides the archived RDKit regression script, its 20 calibration inputs, the input SMILES for compounds 5–8, modified Dimorphite-DL ionization rules, full-precision reference results, and a portable reproduction script. The scientific scope is the SI protonation assessment and Supplementary Tables S4A–S4C.

## Reproduce the numerical results

The SI reports a verified environment of **Python 3.11.16, RDKit 2026.03.5, and NumPy 2.4.6**. The numerical calculation requires NumPy and RDKit; Dimorphite-DL is used separately for state enumeration.

```bash
python -m pip install -r requirements.txt
python scripts/reproduce.py --output-dir results/reproduced
python -m unittest discover -s tests -v
```

Choose a new output directory for each run. The reproduction script verifies the archived source checksum, retains its descriptors, calibration inputs, fitting, and prediction calculations, and writes new outputs only inside the requested directory. It captures the original calculation output and records the versions actually used. It does not run the original script's final hard-coded file-writing section.

Optional candidate-state enumeration uses a **separate environment**, because Dimorphite-DL 2.0.2 declares `rdkit<2026`. This leaves the numerical pKa reproduction environment unchanged:

```bash
python -m venv .venv-enumeration
source .venv-enumeration/bin/activate
python -m pip install -r requirements-enumeration.txt
python scripts/enumerate_states.py --inputs data/molecular_inputs.smi --output results/enumerated_states
```

The tested optional environment uses RDKit 2025.09.6. With the archived customized rules and SI settings, it returns one unique +1 state for each of compounds 5–8. A compatibility run using RDKit 2026.03.5 returned the same state records. Enumeration does not supply the numerical pKa values or assign state probabilities; see [ENUMERATION.md](docs/ENUMERATION.md).

## What is calculated?

1. RDKit identifies the selected aliphatic amine from molecular connectivity. Compounds 5 and 6 use the azetidine nitrogen; compounds 7 and 8 use the piperidine nitrogen.
2. An ordinary least-squares regression is fitted to the 20 literal pKa inputs embedded in the archived script.
3. The model predicts the selected amine's conjugate-acid pKa using a distance-weighted inductive descriptor, nitrogen substitution indicators, and a ring descriptor.
4. The unrounded pKa is converted to the **selected-amine protonated percentage at pH 7.4** using the single-site Henderson–Hasselbalch equation.
5. Training and leave-one-out results are exported alongside the predictions.

| Compound | Selected nitrogen; zero-based input atom index | Estimated pKa | Selected-amine protonation at pH 7.4 | Working net charge |
| --- | --- | ---: | ---: | ---: |
| 5 | Azetidine N; 26 | 9.59 | 99.4% | +1 |
| 6 | Azetidine N; 27 | 10.01 | 99.8% | +1 |
| 7 | Piperidine N; 17 | 9.51 | 99.2% | +1 |
| 8 | Piperidine N; 27 | 9.93 | 99.7% | +1 |

These are the rounded values in Supplementary Table S4A. Its other compounds retain working neutral assignments and have no reported numerical pKa estimates. Atom indices refer to the unmodified strings in `data/molecular_inputs.smi`.

## Repository contents

| Path | Contents |
| --- | --- |
| `source/pka_rdkit_lfer.py` | Unchanged archived calculation script |
| `scripts/reproduce.py` | Portable adapter that reproduces the numerical calculation and records provenance |
| `scripts/enumerate_states.py` | Optional replay of candidate-state enumeration in a temporary package copy |
| `data/calibration_inputs.csv` | The script's 20 literal calibration entries; corresponding to SI Table S4B |
| `data/molecular_inputs.smi` | Input SMILES for compounds 5–8 |
| `data/si_table_s4a.csv` | SI charge assignments and rounded selected-amine results |
| `rules/` | Original and modified SMARTS rule tables, their reconstructed comparison, and the documented enumeration call |
| `results/reference/` | Supplied full-precision numerical reproduction, reference-set leave-one-out results, and provenance |
| `docs/METHODS.md` | Equations, parameters, site selection, calculation steps, and interpretation |
| `docs/ENUMERATION.md` | Optional enumeration environment, execution, outputs, and validation |
| `tests/` | Reproducibility checks |

## Method provenance and interpretation

The numerical pKa values come from the archived regression's predictions. The 20 calibration pKa values come from its literal `REF` list, also reproduced in SI Table S4B; the supplied material does not assign a primary bibliographic citation to each entry. This repository preserves those inputs rather than replacing them or assigning new sources.

The pKa regression uses molecular connectivity, without MD coordinates. Dimorphite-DL candidate-state enumeration and the regression are separate calculations. The Henderson–Hasselbalch percentage describes the selected amine under a single-site model; a working ligand charge is an assigned state, not a population-weighted molecular charge. Leave-one-out error summarizes performance within the calibration set and is not a confidence interval for a PROTAC prediction.

The archived script retains historical comparison fields and printed labels. `OLD_ANCHOR` values do not enter the fit or predictions; Gasteiger nitrogen charge and ring size are diagnostic outputs rather than regression features. Full details and the source checksum are in [METHODS.md](docs/METHODS.md).
