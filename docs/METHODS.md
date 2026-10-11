# Calculation procedure and provenance

## 1. Scope and source materials

This procedure implements the protonation assessment reported in the latest supplied `Supplemental Information-2.docx`, specifically the Protonation Assessment, Calibration Inputs and Model Definition, and Model Performance and Reproducibility sections and Supplementary Tables S4A–S4C. The repository includes the archived calculation and rule files referenced by that method.

Three operations should be distinguished:

| Operation | Input | Output |
| --- | --- | --- |
| Candidate-state enumeration | Molecular SMILES and modified Dimorphite-DL rules | Candidate protonation states under the specified enumeration settings |
| Selected-amine pKa regression | Molecular connectivity and 20 literal calibration entries | A numerical conjugate-acid pKa estimate for the selected amine |
| Henderson–Hasselbalch conversion | Unrounded pKa and pH 7.4 | Selected-amine protonated fraction under a single-site model |

Dimorphite-DL's generic-rule pKa parameters do not supply the reported 9.59, 10.01, 9.51, and 9.93 pKa values. Those numbers are the regression predictions described below. The regression does not use MD frames, coordinates, energies, or simulation trajectories.

The original `source/pka_rdkit_lfer.py` is preserved byte-for-byte. Its SHA256 is:

```text
305e0cf22629f5e0afbc9c1e14d3cba1d8e3b73a99369ffbf63998c8c622c650
```

The supplied full-precision reference files are numerical reproductions, rather than recovered historical output files. The new `scripts/reproduce.py` adapter retains the archived calculation while providing portable paths, structured outputs, and current SI table labels. Reproduction provenance records identify the actual execution environment.

## 2. Molecular inputs and selected sites

The unmodified SMILES for compounds 5–8 are provided in `data/molecular_inputs.smi`. Their character order matters for reproducing the zero-based RDKit atom indices quoted in Table S4A: rewriting a SMILES can change atom numbering while preserving the molecule.

The numerical regression selects aliphatic amines using this archived SMARTS:

```text
[NX3;!$(N-a);!$(N[CX3]=[OX1,SX1]);!$(N[SX4](=O)=O);!$(N=*);!$(N#*);!$(N[NX3,NX2])]
```

This excludes aryl-bound, acyl/thioacyl-bound, sulfonyl-bound, multiply bonded, and specified nitrogen-adjacent nitrogen environments. Each of compounds 5–8 has exactly one matching nitrogen under this rule. In 7 and 8 the aryl-bound nitrogen is excluded.

| Compound | Selected site | Zero-based N index | Number of matching N atoms | Heavy-atom neighbours of selected N | Smallest N-containing ring |
| --- | --- | ---: | ---: | ---: | ---: |
| 5 | Azetidine N | 26 | 1 | 3 | 4 atoms |
| 6 | Azetidine N | 27 | 1 | 3 | 4 atoms |
| 7 | Piperidine N | 17 | 1 | 3 | 6 atoms |
| 8 | Piperidine N | 27 | 1 | 3 | 6 atoms |

For calibration molecules with multiple matches, the original script selects the nitrogen with the smallest inductive descriptor. This applies to its diamine entry. The script adds hydrogens when calculating a diagnostic Gasteiger charge; the selected heavy-atom index is retained.

## 3. Calibration inputs

The regression uses the following 20 literal entries from the original script's `REF` list, reproduced in `data/calibration_inputs.csv` and SI Table S4B.

| Entry | Input SMILES | Input pKa |
| --- | --- | ---: |
| Ammonia | `N` | 9.25 |
| Methylamine | `CN` | 10.66 |
| Ethylamine | `CCN` | 10.70 |
| Dimethylamine | `CNC` | 10.73 |
| Trimethylamine | `CN(C)C` | 9.80 |
| Triethylamine | `CCN(CC)CC` | 10.75 |
| Pyrrolidine | `C1CCNC1` | 11.27 |
| Piperidine | `C1CCNCC1` | 11.12 |
| Azetidine | `C1CNC1` | 11.29 |
| N-Methylpyrrolidine | `CN1CCCC1` | 10.32 |
| N-Methylpiperidine | `CN1CCCCC1` | 10.08 |
| Morpholine | `C1COCCN1` | 8.36 |
| N-Methylmorpholine | `CN1CCOCC1` | 7.38 |
| Ethanolamine | `NCCO` | 9.50 |
| 2-Methoxyethylamine | `NCCOC` | 9.40 |
| Diethanolamine | `OCCNCCO` | 8.88 |
| Triethanolamine | `OCCN(CCO)CCO` | 7.76 |
| Ethylenediamine (pKa1) | `NCCN` | 9.92 |
| Glycinamide | `NCC(N)=O` | 8.20 |
| 3-Aminopropanamide | `NCCC(N)=O` | 9.30 |

The ethylenediamine label is retained as written in the source. The supplied script describes these as tabulated inputs, but it and the supplied SI do not assign primary bibliographic references to individual entries. Their computational provenance is the literal `REF` list; this repository does not assert independently verified experimental provenance for each input.

## 4. Descriptors and regression equation

For the selected nitrogen N, the design vector is:

```text
x = [1, D, I2, I3, R]
pKa = b0 + b1*D + b2*I2 + b3*I3 + b4*R
```

### Distance-weighted inductive descriptor D

```text
D = sum[(chi_i - 2.55) * 0.40**(d(N,i) - 1)]
    over eligible heavy atoms i with 1 <= d(N,i) <= 8
```

`d(N,i)` is the shortest-path bond distance in the RDKit molecular graph. The selected nitrogen is excluded. Carbon contributes zero because its electronegativity is the reference value. Only elements encoded in the source dictionary contribute; other elements are skipped. No descriptor scaling is applied.

| Element | Encoded Pauling electronegativity |
| --- | ---: |
| C | 2.55 |
| N | 3.04 |
| O | 3.44 |
| S | 2.58 |
| F | 3.98 |
| Cl | 3.16 |
| Br | 2.96 |
| I | 2.66 |
| H | 2.20 |

The dictionary includes H, but D is calculated on the input molecular graph before adding explicit hydrogens. The supplied SMILES therefore yield the heavy-atom descriptor described in the SI. The attenuation parameter is fixed at `ALPHA = 0.40`; it is a supplied model parameter, not refitted during this procedure.

### Nitrogen substitution indicators

`I2 = 1` if the selected nitrogen has two heavy-atom neighbours, and zero otherwise. `I3 = 1` if it has three heavy-atom neighbours, and zero otherwise. Thus ammonia and primary amines use `I2 = I3 = 0`. All four PROTAC selected sites have `I2 = 0` and `I3 = 1`.

### Ring descriptor R

For every ring returned by RDKit's `AtomRings()` that contains the selected nitrogen, the script sums `(chi_i - 2.55) * 0.40**(d(N,i) - 1)` over the other encoded atoms in that ring. It uses the molecular shortest-path distance and has no eight-bond cutoff. Contributions are accumulated separately for each qualifying ring. This additional ring term is distinct from D. The four PROTAC selected sites have `R = 0` in the supplied calculation.

### Diagnostic fields

The script also calculates the smallest nitrogen-containing ring size and the Gasteiger charge on the selected nitrogen. Neither is a regression feature. Historical `OLD_ANCHOR` numbers are used only for printed comparisons and do not enter the design matrix, regression, predicted pKa, or Henderson–Hasselbalch fractions.

## 5. Ordinary least-squares fit

For each calibration molecule, the script constructs its design vector, stacks these vectors into the 20 × 5 matrix X, and stacks the literal pKa inputs into y. It fits:

```python
beta, *_ = numpy.linalg.lstsq(X, y, rcond=None)
prediction = X @ beta
residual = y - prediction
```

No regularization, feature scaling, or external fitted coefficients are used. The SI coefficients rounded to four decimal places are:

| Coefficient | Term | Value |
| --- | --- | ---: |
| b0 | Intercept | 10.1754 |
| b1 | D | −6.6536 |
| b2 | I2 | 0.8827 |
| b3 | I3 | 0.1139 |
| b4 | R | −13.0356 |

The supplied reference reproduction retains full precision:

```text
beta = [10.175419379262436,
        -6.6536099388676595,
         0.8826948104437108,
         0.11394048632617047,
       -13.035554581127965]
```

Predictions use the fitted full-precision coefficients, rather than the rounded table values.

## 6. Fit statistics and leave-one-out validation

Training statistics use the 20 fitted values:

```text
R2 = 1 - sum(residual**2) / sum((y - mean(y))**2)
RMSE = sqrt(mean(residual**2))
```

For each entry k, leave-one-out validation removes that row, refits the same five-term model to the remaining 19 entries, and predicts the held-out entry. It does not select new descriptors or change the fixed attenuation parameter. The exported leave-one-out residual is `input pKa - held-out predicted pKa`.

| Statistic | Reference reproduction | SI rounding |
| --- | ---: | ---: |
| Training R2 | 0.8928667194722474 | 0.893 |
| Training RMSE | 0.3733056700097274 | 0.373 |
| Leave-one-out RMSE | 0.49324256829621593 | 0.493 |
| Leave-one-out MAE | 0.42233338240011575 | 0.422 |

These statistics describe the supplied calibration set. They do not establish experimental accuracy for compounds 5–8 or compound-specific confidence intervals. The archived script's historical `pka_lo` and `pka_hi` fields are the prediction minus/plus one leave-one-out RMSE, respectively; they are not statistical confidence bounds. Captured original console output retains its original wording.

## 7. PROTAC predictions and protonated fractions

For each input SMILES in compounds 5–8, the script selects the nitrogen, calculates the same descriptors, constructs x, and predicts `pKa = x @ beta`. It then evaluates the single-site equilibrium `BH+ ⇌ B + H+` at pH 7.4:

```text
Protonated amine (%) = 100 / (1 + 10**(7.4 - pKa))
```

The calculation uses unrounded pKa. Only the final display is rounded to two decimal places for pKa and one decimal place for the percentage.

| Compound | D | R | Unrounded pKa | Unrounded protonated amine (%) |
| --- | ---: | ---: | ---: | ---: |
| 5 | 0.10479667200000006 | 0 | 9.592083687209152 | 99.3615386368367 |
| 6 | 0.041276416000000024 | 0 | 10.01472269385017 | 99.75777216480948 |
| 7 | 0.11734067200000008 | 0 | 9.508620804135996 | 99.22730095710408 |
| 8 | 0.053820416000000024 | 0 | 9.931259810777014 | 99.70659731452041 |

These estimates support the working +1 assignments for compounds 5–8 described in the SI. Compounds 1, 2, 3, 4, 9, 10, 11, and 12 retain working neutral assignments; the SI reports dashes rather than numerical pKa estimates for them. A dash is not a predicted pKa of zero.

The percentage concerns the selected amine under a single-site model. It is not a complete distribution over every potentially ionizable site in the molecule. Working net charge identifies the assigned ligand state and is not a population-weighted molecular charge.

## 8. Dimorphite-DL candidate-state enumeration

The SI specifies locally modified Dimorphite-DL **2.0.2**, with:

```python
protonate_smiles(
    smi,
    ph_min=7.4,
    ph_max=7.4,
    precision=0.0,
    max_variants=32,
)
```

The archived original and modified rule tables are supplied in `rules/`. The reconstructed difference file compares those actual rule tables; it is not a recovered historical patch.

The generic amine rule changed from:

```text
[C:1]-[NX3+0:2]
```

to:

```text
[C:1]-[NX3+0;!$(N[CX3]=[OX1,SX1]);!$(N[SX4](=O)=O);!$(N[CX3]=[NX2]);!$(Na);!$(N[NX3,NX2]);!$(N=*);!$(N#*):2]
```

The rule's ionization-site index is 1, zero-based within the SMARTS match. Its recorded pKa mean (`8.159107682388349`) and standard deviation (`2.5183597445318147`) are unchanged by the modification. These rule parameters govern enumeration; they are distinct from the numerical regression pKa values.

| Rule-table record | SHA256 |
| --- | --- |
| Archived original `site_substructures.smarts.orig` | `32afd9a2de7c55bab70e63f9760b4d731ea26545587c122736d80f52228b1bac` |
| Archived customized table | `e2d3254ced2a2a744ca8f11d7d2eedde7c0c0267c2692bb539bcf78f306bf580` |
| Deposited `rules/site_substructures.smarts`, with attribution comments | `1fda9f2ba85a36437c58e2c4665e932a52d95153d328644912794471c0c0c7e3` |

The deposited customized table adds three attribution/modification-comment lines. Its non-comment rule rows are unchanged from the archived customized table. Enumeration output records the checksum of the actual deposited rule table used.

The numerical reproduction command does not require Dimorphite-DL and does not claim to recover the historical enumerated state list. The separate optional runner replays the recorded enumeration method for the available inputs, compounds 5–8. Use a dedicated environment so that Dimorphite-DL 2.0.2's declared `rdkit<2026` dependency does not alter the SI numerical pKa environment:

```bash
python -m venv .venv-enumeration
source .venv-enumeration/bin/activate
python -m pip install -r requirements-enumeration.txt
python scripts/enumerate_states.py --inputs data/molecular_inputs.smi --output results/enumerated_states
```

The wrapper temporarily copies the installed package and substitutes the deposited rule table in that copy, leaving the installation unchanged. It writes `enumerated_states.json` and `enumerated_states.csv`, recording parameters, runtime versions, input and rule checksums, unique canonical isomeric state SMILES, and net formal charges. Output ordering is not a probability ranking.

The tested Dimorphite-DL 2.0.2 / RDKit 2025.09.6 / Loguru 0.7.3 environment returned one unique +1 state for each of 5, 6, 7, and 8. A separate compatibility run retaining RDKit 2026.03.5 returned identical state records; that RDKit version is outside Dimorphite-DL 2.0.2's declared dependency bound and is not the pinned optional environment. These new enumeration results are not recovered historical outputs, experimental states, or numerical pKa estimates. See [ENUMERATION.md](ENUMERATION.md) for the complete optional workflow.

## 9. Running the portable numerical reproduction

From the repository root:

```bash
python -m pip install -r requirements.txt
python scripts/reproduce.py --output-dir results/reproduced
python -m unittest discover -s tests -v
```

The environment reported in the SI is Python 3.11.16, RDKit 2026.03.5, and NumPy 2.4.6. The runtime versions for each new reproduction are written to its JSON provenance record. Minor floating-point differences across platforms should be assessed against the full-precision reference and the SI rounding, without silently changing the model inputs.

The adapter:

1. Reads the archived script and verifies its SHA256.
2. Executes the original calculation in memory up to, but not including, its final `json.dump` file-writing section.
3. Retains the original 20 calibration inputs and the four molecular inputs.
4. Recomputes the unrounded predictions, selected atom indices, descriptors, and protonated fractions.
5. Exports reference-set fitted and held-out predictions.
6. Writes new files to the requested new output directory; it refuses an existing directory.

| Output | Meaning |
| --- | --- |
| `pka_model_reproduction.json` | Source checksum, software versions, coefficients, statistics, calibration records, compound records, and captured original console output |
| `calibration_inputs_and_loo.csv` | All 20 inputs, descriptors, training predictions/residuals, and held-out predictions/residuals |
| `compound_predictions.csv` | Selected atom indices, descriptors, full-precision pKa and fractions, and rounded SI values |
| `molecular_inputs.smi` | Exact compound 5–8 SMILES used by the executed source |
| `reproduction_stdout.txt` | Captured original calculation console output |

The records in `results/reference/` provide comparison values. Code provenance, numerical reproducibility, and bibliographic provenance of calibration inputs are separate properties; the supplied data establish the first two without assigning undocumented citations to individual reference values.
