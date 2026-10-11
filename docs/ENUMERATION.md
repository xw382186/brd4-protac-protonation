# Optional replay of protonation-state enumeration

The archived molecular inputs cover compounds **5–8 only**. This optional runner
applies the SI settings to those inputs with Dimorphite-DL 2.0.2:

```python
protonate_smiles(smiles, ph_min=7.4, ph_max=7.4, precision=0.0, max_variants=32)
```

The numerical pKa adapter, state enumeration, and the selection of a +1 working
state are separate operations. This runner does **not** calculate the SI's
compound-specific numerical pKa values, derive them from Dimorphite-DL, or choose
a working state. Enumerated states are not assigned equilibrium probabilities.
`precision=0.0` specifies the pKa-window expansion factor used by the package;
it does not imply that numerical pKa values have zero uncertainty.

## Installation and execution

Use Python 3.10 or newer in a dedicated environment. For example:

```bash
python -m venv .venv-enumeration
source .venv-enumeration/bin/activate
python -m pip install -r requirements-enumeration.txt
python scripts/enumerate_states.py
```

The wrapper copies the installed `dimorphite_dl` package to a temporary
directory, replaces `smarts/site_substructures.smarts` **in that temporary copy**
with the archived project table, and invokes a fresh Python subprocess whose
`PYTHONPATH` starts with the copy. The installed package is never patched. The
temporary copy is deleted when execution finishes. This also prevents a cached
rule-table singleton from retaining a previous set of rules.

Inputs and output location can be set explicitly:

```bash
python scripts/enumerate_states.py \
  --inputs data/molecular_inputs.smi \
  --rules rules/site_substructures.smarts \
  --output results/enumerated_states
```

## Outputs and interpretation

- `enumerated_states.json`: parameters, software versions, input/rule SHA-256
  checksums, original input SMILES, unique canonical isomeric state SMILES, and
  their net formal charges.
- `enumerated_states.csv`: one row per unique enumerated state.

The state index follows lexical ordering of canonical SMILES. It is an output
identifier, not an energetic or probability ranking. The runner rejects missing,
duplicate, or extra compound identifiers because the archived input set is
restricted to 5–8. It does not reconstruct unavailable inputs for other compounds.

This wrapper is newly written for reproducibility. The supplied rule table and
molecular input file are archived preparation inputs, rather than newly inferred
chemical parameters. The upstream `.orig` rule table is retained for comparison;
the project-specific rule changes are recorded separately in `rules/`.

The archived customized table has SHA-256
`e2d3254ced2a2a744ca8f11d7d2eedde7c0c0267c2692bb539bcf78f306bf580`.
The redistributed table adds three prominent attribution/modification-comment
lines and has SHA-256
`1fda9f2ba85a36437c58e2c4665e932a52d95153d328644912794471c0c0c7e3`.
Its non-comment rule rows are unchanged. The output JSON records the checksum of
the actual redistributed table used for that run.

## Validation

The runner was executed against the archived Dimorphite-DL 2.0.2 / RDKit
2025.09.6 / Loguru 0.7.3 dependency set using the available Python 3.11 runtime.
It returned one unique state with net formal charge **+1** for each of compounds
5, 6, 7, and 8 under the archived customized rules and fixed SI settings. These
are reproduced enumeration outputs, not independently measured protonation states
or numerical pKa predictions.

A separate compatibility run exposed only Dimorphite-DL and Loguru from the
archive, retaining the runtime's RDKit 2026.03.5. It returned identical state
records. RDKit 2026.03.5 lies outside the dependency bound declared by
Dimorphite-DL 2.0.2 (`rdkit<2026`), so it is not the pinned reproduction
environment. The archived installed rule table remained unchanged after the
runner finished.

See `THIRD_PARTY_NOTICES.md` and `third_party/Dimorphite-DL-LICENSE.md` for the
upstream rule-table attribution and license.
