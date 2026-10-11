# Third-party notices

## Dimorphite-DL 2.0.2

The rule tables `rules/site_substructures.smarts.orig` and
`rules/site_substructures.smarts` derive from Dimorphite-DL, Copyright 2025
durrantlab, distributed under the Apache License, Version 2.0.

Upstream project: <https://github.com/durrantlab/dimorphite_dl>

The complete license supplied with the archived Dimorphite-DL 2.0.2 distribution
is reproduced in `third_party/Dimorphite-DL-LICENSE.md`.

The `.orig` file preserves the upstream table from the preparation archive. The
other table contains the project-specific changes recorded by
`rules/dimorphite_rule_record.txt` and
`rules/dimorphite_rule_reconstructed.diff`. These changes are not attributed to
the upstream authors. The new `scripts/enumerate_states.py` wrapper invokes the
upstream API and applies the archived table in a temporary package copy.

No upstream Python implementation, RDKit binaries, or Loguru implementation is
redistributed here. Those optional dependencies are obtained separately using
`requirements-enumeration.txt` and remain subject to their respective licenses.
