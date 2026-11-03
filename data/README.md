# Canonical SAFT acquisitions

These are the archived weld regression datasets. New models use the separate
[geometry/acquisition contract](../saft/FORMAT.md); they do not need this weld
layout or these source/receiver/sample counts. Keep generated data in `output/`
or an external dataset directory and pass it explicitly to `tools/analyze.py`.

`full_weld/` and `partial_weld/` each contain `dataset.json` and 60 CSVs:
five sources × reference/LOF × left/right/weld receiver groups × vx/vy.
Both manifests include the material, geometry, source/receiver coordinates,
original completion evidence, audit hashes and repair intervals.

The full weld has 337 receivers; the partial weld has 331. Histories contain
5,001 samples at 10 ns intervals (0–50 microseconds). These are the canonical
resampled histories used in the published SAFT analysis. Larger native-time
solver tables and solved fields are unnecessary for replay and are omitted.

All 120 CSVs are copied byte-for-byte from the source study and checked against
their original SHA-256 hashes by `python tools/verify_repository.py`. Historical
absolute paths in the manifests document provenance; the engine resolves the
recording paths within these folders.

Combined recording size: 1,572,025,590 bytes. Git LFS tracks the CSVs. The custom
full-weld processed bank is a generated output, not a second stored acquisition.
