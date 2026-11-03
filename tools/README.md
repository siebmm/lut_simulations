# Standalone workflows

Run from the repository root after installing `requirements.txt` and building
the native tools. All paths are discovered relative to this checkout.

| Command | Result |
|---|---|
| `python tools/build_paper.py` | Compile the sole standalone paper to `paper/paper.pdf` with Tectonic; verify references and record source/PDF hashes |
| `python tools/make_example.py` | Generate a non-weld polygon dataset with deterministic synthetic histories in `output/example_dataset` |
| `python tools/validate_dataset.py DATASET` | Validate a schema-2 geometry/acquisition and every CSV, or check an archived schema-1 dataset's hashes |
| `python tools/analyze.py --dataset DATASET --request REQUEST --output OUTPUT` | General reconstruction, images and summaries, with optional user-defined region contrasts |
| `python -m unittest discover -s tests -v` | General-dataset interface, geometry, sampling and independent numerical regressions; requires built tools |
| `python tools/verify_repository.py` | Snapshot, model and acquisition checksums, Python syntax, final figure assets and complete survey-row count |
| `python tools/verify_interfaces.py` | Nine supported pixels, three signal cases, and agreement across all four native exporters |
| `python studies/laser/src/verify_saft_datasets.py` | Independent ray/focusing physics, channel mapping, air masks, muting and dataset isolation; writes `output/physics_check/check.json` |
| `python tools/reproduce_main.py` | Fresh full/partial published grids and matched unmuted control; checks nine signal-case metric sets |
| `python tools/reproduce_main.py --case full` | Full weld only; regenerates its exact processed bank |
| `python tools/reproduce_main.py --case partial` | Two partial-weld grids only |
| `python tools/reproduce_figures.py` | Baseline/diagram artwork and final partial-review figures in an isolated workspace |
| `python tools/reproduce_figures.py --group review` | Only the five revised partial-review figures |

The paper replay commands use `replay_workspace.py` to materialize the inputs
listed in `studies/laser/replay_inputs.json` under `output/replay/`. Shared arrays
are copied from one canonical file; matching records are read from the paper's
existing supplement. Source recordings remain in `data/` and are not copied.
Generated reconstructions can be reused between replay commands.

Set `SAFT_WORKSPACE` to a fresh directory under `output/` or outside the checkout
for a clean replay. A source/input fingerprint prevents reuse of results after
the replay code or definitions change. Do not run the two paper replay commands
simultaneously; their original filter-bank working paths are shared.

For other specimens and parameter studies, use `analyze.py` and the
[public data contract](../saft/FORMAT.md). The retained study modules are internal
helpers for the fixed paper replays. Historical broad-survey and packaging entry
points have been removed; the paper's complete comparison records remain in
`paper/evidence.zip`.

Run the integrity, interface and physics checks before the larger replay.
The integrity check uses only the Python standard library and diagnoses missing
Git LFS content. The native checks require the built tools and `requirements.txt`.
Verification scripts use assertions: run Python normally, without `-O`.

The checked snapshot is `verification/standalone_manifest.json`. After an
intentional source/data change, use `python tools/verify_repository.py
--write-manifest` to record the new snapshot. Runtime reports go to `output/`.
