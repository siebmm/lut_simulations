# Standalone workflows

Run from the repository root after installing `requirements.txt` and building
the native tools. All paths are discovered relative to this checkout.

| Command | Result |
|---|---|
| `python tools/verify_repository.py` | Snapshot, model and acquisition checksums, Python syntax, final figure assets and complete survey-row count |
| `python tools/verify_interfaces.py` | Nine supported pixels, three signal cases, and agreement across all four native exporters |
| `python studies/laser/src/verify_saft_datasets.py` | Independent ray/focusing physics, channel mapping, air masks, muting and dataset isolation; writes `output/physics_check/check.json` |
| `python tools/reproduce_main.py` | Fresh full/partial published grids and matched unmuted control; checks nine signal-case metric sets |
| `python tools/reproduce_main.py --case full` | Full weld only; regenerates its exact processed bank |
| `python tools/reproduce_main.py --case partial` | Two partial-weld grids only |
| `python tools/reproduce_figures.py` | Baseline/diagram artwork and final partial-review figures in an isolated workspace |
| `python tools/reproduce_figures.py --group review` | Only the five revised partial-review figures |
| `python tools/run_study.py STUDY/SCRIPT.py [arguments]` | Runs an included study script against a writable copy of its protocols/evidence |

`run_study.py` seeds `output/workspace/` once, updates Python source on later
runs, and preserves generated results between commands. It reads acquisitions
and executables from this repository, without copying the large CSVs again.
Use this wrapper for the original study scripts, which write result folders.
Do not run multiple reconstruction studies simultaneously: their original
filter-bank working directories are reused between variants.

Example filter study sequence:

```sh
python tools/run_study.py filter_study_20261101/study.py raw
python tools/run_study.py filter_study_20261101/study.py verify
python tools/run_study.py filter_study_20261101/study.py run
python tools/run_study.py filter_study_20261101/study.py summarise
```

The revised partial-weld atlas uses `partial_improvement_20261102/review.py`
with `verify`, `atlas`, then `analyse`; its settings, confirmations, refinements
and source-subset analyses are adjacent. Full surveys can be substantially more
expensive than the main replay. Their original protocols and complete summary
records are included, so running every survey is not required to read the paper.

Run the integrity, interface and physics checks before the larger replay.
The integrity check uses only the Python standard library and diagnoses missing
Git LFS content. The native checks require the built tools and `requirements.txt`.
Verification scripts use assertions: run Python normally, without `-O`.

The checked snapshot is `verification/standalone_manifest.json`. After an
intentional source/data change, use `python tools/verify_repository.py
--write-manifest` to record the new snapshot. Runtime reports go to `output/`.
