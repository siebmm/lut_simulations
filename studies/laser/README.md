# Inputs and helpers for the paper replays

This directory contains only the inputs and Python helpers used by the
repository's main-result replay, paper-figure replay and independent weld
physics check. Use the public commands in `tools/`; new models and parameter
studies use `tools/analyze.py`.

`replay_inputs.json` is the dependency index. It lists the required helper
modules and maps each logical data path to one canonical file or a member of
`paper/evidence.zip`. Hashes verify that shared inputs preserve the exact
scientific bytes. The replay workspace reconstructs the historical relative
paths without retaining duplicate copies in Git.

| Retained content | Purpose |
|---|---|
| `paper/comsol_chapter_20261101/` | Reference arrays for all 384 baseline configurations, calibration exports, figure generators and array-to-metric checks |
| `paper/explanatory_saft_20261101/`, `paper/combined_explanatory_saft_20261101/` | Coupling, coverage, source-count and path diagrams |
| `paper/filter_study_20261101/study.py` | Exact published filtering, native execution and scoring helpers |
| `paper/partial_improvement_20261102/` | Reconstruction, comparison and plotting of the final partial-weld figures |
| `src/` | Repository discovery and independent filtering, boundary-response and ray/focusing checks |
| Remaining indexed data | Canonical inputs not already available elsewhere or in the supplement |

The historical directory names identify the source studies; no additional
manuscripts live here. The sole paper is `../../paper/paper.pdf`. All 20,412
acquisition/path comparison rows, including unsuccessful outcomes, remain in
`../../paper/evidence.zip`.

`tools/verify_repository.py` rejects unindexed study files and byte-identical
copies. Runtime results, logs and temporary workspace copies belong under
`output/` or `build/`, both excluded from Git.
