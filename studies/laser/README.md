# Scientific analysis and compact evidence

The original study names and relative layout are retained so frozen protocols
and cross-study imports remain traceable. The canonical acquisitions are stored
once at the repository's `data/`; all native tools share `saft/src/`.

| Folder under `paper/` | Role in the standalone paper |
|---|---|
| `filter_study_20261101` | Exact FIR/IIR/analytic filter definitions, native parity, full-weld selected result and controls |
| `equal_receiver_study_20261101` | Matched partial-weld baseline comparisons |
| `partial_improvement_20261102` | 324-path atlas, settings/admission tests, finer/wider/mute checks and final weighted partial result |
| `partial_acquisition_20261101` | Source/receiver acquisition comparisons |
| `gate_ablation_20261101` | Physical admission-gate controls |
| `partial_geometry_paths_20261030` | Finite-wall/corner path hypotheses and independent geometry checks |
| `comsol_chapter_20261101` | Portable data for 384 image configurations, figure generators and independent array-to-metric verification |
| `explanatory_saft_20261101`, `combined_explanatory_saft_20261101` | Coupling, coverage, source requirements and path illustrations |
| Other retained study folders | Earlier baseline, geometry and feasibility dependencies used by these analyses |

`results/` contains the frozen scientific summaries, settings and protocols.
It omits reproducible per-grid intermediate arrays and repeated image CSVs.
`paper/*/data/` retains compact arrays used by the figure generators, including
the calibration exports. `src/` contains shared numerical/reference helpers.
The full electronic supplement is also supplied as `paper/evidence.zip`
at repository root, preserving the paths named in the standalone paper.

Run analyses through `tools/run_study.py` or the main/figure replay tools.
The sole manuscript and its PDF live in the root `paper/`. This study tree keeps
historical directory names for computational dependencies; it contains analysis
code and evidence, not additional manuscripts. Historical manuscript assembly
and desktop packaging are outside this subset.
