# Calculations supporting the weld-imaging paper

These calculations compare how filtering, wave paths, receiver coverage and
source selection affect the reconstructed defect. They supply the reference
images and figures used in the [full paper](../../paper/paper.pdf).

| Scientific purpose | Location |
|---|---|
| Compare 384 baseline imaging configurations and check laser calibration | `paper/comsol_chapter_20261101/` |
| Explain ray coverage, source count and reflected paths | `paper/explanatory_saft_20261101/`, `paper/combined_explanatory_saft_20261101/` |
| Apply the published filters and image-quality measures | `paper/filter_study_20261101/study.py` |
| Reconstruct and compare the final partial-weld images | `paper/partial_improvement_20261102/` |
| Independently check travel times, boundary response and focusing | `src/` |

The dates identify the original calculations. These are supporting calculations,
not separate papers. The complete set of 20,412 partial-weld comparisons,
including unsuccessful configurations, is in `../../paper/evidence.zip`.

## Repeat the comparisons

From the repository root:

```sh
python tools/reproduce_main.py
python tools/reproduce_figures.py
```

New results go under `output/`. The first command checks the principal images;
the second regenerates the baseline, explanatory and final-review figures.
For a different specimen or new imaging settings, use
[the general analysis workflow](../../saft/FORMAT.md).

`replay_inputs.json` lists the numerical inputs required for these calculations.
Each is stored once, either here or in the paper supplement. The replay commands
assemble a working copy automatically. See [running the calculations](../../tools/README.md)
for setup and optional checks.
