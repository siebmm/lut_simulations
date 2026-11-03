# Laser-ultrasonic imaging of full and partial welds

**[Read the complete paper (PDF)](paper.pdf).**

This is the repository's only manuscript. It contains the full detailed methods,
source calibration, SAFT definitions, results, negative comparisons, path and
source-subset registers, and reproduction instructions. It is a self-contained
article with an abstract, numbered sections, conclusions and its own references.

- `paper.tex`: complete editable LaTeX source.
- `paper.pdf`: rendered paper, ready to read.
- `references.bib`: the three cited bibliography entries.
- `figures/`: exactly the 27 figures cited by the paper.
- `evidence.zip`: numerical supplement, including all 20,412 acquisition/path
  rows and the settings, controls and negative outcomes described in the paper.

## Build the PDF

Install Tectonic, then run from the repository root:

```sh
python tools/build_paper.py
```

Use `--tectonic PATH` if the executable is not on `PATH`. Tectonic downloads its
standard TeX resources on first use. The build needs no COMSOL, native SAFT tools,
or Python scientific packages. It writes the finished PDF here, intermediate
files to `build/paper/`, and a source/PDF hash receipt to `verification/paper.json`.
Missing references and overfull text boxes fail the build. After an intentional
edit and rebuild, update the repository snapshot with
`python tools/verify_repository.py --write-manifest`.

## Reproduce the evidence

The model inputs are in `../models/`, velocity histories in `../data/`, and
scientific analysis and figure generators in `../studies/laser/`. From the
repository root, `python tools/reproduce_main.py` checks the main numerical
results and `python tools/reproduce_figures.py` regenerates the artwork under
`output/replay/`. These workflows preserve the supplied manuscript assets.
Extract `evidence.zip` here to obtain the `evidence/` paths cited in the paper.
