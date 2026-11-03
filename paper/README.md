# Laser-ultrasonic imaging of full and partial welds

**[Read the complete paper (PDF)](paper.pdf).**

The paper follows laser generation, source calibration, propagation through the
welds and SAFT imaging. It includes the detailed equations, assumptions,
comparisons and negative results, together with the full source/path registers.
It is the repository's only manuscript.

| Material | Contents |
|---|---|
| `paper.pdf` | Complete rendered paper |
| `paper.tex` and `references.bib` | Editable text and references |
| `figures/` | The 27 figures cited in the paper |
| `evidence.zip` | Numerical supplement, including all 20,412 acquisition/path comparisons |

To explore the work, start with the
[laser-generation model](../models/laser_generation/README.md),
[specimen models](../models/README.md), [recordings](../data/README.md), or
[imaging method](../saft/README.md).

To repeat the principal images and figures, follow
[the reproduction instructions](../tools/README.md#reproduce-the-paper).
Extract `evidence.zip` here to obtain the `evidence/` paths cited in the paper.
To compile the PDF from its source, install Tectonic and run
`python tools/build_paper.py` from the repository root.
