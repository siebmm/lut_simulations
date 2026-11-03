# Final chapter and SAFT appendix

- `chapter_comsol.tex`: the user's final COMSOL report chapter, unchanged.
- `appendix_saft_details.tex`: the detailed SAFT appendix, unchanged.
- `appendix_preview.pdf`: the existing compiled 67-page appendix.
- `appendix_preview.tex`: the standalone appendix wrapper.
- `standalone.tex`: a complete local driver for the final chapter and appendix.
- `standalone.pdf`: the compiled 83-page chapter and appendix, with all citations resolved.
- `figures/`: all figures referenced by the chapter and appendix.
- `chapter_references.bib` and `appendix_references.bib`: cited bibliography
  entries available in the original document packages. The chapter file restores
  the original `levesque` key from its documented `levesque-2016` remapping.

## Integrate into the report

Copy the two chapter sources into the directory containing your report's
`main.tex` and merge `figures/` into its figures directory. Include the chapter
in the main text and the detailed appendix after `\appendix`:

```latex
\input{chapter_comsol}
% ... remaining report chapters ...
\appendix % once, before the first appendix
\begin{spacing}{1.25}
\input{appendix_saft_details}
\end{spacing}
```

Place this appendix first to match the chapter's Appendix A references. Merge
only missing bibliography keys into the report's existing bibliography.
The report supplies its document class, packages, main bibliography, and any
other chapters/appendices referenced by these sources. The existing worked
SAFT example is a separate report appendix.

All cited bibliography entries are included. The five entries absent from the
initial companion repository were copied from the user's supplied `bibfile.bib`;
their provenance is recorded in `../verification/bibliography_import.json`.
`chapter_references.bib` and `appendix_references.bib` have disjoint keys and
can be used together, as in `standalone.tex`.

## Compile the appendix alone

From this directory, with a suitable TeX installation:

```sh
tectonic appendix_preview.tex
# or
latexmk -xelatex appendix_preview.tex
```

The source wrapper contains the required packages. Its existing PDF is retained
without changes; the final chapter is an inclusion source, not a standalone
document.

For the chapter and appendix together, run `tectonic standalone.tex` or
`latexmk -xelatex standalone.tex` from this folder. The final manuscript text
is preserved, including its literal `XX` pointer to the user's separately
written worked-example appendix; that editorial pointer is not a build input.

## Supporting research

`evidence.zip` is the complete electronic supplement described in the appendix.
Extract it here to obtain the named `evidence/` hierarchy, including all 20,412
partial-weld atlas rows, failed/negative outcomes, source subsets and settings.
The archive is compressed to avoid hundreds of repeated settings records in
the working tree.

Both raw SAFT acquisitions are included in `../data/`; scientific scripts,
compact figure arrays and frozen comparison records are under
`../studies/laser/`. From the repository root, `python tools/reproduce_main.py`
checks the main results and `python tools/reproduce_figures.py` regenerates the
artwork in `output/workspace/` without overwriting these final PDFs or sources.
