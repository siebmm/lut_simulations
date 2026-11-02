# Final chapter and SAFT appendix

- `chapter_comsol.tex`: the user's final COMSOL report chapter, unchanged.
- `appendix_saft_details.tex`: the detailed SAFT appendix, unchanged.
- `appendix_preview.pdf`: the existing compiled 67-page appendix.
- `appendix_preview.tex`: the standalone appendix wrapper.
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

Five chapter citation keys were not supplied in the original package:
`comsol-elastic-mesh`, `kaiser-1974`, `scipy-kaiserord`, `wan-2017`, and
`yu-2021`. These must come from the parent report's bibliography. No replacement
entries have been invented. The standalone appendix has all of its cited entries.

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

## Supporting research

The appendix describes an electronic evidence supplement. That supplement,
raw acquisition data, figure-generation scripts, study comparisons, and
verification reports belong to the full AllseasNDT study and are deliberately
not duplicated here. References to those records in the scientific text are
preserved. This folder contains the documents and their compilation assets.
