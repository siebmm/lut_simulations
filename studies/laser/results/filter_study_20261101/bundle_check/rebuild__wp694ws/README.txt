COMSOL CHAPTER: FULL AND PARTIAL LASER-WELD IMAGING
1 November 2026

INTEGRATION INTO THE SUPPLIED MAIN REPORT

Unzip comsol_chapter_latex.zip into the folder containing main.tex, preserving
figures/comsol/ and tables/comsol/. Review existing files before replacing them.
The supplied main.tex already contains the required line:

    \input{chapter_comsol}

The chapter has no document class, standalone title page or bibliography.
It uses the report's chapter, equation, table and figure numbering. All labels
and BibTeX keys begin with comsol: to avoid collisions with other chapters.

Heading spacing is defined once in comsol_heading_style.tex. Each COMSOL
chapter loads this style locally, including the chapter title and all section
and paragraph heading levels, paragraph gaps and surrounding float spacing.
Headings have more space before them than after them; chapter titles use a
normal block layout without the original negative top-spacing adjustment.
Keep this file with the other TeX inputs.
No change to main.tex is required for the corrected COMSOL heading spacing;
the surrounding report chapters retain their existing heading style.

Two existing preamble conflicts must be corrected:

1. Load xcolor with its table option immediately after \documentclass, before
   svg or any package that may load xcolor indirectly:

       \usepackage[table]{xcolor}

   Remove xcolor from the later grouped \usepackage list and remove its other
   standalone load. Loading the table option only after svg can still clash.

2. Replace all three caption loads with one line, before subcaption:

       \usepackage[font=footnotesize,skip=0.333\baselineskip]{caption}

After chapter_conclusion, inside the existing 1.25 spacing environment, add:

    \appendix
    \input{appendix_comsol_full}
    \input{appendix_comsol_partial}

If the report already has appendices, call \appendix only once and insert these
two inputs among them. References to appendix letters update automatically.

Replace the bibliography line with:

    \bibliography{bibfile,comsol_references}

Alternatively merge the six entries from comsol_references.bib into bibfile.bib
and retain the report's original bibliography command.

main_integration_example.tex shows these exact modifications to the supplied
main.tex. It is an integration example: the other chapters, title-page logos
and original bibfile were not provided. They have not been fabricated.

PREVIEW AND COMPILATION

preview.tex reproduces the supplied report's 11 pt report class, A4 margins
and 1.25 body spacing, with consistent heading gaps and package conflicts fixed.
It contains the chapter, both appendices and the bibliography.

The appendices are complete, self-contained simulation papers, not abbreviated
workflow notes. Their model, source calibration, equations, acquisition,
controls, complete path results, sensitivity registers and conclusions are
included in full. The same content is used in both the report and the separate
standalone builds, so they cannot silently drift apart.

    tectonic --keep-logs preview.tex
    tectonic --keep-logs full_standalone.tex
    tectonic --keep-logs partial_standalone.tex

The canonical delivery PDFs are comsol_chapter_and_appendices.pdf,
full_weld_simulation_paper.pdf and partial_weld_simulation_paper.pdf in the
repository's studies/laser/results/reports folder. The previous shorter
manuscript is preserved in archive/comsol_chapter_before_expansion_20261101.zip.

Or run pdflatex, bibtex, pdflatex and pdflatex on preview. No SVG conversion,
shell escape or Python is required to compile the delivered PDF figures.
The companion preview is a chapter/appendix integration check, not a build
of the unavailable remainder of the main report.

SCIENTIFIC SCOPE

The main chapter selects supported examples, retaining unsuccessful paths and
processing variants in the appendices and numerical evidence. The full-weld
direct-S example uses a 257-tap 0-5 MHz filter with a conservative 3.28 us guard.
The partial-weld S|SS example uses 657 taps and 0.5-5 MHz. Both studies use
equal receiver averaging, with no modeled receiving/wall response correction. These are exploratory findings for a known simulated flaw, not
independently validated optimal settings or experimental detection limits.

The LOF contour is drawn explicitly only in the geometry illustration. In SAFT
images a small tick outside the upper edge shows the mean LOF x coordinate;
an external grey bracket shows the tip-depth interval. No bright line or large
tip circles obscure the image. The bracket is a location guide, not measured
flaw length. The comparison of marker alternatives is in the evidence package.

NUMERICAL EVIDENCE AND REGENERATION

comsol_chapter_evidence.zip is the larger, optional research package. It adds:
384 configurations / 1152 signal-case images, full metrics and fixed protocols,
input fingerprints, calibration data, geometry-extension summaries, figure
generation sources, verification records and a snapshot of the research driver.
Raw acquisition CSVs and COMSOL solutions are not duplicated.

From the unpacked evidence directory, with NumPy, SciPy and Matplotlib:

    python make_tables.py
    python geometry_artwork.py
    python make_artifacts.py
    python source_calibration.py
    python expand_reports.py
    python recording_artwork.py
    python explanation_artwork.py
    python verify_data.py

The geometry/source-response artwork has a portable generator and vector assets.
The added calibration plots, complete control atlases and all numerical tables
are generated from the same verified evidence. Recording-level figures are
regenerated from compact plotting arrays in data/recordings/, with CSV input
fingerprints and checks against the original arrival times, amplitude references
and native filtered trace. recording_artwork.py --collect refreshes these arrays
from repository acquisition files. data/recording_figure_provenance/ identifies
the historical assets, before the legend-layout revision; it does not describe
the current figure hashes. The historical expand_reports.py --assemble command
is retired so it cannot overwrite the revised manuscripts. Default regeneration
updates artwork and numerical tables; the .tex manuscripts are edited directly.
The numerical research runner is intended to run in the NDT repository with
its verified saft_data folders and C++/Qt toolchain; see research/README.txt.
The standard workbench frequency controls do not encode the shorter filter's
tap count and separate timing guard. Existing workbench defaults are unchanged.

Figures are vector PDFs with 600 dpi PNG companions. Numerical scores always
use untruncated arrays; display range, annotations and colour do not affect them.
Legends occupy dedicated space outside the plots. figure_layout.py checks their
rendered bounds against the canvas, data panels and axis text during generation;
verification/figure_layout/ stores the results. This check complements visual
inspection of every figure and of the compiled chapter/standalone papers.
With PyMuPDF installed, python verify_figures.py also checks every figure PDF
for text outside the exported canvas and missing embedded fonts.

READING STRUCTURE

The chapter develops the physical measurement, independent source calibration,
a worked SAFT example, route/coupling/time selection, receiver/source combination,
and image interpretation before the results. The results distinguish common-rule
path comparisons, selected processing and acquisition-group reductions.

Each standalone paper remains complete and self-contained. Its principal findings
and discussion precede a complete-evidence section containing the full atlases,
filter definitions, timing details and configuration registers. The same content
is used in the chapter appendices. signal_cases.tex, saft_walkthrough.tex,
image_reading.tex and processing_controls.tex supply consistent explanations.

The two new figures in explanation_artwork.py are native vector artwork. The
worked example samples the existing S0/R94 waveform at the explicitly calculated
S|S time for an illustrative pixel (0,20) mm; it uses an existing common-protocol
image and performs no new reconstruction. verification/walkthrough_check.json
records the geometry, delay and sampled complex value. The primary study now
uses equal receiver means throughout. All 280 partial-weld configurations have
been reconstructed with that declared rule, and all resulting tables and figures
use those arrays. The 104 full-weld configurations retain their equal means.

EQUAL RECEIVER RECONSTRUCTION

The active partial-weld driver is paper/equal_receiver_study_20261101/run.py;
its output is results/equal_receiver_study_20261101. The full-weld results remain
in results/filter_study_20261101. make_artifacts.py --collect merges these two
explicit sources and rejects any publication configuration with non-unit
receiver weights. The data and narrative do not mix processing modes.
The supplementary bevel/corner outcomes use the equal-mean 0.2 mm tests;
the six conventional supplementary controls use 0.1 mm.

GATE AUDIT AND SURFACE TERMINOLOGY

paper/gate_ablation_20261101/run.py compares wall_gate=0.05 versus zero and
min_offset_mm=1.5 versus zero, independently and jointly, across 98 declared
configurations. The separate audit supplies 294 ablated configurations / 882
paired signal-case comparisons. Both gates remain in the common study protocol:
the wall gate is inactive for the principal images but affects multiple-reflection
partial-weld paths; the separation rule measurably changes both principal images.
data/gate_audit/ supplies its fixed protocol, all paired metrics and verification.
The research/ snapshot includes its reproduction and reporting drivers. Reconstructed
audit arrays are kept in results/gate_ablation_20261101 in the repository; they are
not duplicated in the portable evidence archive.

Weld surface means the exposed top of the weld fill, at depth zero for the full
weld and at the 16.4 mm groove floor for the partial weld. Side surfaces means
the two flat top surfaces beside the weld at depth zero. The manuscript, tables
and figure legends use these terms consistently. Historical machine identifiers
plate/floor remain unchanged in the evidence so that existing files stay traceable.
