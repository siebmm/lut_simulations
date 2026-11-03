LASER-ULTRASONIC WELD IMAGING: PAPER AND NUMERICAL EVIDENCE
30 October 2026

Read paper.pdf first. It is a concise synthesis of the full/partial comparison
and a new, explicitly exploratory partial-weld improvement study.

MAIN FOLLOW-UP FINDING
For partial-weld S | SS at 0.5-5 MHz, change Receiver sum from Timing-only mean
to Known-response weighted. Keep the other physical controls unchanged.
The near-flaw / strongest distant background contrast increases from 4.36 to
8.05 dB, with 95.8% near-flaw support in both cases. This is an improvement in
relative contrast on one simulated flaw, not uniform background suppression,
better sizing, or validated experimental sensitivity.

LOAD IN THE C++ WORKBENCH
1. Open studies/laser/saft_data/partial_weld (or your relocated copy).
2. Load workbench/saft_settings_partial_S_SS_weighted.json and reconstruct.
   This is the main broad-band candidate with all source and receiver groups.
3. Optional: saft_settings_partial_S_SS_weighted_plate.json keeps all source
   heights but uses only the two plate receiver groups.
4. Optional: saft_settings_partial_S_SS_band1_weighted.json uses 1-5 MHz.
5. Open saft_partial_response_comparison.saftblend in Combine paths to compare
   four precomputed alternatives. Only broad-band weighting is enabled by
   default: the layers are alternative reconstructions of the same path,
   so adding them together is not independent multiview evidence.

For the full weld, open the full_weld dataset first. Load
saft_settings_full_shoulders_plates.json for direct S with two shoulder sources
and the 314 plate receivers. The four-path saft_full_shoulders_plates.saftblend
contains the exact corresponding direct and reflected core views.
All new files are distinct from the user's original settings and blend files.

DISPLAY
Image panels retain the normal magma heatmap on black for low signal.
Known LOF is a thin cyan dashed segment with open tip markers. Grey is missing
support, not zero amplitude; white in geometry overviews is air.
The common physical display reference is 80 micrometres/second. The core atlas
uses 60 dB range; partial improvement and wider-field panels use 40 dB.
Metrics are always calculated from untruncated float64 arrays.
The two blend files use pastel display colours, not statistical confidence.

FILES AND REPRODUCTION
paper.tex, references.tex, table_*.tex: editable paper sources.
figures/: seven vector PDFs and 600 dpi PNGs; no GUI screenshots.
  F01 shows the complete specimen width, the actual saved full-weld box
  (x +/-6 mm, depth 5-25.25 mm) and the partial-weld study box separately.
  Its J-root outline is closed, with equal coordinate scales in detail (c).
  F07_whole_specimen is a standalone full-width overview at equal x/z scales.
  The paper's numerical full/partial comparisons still use the matched lower
  window specified in its methods; changing the geometry annotation does not
  change or relabel those reconstruction arrays.
data/index.json: 138 archived reconstruction configurations, 270 case images,
  settings, metrics, stage labels and hashes of the compressed numeric arrays.
  This includes all 82 new configurations / 142 new case images, 36 prior
  baseline configurations / 108 cases, and 20 prior verification controls.
data/*.npz: x_mm/depth_mm coordinates in mm, points_m in metres; LOF, pristine
  reference and signed-waveform-difference amplitude scores in m/s. Only the
  processed variant (native variant 2) is retained. Per-source complex means
  and validity arrays permit independent RMS checks. NaN means unsupported.
all_results.csv: one row per archived signal-case result.
followup_protocol.json: initial follow-up plan, preserved unchanged.
confirmation_protocol.json: subsequent plan after response weighting showed
  a promising gain. This is deliberately not represented as blind validation.
data/provenance.json: input fingerprints and native executable/code hashes.
data/source_calibration/: compact, unfiltered exports from the existing 10 and
  100 J/m2 thermal/dipole calibration runs, with original input SHA256 hashes.
  Source amplitudes, fitted and withheld-angle errors, temperature, thermal
  energy and displacement/velocity checks are reproduced independently.
  Array units and the production parameter values are in provenance.json.
requests/: complete native CLI requests for each archived configuration,
  plus request-index-to-array mappings. Raw FE records are not included.
workbench/: five settings files, two native blends and exact Qt readback checks.
verification/: final PDF checks and independent source-bundle rebuild result.

Rebuild the PDF from this directory with Tectonic or another modern LaTeX engine:
  tectonic --keep-logs paper.tex
Rebuild figures and tables with Python (NumPy and Matplotlib):
  python figures.py
  python tables.py
Regenerate workbench exports with Python and PySide6:
  python export_workbench.py
These three operations read only the portable data in this bundle.
Verify the source-calibration values in Section 2.2 with NumPy and SciPy:
  python source_calibration.py
This also uses only portable data; --archive is for repository maintainers
collecting the original raw exports, not for a standalone bundle rebuild.

To reconstruct scores from the original raw records, use the native executable
identified in data/provenance.json and the matching dataset folders:
  SAFTWorkbenchCpp.exe --dataset PATH_TO_partial_weld --compute requests/partial.json --output OUTPUT_DIRECTORY
  SAFTWorkbenchCpp.exe --dataset PATH_TO_full_weld --compute requests/full.json --output OUTPUT_DIRECTORY
Each command verifies its dataset inputs. Request indices are listed in the
corresponding mapping JSON. Eight worker threads are requested. Existing raw
data and model files must be supplied separately; they are not copied here.
The repository-oriented analyse.py / confirm.py / evidence.py record how this
follow-up was prepared; they additionally depend on the preceding study folder.

INTERPRETATION
The common primary geometry-led controls were established before the earlier
18-path comparison. This follow-up knows that comparison and the prescribed
flaw. The small set of factors is physically motivated, but the preferred
weighting result is selected with knowledge of its outcome. Fine/wide checks
are consistency checks on the same case, not held-out flaw validation.
Source/receiver group reductions establish useful combinations for this model.
They do not establish a minimum receiver count or a universally required side.
No new FE solves, mesh-convergence claims, instrument-noise simulations, or
modifications to the application or raw records are part of this report.
