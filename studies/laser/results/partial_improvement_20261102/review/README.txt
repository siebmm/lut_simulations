COMPREHENSIVE PARTIAL-WELD REVIEW - 2 NOVEMBER 2026

CURRENT WORKING CHOICE
All five sources/all 331 receivers, S | SS, vy, complex response weighting,
1-5 MHz, 657 taps, beta 5.10226, 15% angular gate, 5% wall gate, 1.5 mm exclusion,
eight real/effective receivers, 25% eligible retention, two qualifying sources,
CF exponent one, RMS, direct L/S/R pulse screens, and only the specimen-only
backwall-SS fixed mute with 0.30 microsecond half-width. Its central +/-0.225 us
is zero; it rises smoothly to one at +/-0.30 us. Full FIR support remains required.

At 0.025 mm: Cmax 14.44694 dB, near support 96.90173%, near peak 58.63628 um/s,
background 11.11270 um/s, peak distance 0.07111 mm. Matched unmuted Cmax 9.08623 dB.
Aligned 0.05 mm gives 15.76544 dB; shifted 0.05 mm gives 14.45208 dB.

PRESETS (schema 4, original partial-weld dataset fingerprint)
settings/saft_partial_weighted_SS_mute_025.json - 207,345-pixel publication grid
settings/saft_partial_weighted_SS_mute_05.json - quicker aligned 0.05 mm view
settings/saft_partial_right_S_LL_mute_05.json - S3/right single-shot alternative
Selections derive from native-confirmed runs; live GUI loading was not tested.

EVIDENCE MAP
protocol.json, atlas_metrics.csv, atlas_summary.json: 972 native group exports
combine into 20,412 fixed acquisition/path views, with three signal controls each.
verification.json: 16 independent native requests / 48 signal-case checks.
settings/: 252 native configurations, offline comparisons, protocols, metrics.
confirmation/: 36 native finer/wider/mute-width/receiver-removal checks.
refinement/: muted/unmuted 0.025 mm and shifted 0.05 mm; all 80 source-subset views.
path_confirmation/: 15 fixed atlas acquisitions, each fine and wide (30 native).
final_checks/: changed offline controls verified natively; nested-grid checks;
0.1 mm wider-field result; exported settings hashes.
publication_atlas_summary.json, publication_inputs.json: paper/register provenance.
pdf_qa/: independent Poppler page renders, font/bounds/reference checks.

INTERPRETATION
The original 25% fraction and receiver-count floors are inactive in both original
all-aperture anchors (600 threshold combinations). Missing support is governed
by visibility, launch sensitivity and source qualification. Near support is not
fraction of flaw length detected; the earlier exact segment had two-source support
at all 101 sampled points under the nominal 20% angular gate.
S3/right S | LL with narrow LL mute gives 5.93 dB locally / 6.82 dB wider.
S3/left unmuted S | SS gives 11.49 / 11.03 dB. S1/left+weld SLS | S gives
15.80 / 4.01 dB. High-scoring pristine lookalikes and remote false indications
remain in the registers. Corner amplitudes/phases are not calibrated.
The principal wider 0.1 mm field has only 3.76% whole-field support: high local
support does not mean the entire plate is inspected. No noise, mesh/absorber
convergence, independent-flaw validation or experimental detection claim.

REPRODUCTION
Source scripts: studies/laser/paper/partial_improvement_20261102/
  build_review.ps1 builds the exporter against unchanged production native sources.
  review.py verify / atlas performs the catalogue and group checks.
  settings_review.py run performs settings and offline comparisons.
  review_confirm.py run, review_refine.py, review_paths_confirm.py run follow-ups.
  review_verify.py checks consistency and exports presets.
  review_publication.py generates only the two requested manuscripts.
  review_qa.py checks/renders compiled PDFs; review_delivery.py packages them.
Large arrays stay locally as arrays.npz. The compact evidence ZIP includes all
metrics, protocols, scripts, figures and settings but not source recordings.
One failed atlas batch was resumed from completed exports; final checks cover
all 972 completed configurations. Earlier study scripts/archives remain intact.
