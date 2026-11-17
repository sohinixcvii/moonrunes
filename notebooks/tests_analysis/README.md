# tests_analysis — the analyses and tests around the pipeline, one notebook per Part

Every analysis whose numbers or caches the paper uses, and the tests run on them: the beam-matched ring TODs, the
spectral fits along the TRIS ring, the model comparison, the 408 MHz sky models, the T–T indices, the Haslam
scale/offset fits, the zero-level and survey-subset robustness tests, the prediction at other frequencies, the beam
tests and the Fisher information of the 408 MHz prediction. `10_paper_plots.ipynb` reads the caches and CSVs written
here.

The Parts were split out of `notebooks/test_results.ipynb` (2026-11-17), cells copied as they were, with their stored
outputs; Part L is new (ported from the old notebook 08). `test_results.ipynb` stays in `notebooks/` as the record of
the merged notebook and is no longer maintained.

**The pipeline comes first.** `python run_pipeline.py --step all` (`moonrunes.pipeline`, `docs/PIPELINE.md`) makes the
beam-matched TODs, the model comparison (TRIS excluded), the model selection, our spectral model (the best model fitted
with TRIS), the 408 MHz and the RHINO predictions in `paper_outputs/pipeline/`. The Parts read those products and do
not rebuild them; a Part's own results are checks (`test_outputs/`) unless the paper quotes them.

## The Parts

"Cell offset" converts a cell number of `test_results.ipynb` (as quoted in `HANDOVER.md` and older notes) to the Part
notebook: cell in the Part = cell in `test_results` − offset.

| Notebook | What | Was | Cell offset | Writes (main files) | Feeds the paper |
|---|---|---|---|---|---|
| `Part_A_Beam_matched_ring_TODs` | beam-matched ring TODs (pipeline step 2) against the old TODs | `18_beam_matched_tods` | 1 | `test_outputs/beam_matched_vs_old/` | — |
| `Part_B_Spectral_fits_along_the_ring_and_the_408_MHz_sky_models` | curved power-law fits per ring sample; the 408 MHz sky models (the power law is step 5) | `12_sky_models_vs_fit_408` | 13 | `test_outputs/sky_models_408/` | — (a check of steps 5–6) |
| `Part_C_Model_comparison_diagnostics` | diagnostics of the step 3 model comparison | `11_model_comparison_diagnostics` | 31 | `paper_outputs/model_comparison/diagnostic_plots/` | Section 5 counts (main.tex:501) |
| `Part_D_Model_choice_decision_sheet` | model-choice decision sheet (five criteria) | `15_model_choice_decision_sheet` | 93 | `paper_outputs/model_choice/` | stacked residuals (main.tex:503) |
| `Part_E_T_T_indices_and_the_sky_models_against_TRIS` | T–T indices along the ring; sky models fitted to TRIS at 600/820 MHz | `14_ring_tt_betas` | 122 | `paper_outputs/phil_plots/` | `tab:step1` (Plot 6a) |
| `Part_F_Haslam_and_the_sky_models_against_the_power_law_prediction` | scale–offset of Haslam and the sky models against the prediction (York; a check of step 6) | `16_haslam_vs_power_law_408` | 160 | `test_outputs/plot3_haslam_vs_fit/` | main.tex:647 ratio |
| `Part_G_Haslam_calibrated_with_one_gain_and_offset` | Haslam calibrated with one gain and offset | `17_haslam_gain_offset_calibration` | 176 | `test_outputs/plot3_haslam_calibration/` | possibly main.tex:647 |
| `Part_H_Zero_level_offsets_as_fit_parameters` | zero-level offsets as fit parameters | `20_zero_level_offsets` | 194 | `paper_outputs/zero_level_offsets/` | main.tex:403; input of I |
| `Part_I_Robustness_of_the_Haslam_gain_and_offset_survey_subsets_and` | robustness of the Haslam gain and offset | `21_robustness_check` | 219 | `paper_outputs/robustness_check/` | main.tex:396 and others |
| `Part_J_Fits_on_13_points_and_the_TOD_prediction_at_any_frequency` | fits on 13 points and the TOD prediction at any frequency (Tests 1–5) | `19_power_law_tod_prediction` | 232 | `test_outputs/tod_prediction/` | — |
| `Part_K1_Beam_widths_FWHM_and_the_smoothing` | beam widths (FWHM) and smoothing | `09_beam_fwhm_tests` | 257 | — | the 19.155°/23.366° beam |
| `Part_K2_Beam_matching_harmonic_checks_and_beam_specification` | map-route elliptical beam, harmonic checks, convolve→deconvolve, two-step vs one-step matching, beam variants | `03_diagnostics` cells 39–114 | 277 | `test_outputs/beam_variants/` | main.tex:267 |
| `Part_K3_…`, `Part_K4_…`, `Part_K5_…` | saved results of the script-based beam tests (`scripts/beam_averaging_blobs/`, `scripts/fwhm_scaling_test.py`, `scripts/beam_tests_lwa35_stockert1420/`) | scripts | 355, 357, 359 | — | main.tex:273–275 |
| `Part_L_Fisher_information_in_the_408_MHz_prediction` | Fisher information on ln A(408) and σ(ln A) with data sets left out; g at the median index | notebook 08 cell 54 (new) | — | `paper_outputs/fisher_share/` | main.tex:625–626 |

**Part K** (the banner of the beam tests in `test_results`): how the TRIS beam is described and applied, and how much
it matters. K1 checks the beam widths and the smoothing code; K2 is the beam-matching and harmonic test suite from the
old diagnostics notebook; K3–K5 are the script-based tests behind the beam paragraphs of the paper (main.tex:271–275),
shown from their saved results.

## Data sets

Parts B, D, F, G, H, I and L read step 5's data sets: every map `sed.exclude_maps` leaves in with data on the ring,
plus TRIS 600 and 820 — 13 points with the current config (`moonrunes.pipeline.fit_data_sets`). Part I also fits
notebook 08's six ring surveys as a named subset and checks them against notebook 08's numbers; Part L does the same
check for the Fisher numbers. Part C and criterion 1 of Part D use step 3's data sets (TRIS excluded, 11 points).
Part E keeps its own T–T set (the six surveys plus Haslam). Part J chooses its own 13 points.

## Running

* **Kernel:** the `limtod` conda environment (with `joblib`). Any working directory inside the repo: each Part finds
  the repo root from `configs/run_pipeline.yaml`.
* **Order:** the pipeline first; then A; B before D–I and K2 (B's curved fits); H before I. Each Part starts with its
  own setup cell and runs on its own once the Parts or pipeline steps it reads have run.
* **Run time:** the first code cell records the start time and the last cell prints the total run time.
* **Parallel fits:** the UltraNest runs of B, C and D and the maximum-likelihood fits of C run in parallel (joblib,
  one process per core). Each run seeds numpy's random generator itself, or its random starting points are drawn
  in the serial order first, so the results are the ones a serial run gives (checked bit for bit on the real models).
  Part C's `run_nested` is unseeded, as it always was.
* **Caches:** a cached fit is reused only if it was made from the same ring TODs and on the same data sets
  (Part B's curved fits against step 5's data sets; Part D's runs record their data sets).
* **Stored outputs** are those of `test_results.ipynb` (not re-run), except Part L's.
