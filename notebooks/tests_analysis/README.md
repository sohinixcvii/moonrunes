# tests_analysis — the analyses and tests around the pipeline, one notebook per Part

The analyses whose numbers or caches the paper uses, and the tests run on them: the beam-matched ring TODs, the
spectral fits along the TRIS ring, the model comparison, the 408 MHz sky models, the T–T indices, the Haslam
scale/offset fits, the zero-level and survey-subset robustness tests, the prediction at other frequencies, the beam
tests and the Fisher information of the 408 MHz prediction. `10_paper_plots.ipynb` reads the caches and CSVs written
here.

**The pipeline comes first.** `python run_pipeline.py --step all` (`moonrunes.pipeline`, `docs/PIPELINE.md`) makes the
beam-matched TODs, the model comparison (TRIS excluded), the model selection, our spectral model (the selected model
fitted with TRIS), the 408 MHz and the RHINO predictions in `paper_outputs/pipeline/`. The Parts read those products
and do not rebuild them; a Part's own results are checks (`test_outputs/`) unless the paper quotes them.

## The Parts

| Notebook | What | Writes (main files) | Feeds the paper |
|---|---|---|---|
| `Part_A_Beam_matched_ring_TODs` | step 2's beam-matched ring TODs against TODs observed at each map's native resolution | `test_outputs/beam_matched_vs_old/` | — |
| `Part_B_Spectral_fits_along_the_ring_and_the_408_MHz_sky_models` | curved power-law fits per ring sample; the 408 MHz sky models against the power law (step 6) and the curved power law | `test_outputs/sky_models_408/` | — (a check of steps 5–6) |
| `Part_C_Model_comparison_diagnostics` | diagnostics of the step 3 model comparison | `paper_outputs/model_comparison/diagnostic_plots/` | Section 5 counts |
| `Part_D_Model_choice_decision_sheet` | model-choice decision sheet (five criteria) | `paper_outputs/model_choice/` | stacked residuals |
| `Part_E_T_T_indices_and_the_sky_models_against_TRIS` | T–T indices along the ring; sky models fitted to TRIS at 600/820 MHz | `paper_outputs/phil_plots/` | `tab:step1` |
| `Part_F_Haslam_and_the_sky_models_against_the_power_law_prediction` | figures of step 6's scale–offset fits (York) | `test_outputs/plot3_haslam_vs_fit/` | Haslam/prediction ratios |
| `Part_G_Haslam_calibrated_with_one_gain_and_offset` | Haslam calibrated with one gain and offset | `test_outputs/plot3_haslam_calibration/` | — |
| `Part_H_Zero_level_offsets_as_fit_parameters` | one zero-level offset per survey as a fit parameter | `paper_outputs/zero_level_offsets/` | zero levels; input of I |
| `Part_I_Robustness_of_the_Haslam_gain_and_offset_survey_subsets_and` | robustness of the Haslam gain and offset (survey subsets, sample subsets, zero-level treatment) | `paper_outputs/robustness_check/` | g, z robustness |
| `Part_J_Fits_on_13_points_and_the_TOD_prediction_at_any_frequency` | least-squares fits on 13 points and the TOD prediction at any frequency (Tests 1–5); β against Irfan & Puglisi's maps | `test_outputs/tod_prediction/` | — |
| `Part_K1_Beam_widths_FWHM_and_the_smoothing` | beam widths (FWHM) and the smoothing code | — | the 19.155°/23.366° beam |
| `Part_K2_Beam_matching_harmonic_checks_and_beam_specification` | elliptical beam on the map route, harmonic checks, convolve→deconvolve, two-step vs one-step matching, beam variants | `test_outputs/beam_variants/` | beam paragraphs |
| `Part_K3_…`, `Part_K4_…`, `Part_K5_…` | saved results of the script-based beam tests (`scripts/beam_averaging_blobs/`, `scripts/fwhm_scaling_test.py`, `scripts/beam_tests_lwa35_stockert1420/`) | — | beam paragraphs |
| `Part_L_Fisher_information_in_the_408_MHz_prediction` | Fisher information on ln A(408) and σ(ln A) with data sets left out; g at the median index | `paper_outputs/fisher_share/` | Fisher shares |

**Part K**: how the TRIS beam is described and applied, and how much it matters. K1 checks the beam widths and the
smoothing code; K2 is the beam-matching and harmonic test suite; K3–K5 show the script-based tests behind the beam
paragraphs of the paper from their saved results.

## Data sets

Parts B, D, F, G, H, I and L read step 5's data sets: every map `sed.exclude_maps` leaves in with data on the ring,
plus TRIS 600 and 820 — 13 points with the current config (`moonrunes.pipeline.fit_data_sets`). Part I also fits the
six-survey subset (TRIS 600, 820, LWA1 35, Maipu 45, LWA1 80, Stockert 1420) as one of its survey subsets; Part L
gives its Fisher numbers for that subset too. Part C and criterion 1 of Part D use step 3's data sets (TRIS excluded,
11 points). Part E's T–T fits use step 5's data sets plus Haslam (14 TODs, 90 pairs). Part J lists its own 13 points.

## Running

* **Kernel:** the `limtod` conda environment (with `joblib`). Any working directory inside the repo: each Part finds
  the repo root from `configs/run_pipeline.yaml`.
* **Order:** the pipeline first; then A; B before D–I and K2 (B's curved fits); H before I. Each Part starts with its
  own setup cell and runs on its own once the Parts or pipeline steps it reads have run.
* **Run time:** the first code cell records the start time and the last cell prints the total run time.
* **Parallel fits:** the UltraNest runs of B, C and D and the maximum-likelihood fits of C run in parallel (joblib,
  one process per core). Each run seeds numpy's random generator itself, or its random starting points are drawn
  in the serial order first, so the results are the ones a serial run gives. Part C's `run_nested` is unseeded.
* **Pipeline products, not copies:** Parts B, C and J take their SED points, Haslam's ring TOD and the TRIS rings from
  step 2's `ring_tods.npz` (`moonrunes.pipeline.load_ring_tods`, `points`). Part F reads step 6's
  `scale_offset_408.csv` and Part B the power law's rows of `residuals_408.csv`.
* **Caches:** a cached fit is reused only if it was made from the same ring TODs and on the same data sets
  (Part B's curved fits against step 5's data sets; Part D's runs record their data sets).
