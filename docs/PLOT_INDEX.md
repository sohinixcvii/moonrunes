# Plot index

Every plot the code writes, where it goes and what makes it. Cell numbers are 0-based (as `nbformat` counts them); the
section heading in the notebook is the stable way to find a plot if cells are added.

| Where | What is there | Made by |
|---|---|---|
| paper repository `figures/` | the paper figures `fig1_…`–`fig10_…` (PDF + PNG; the paper includes the PDFs) | `notebooks/10_paper_plots.ipynb` |
| `paper_outputs/pipeline/step*/` | the pipeline's products (CSV, NPZ; no figures) and `PIPELINE_REPORT.md` | `run_pipeline.py` |
| `paper_outputs/model_comparison/diagnostic_plots/` | diagnostics of the step 3 model comparison | Part C |
| `paper_outputs/model_choice/` | the model-choice decision sheet | Part D |
| `paper_outputs/phil_plots/` | T–T β; sky models against the TRIS rings | Part E |
| `paper_outputs/zero_level_offsets/` | zero-level offsets as fit parameters | Part H |
| `paper_outputs/robustness_check/` | robustness of the Haslam gain and offset | Part I |
| `paper_outputs/fisher_share/` | Fisher information on ln A(408) (CSVs) | Part L |
| `test_outputs/beam_matched_vs_old/` | step 2's TODs against native-resolution TODs | Part A |
| `test_outputs/sky_models_408/` | sky models at 408 MHz against the power-law and curved fits | Part B |
| `test_outputs/plot3_haslam_vs_fit/` | Haslam and the sky models in the scale–offset plane (step 6's fits) | Part F |
| `test_outputs/plot3_haslam_calibration/` | Haslam calibrated with one gain and offset | Part G |
| `test_outputs/tod_prediction/` | least-squares fits and predicted TODs (Tests 1–5) | Part J |
| `test_outputs/beam_variants/` | beam-specification variants | Part K2 |
| `test_outputs/fwhm_scaling_test/`, `beam_averaging_blobs/<map>/`, `lwa35_stockert/` | the script-based beam tests | `scripts/` (shown in Parts K3–K5) |
| inside the notebook only | input, frame and TRIS-map checks; beam and harmonic tests | `diagnostics.ipynb`; Parts K1, K2 |

"Parts" are the notebooks in `notebooks/tests_analysis/`. The TODs are always the pipeline's beam-matched ring TODs
(`tris_cuts` beam, 120 samples). Data sets: step 5's (13 points with TRIS) unless a row says otherwise; step 3's
comparison has 11 (TRIS excluded).

## Paper figures — `notebooks/10_paper_plots.ipynb`

| File | Shows | Label | Cell |
|---|---|---|---|
| `fig1_data_maps_tris_slice` | the twelve maps with data on the ring (step 5's survey maps and Haslam), native resolution, over the TRIS declination band, two columns | `fig:maps-slice` | 5 |
| `fig2_method_tris_beam` | the TRIS beam on the sky (dB, north up, east left) with the E/H axes, and the E- and H-plane cuts | `fig:tris-beam` | 20 |
| `fig3_method_tris_maps` | the TRIS 600 and 820 MHz maps from step 0 | `fig:tris-maps` | 7 |
| `fig4_method_models_vs_tris_tods` | the sky models and their corrections at 600.5 / 817.8 MHz against the TRIS rings; model − TRIS | `fig:tods-600-820` | 10 |
| `fig5_tt_beta_matrix` | T–T β ± σ and χ²/ν for every pair of the 14 ring TODs | `fig:beta-matrix-ring` | 14 |
| `fig6_tt_beta_vs_frequency` | T–T β against √(ν₁ν₂), bar from ν₁ to ν₂ | `fig:beta-vs-freq` | 16 |
| `fig7_modelcomp_fits_vs_tris` | maximum-likelihood power law, curved power law and moment expansion at 600.5 / 817.8 MHz against TRIS | `fig:modelcomp-tris` | 18 |
| `fig8_skymodel_models_vs_prediction_408` | our 408 MHz prediction (±1σ) against the 408 MHz TODs of Haslam and the sky models; model − prediction | `fig:models-vs-prediction-408` | 28 |
| `fig9_skymodel_tris_anchor_gz` | the g–z plane (Haslam = g P + z) with the 1σ ellipses for TRIS only, without TRIS and all data sets, and TRIS only against a fixed extrapolation β | `fig:tris-anchor` | 23 (fits in 22) |
| `fig10_rhino_tods_rhino_beam` | our prediction at 60–80 MHz through the RHINO beam against the sky models | `fig:rhino-rhino-beam` | 26 (function in 25) |
| `rhino_tods_tris_beam` | the same through the TRIS beam (not in the paper) | — | 25, last line |

## Part A — beam-matched ring TODs → `test_outputs/beam_matched_vs_old/`

| File | Shows | Cell |
|---|---|---|
| `fig1_beam_matched_minus_old` | beam-matched − native-resolution TOD [K] for the surveys and the 408 MHz models; table `beam_matched_minus_old.csv` | 10 |

## Part B — sky models against the fits at 408 MHz → `test_outputs/sky_models_408/`

| File | Shows | Cell |
|---|---|---|
| `fig1_T408_vs_ra` | the power-law and curved 408 MHz predictions ±1σ and the seven sky models along the ring | 10 |
| `fig2_residual_vs_ra` | sky model − fit [K], against each fit | 12 |
| `fig3_chi2_per_sample` | χ²/N per model and fit; table `chi2_sky_models_vs_fits.csv` (cell 14) | 16 |

## Part C — model comparison diagnostics → `paper_outputs/model_comparison/diagnostic_plots/`

| File | Shows | Cell |
|---|---|---|
| `fig1_dlogZ_vs_ra` | Δln Z of each model against the power law per ring sample, tie band, Cygnus A | 7 |
| `fig2_logZ_vs_ra` | ln Z of each model along the ring | 9 |
| `fig3_dlogZ_histograms` | histograms of Δln Z | 11 |
| `fig4_preference_counts` | samples per evidence category | 13 |
| `fig5_sed_fits_per_unit` | the fitted models with 68 % bands at six ring samples; Haslam plotted, not fitted | 18 |
| `fig6_maxlogl_vs_logZ` | fit gain (Δ max ln L) against evidence gain: the Occam penalty | 20 |
| `fig7_C_posterior` | posterior of the curvature C with its prior | 22 |
| `fig8_prior_sensitivity_C` | curved − power-law Δln Z against the width of the C prior | 26 |
| `fig9_convergence_dlogZ` | Δln Z from repeated and re-tuned UltraNest runs against the exact grid integral; `convergence_check1`–`4*.csv` (cells 30–36) | 38 |
| `fig10_models_at_tris_frequencies` | each model's maximum-likelihood fit at 600.5 / 817.8 MHz with the TRIS rings; `part5_*.csv` | 43 |
| `fig11_chi2_per_sample` | χ² and χ²/ν of every fit along the ring; `part5_chi2_*.csv` | 45 |
| `fig12_dscore_vs_ra_masked` | Figure 1 with the Cygnus A mask; `part6_counts_masked.csv` | 49 |
| `fig13_occam_penalty` | how much of each Δln Z is the Occam penalty; `part6_occam_penalty*.csv` | 51 |
| `fig14_information_criteria` | AIC and BIC; `part6_information_criteria_*.csv` | 53 |
| `fig15_tris_prediction_z` | held-out prediction of TRIS 600 and 820; `part6_tris_prediction_z*.csv` | 55 |
| `fig16_stacked_residuals` | stacked signed residuals against frequency; `part6_stacked_residuals.csv` | 57 |
| `fig17_error_scale_sensitivity` | the verdict under scaled error bars (maximum likelihood); `part6_error_scale_counts.csv` | 59 |

## Part D — model-choice decision sheet → `paper_outputs/model_choice/`

| File | Shows | Cell |
|---|---|---|
| `decision_sheet.png` | five criteria × three sample sets (all, Cygnus A ± 15°, ± 25° out) | 21 |
| `stacked_residuals_vs_freq.png` | stacked (data − model)/σ per survey; `stacked_residuals.csv` (cell 17) | 23 |
| `dlnZ_vs_RA_with_windows.png` | Δln Z per sample with the exclusion windows | 25 |
| tables | `decision_summary.csv`, `per_sample_values.csv` (cell 27), `loo_standardised_errors.csv` (cell 15) | — |

## Part E — T–T β and the sky models against TRIS → `paper_outputs/phil_plots/`

| File | Shows | Cell |
|---|---|---|
| `plot4_tt_ring_all_pairs` | T–T plots for every pair, York fits; table `plot4_tt_betas_data.csv` (cell 7) | 9 |
| `plot4a_beta_vs_pair_frequency` | β of every pair against √(ν₁ν₂), pairs with Maipu set apart | 11 |
| `plot4c_beta_vs_mean_frequency_span` | β against √(ν₁ν₂), bar from ν₁ to ν₂ | 13 |
| `plot4b_beta_matrix_ra_split` | β matrices for RA < 200° and > 200°; `plot4b_tt_betas_ra_split.csv` | 15 |
| `plot4d_…`, `plot4e_…`, `plot4f_beta_matrix_whole_ring` | Plots 4a, 4c and the whole-ring matrix without the Maipu or RA splits | 17–19 |
| `plot6_tods_600_820_data_vs_models` | TRIS 600/820 MHz TODs with the sky models and their corrections | 24 |
| `plot6a_tris_minus_calibrated_models` | TRIS − (a·model + b) after fitting each model to TRIS; `plot6a_models_fitted_to_tris.csv` (`tab:step1`) | 26 |
| `plot6b_gsm_correction` | a suggested GSM2008 / GSM2016 correction from the two rings; `plot6b_gsm_corrections.csv` | 28 |
| `plot6c_predicted_408_vs_haslam` | the corrected GSMs at 408 MHz against Haslam and its corrections; `plot6c_*.csv` | 30 |
| `plot6d_corrected_gsm_vs_tris_tods` | the corrected GSMs at 600 / 820 MHz with the TRIS rings; `plot6d_*.csv` | 32 |
| `plot5_beta_matrix_data_vs_models` | β of the data pairs and the β each sky model implies; `plot5_model_implied_betas.csv` (cell 34) | 36 |

## Part F — Haslam and the sky models against the prediction → `test_outputs/plot3_haslam_vs_fit/`

| File | Shows | Cell |
|---|---|---|
| `fig1_haslam_vs_power_law_tt` | Haslam (uncorrected, Carter, Irfan) against the prediction, York fit | 8 |
| `fig2_scale_offset_plane` | every 408 MHz model in the (a, b) plane | 10 |
| `fig3_ratio_along_ring` | model / prediction along the ring | 12 |

## Part G — Haslam with one gain and offset → `test_outputs/plot3_haslam_calibration/`

| File | Shows | Cell |
|---|---|---|
| `fig1_haslam_vs_power_law_tt` | Haslam against the prediction with the gain–offset fit | 10 |
| `fig2_calibration_plane` | the fitted (gain, offset) with its ellipse, and the published corrections | 12 |
| `fig3_along_ring` | the Haslam versions and the fitted correction along the ring; `haslam_calibration.csv` (cell 16) | 14 |

## Part H — zero-level offsets → `paper_outputs/zero_level_offsets/`

| File | Shows | Cell |
|---|---|---|
| `figA_offset_profiles` | the profile of each survey's offset; tables `table1_offsets.csv`, `table2_chi2.csv`, `table3*.csv` (cells 12–16) | 18 |
| `figB_spectra_with_without_offsets` | spectral fits at three samples with and without the offsets | 20 |
| `figC_residuals_by_survey_K` | data − model per survey, with and without offsets; `figC_residuals_by_survey.csv` | 22 |

## Part I — robustness of the Haslam gain and offset → `paper_outputs/robustness_check/`

| File | Shows | Cell |
|---|---|---|
| `robustness_gz_plane` | the g–z plane for survey subsets, spectral models, Cygnus A masks, RA halves and zero-level treatment; `robustness_table.csv` (cell 9) | 11 |

## Part J — least-squares fits and the TOD prediction → `test_outputs/tod_prediction/`

| File | Shows | Cell |
|---|---|---|
| `test1_fit_parameters` | the fitted parameters along the ring, with Irfan & Puglisi's β (`irfan_puglisi_2025_ring.csv`, cell 11) | 13 |
| `test2_tris_residuals_K` | the models against TRIS, residuals in K | 15 |
| `test3_chi2_per_frequency`, `test3_residuals_K_all_frequencies` | χ² and residuals at every frequency; `test3_*.csv` | 17 |
| `test4_haslam_408_heldout` | the held-out check at 408 MHz | 19 |
| `test5_predicted_tods` | predicted ring TODs at other frequencies | 21 |

## Part K — the beam

- **K1** (inside the notebook): the σ ↔ FWHM constant, elliptical Gaussians, the archive-cut beam's FWHM measured
  numerically, the E/H cut orientation, a round trip on a beam with a known answer.
- **K2** (inside the notebook, except the beam variants): the three kernels on the sky and along E/H; which kernel
  each survey gets; a point source through the smoothing code; harmonic checks H1–H6 (C_ℓ, limTOD's pointed beam,
  TODs, Haslam's harmonic rotation, harmonic vs real-space smoothing); convolve → deconvolve with limTOD; two-step vs
  one-step beam matching; beam-specification variants → `test_outputs/beam_variants/`: `figA_T_vs_RA_beam_variants`
  (cell 70), `figB_difference_from_V0` (72), `table3_difference_from_V0.csv` (74), `table4_model_vs_TRIS_score.csv`
  (76).
- **K3–K5**: the saved plots of the script tests below.

## `diagnostics.ipynb` (plots inside the notebook)

- **Part A — the inputs:** the TRIS products, temperature convention, beam and prior template; ARCADE 2 and its mask;
  Haslam; the TRIS maps and how unobserved sky is marked; every map in `res/` side by side and along the TRIS stripe;
  the TRIS TOD each map predicts; the nside of every file.
- **Part B — coordinate frames:** what each input declares, pixel ordering, the rotation per map, known sources against
  literature positions, ARCADE 2's mask through the rotation.
- **Part C — the TRIS maps (step 0):** convergence, residuals along the ring, the zero level and its degeneracy, prior
  → posterior shrinkage; prior-only and high-noise runs.

## Script tests (saved files)

**FWHM scaling** — `test_outputs/fwhm_scaling_test/`, `scripts/fwhm_scaling_test.py`

| File | Shows |
|---|---|
| `fwhm_scaling_comparison.png` | two-step / one-step − 1 [%] for each FWHM scaling, one line per kernel |
| `fwhm_scaling_residual_maps.png` | the residual two-step − one-step [K] on the sky, 3 kernels × 4 scalings |

**Beam averaging with spectral-index blobs** — `test_outputs/beam_averaging_blobs/<map>/`, `scripts/beam_averaging_blobs/`

| File | Shows |
|---|---|
| `plotA_shape_change_vs_lst.png` | change in the shape of T vs LST [%] at 600 and 820 MHz, blobs marked |
| `plotB_beam_dilution.png` | largest change in shape without and with the TRIS beam |
| `t_vs_lst_comparison.png` | normalised T, T_ν/T₄₀₈ and shape change vs LST, blobs against a uniform-β control |
| `beta_map.png` | the β map in equatorial with the TRIS ring |

**LWA1 35 MHz and Stockert 1420 MHz** — `test_outputs/lwa35_stockert/`, `scripts/beam_tests_lwa35_stockert1420/`
(write-up in its `README.md`): the FWHM-scaling and blob tests on those two maps, compared with Haslam.

## `paper_outputs/nb08_paper_plots/`

Two of its tables are read by `10_paper_plots` cell 22 for Figure 9: `fig15_part3_g_z_vs_beta.csv` (the β grid and
reference points) and `fig15_part1_haslam_gain_offset.csv` (a check of the TRIS-only fits). Nothing in the repo writes
them (`QUICKSTART.md` open item 6); the other files in the folder are not read.
