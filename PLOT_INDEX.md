# moonrunes — plot index

As of 2026-10-28 (notebooks 11–18 and the five project steps added; ring TODs beam-matched in the pipeline; notebook 11 with the order-3 moment expansion); notebook 08 Figs 15–16 added 2026-11-04; notebooks 19, 20 and 21 and the paper figure `skymodel_tris_anchor_gz` (notebook 10 Figure 8) added 2026-11-05; notebook 19 Test 6 (RHINO) added 2026-11-05, and its paper versions `rhino_tods_tris_beam`, `rhino_tods_rhino_beam` (notebook 10 Figures 9–10) 2026-11-06; paper figures renamed to `<section>_<content>` and all ten placed in `main.tex` 2026-11-06. The online copy <https://claude.ai/code/artifact/a5a7762b-b179-4813-8aa6-d3ea8d198d39> is the 2026-10-15 version.

## How to use this index

Every plot in the project lives in one of three places. Saved files are listed by name; plots that exist only inside a notebook are listed by notebook, cell number and section heading.

| Where | What is there | Made by |
| --- | --- | --- |
| `moonrunes_paper/figures/` (in `Papers_and_reports`) | the paper figures, PDF + PNG (the paper includes only the PDFs), named `<section>_<content>` after the paper section they sit in: all maps over the TRIS slice §2.2; the TRIS beam §3.2; TRIS maps §3.3; sky models vs TRIS §3.4.2; β matrix and β vs frequency §4.1; model-comparison fits vs TRIS §5.1; the 408 MHz models vs the power-law prediction §6.2; the TRIS anchor of the Haslam gain and offset §6.3; the RHINO prediction through the RHINO beam §7.2 (the TRIS-beam version is made but not in the paper). Rebuilt 2026-10-28, anchor figure added 2026-11-05, RHINO figures, the renaming and the 408 MHz figure 2026-11-06 | `notebooks/10_paper_plots.ipynb` |
| `moonrunes/outputs/paper_plots/` | the working figure set (Figures 1–16), PDF + PNG | `notebooks/08_update_plots.ipynb` |
| inside each notebook | diagnostic and exploration plots, not saved to files | notebooks 01–09 |
| `moonrunes/outputs/fwhm_scaling_test/` | the FWHM scaling test (added 2026-10-19) | `scripts/fwhm_scaling_test.py` |
| `moonrunes/beam_averaging_blobs/<map>/` | the spectral-index blob test (added 2026-10-19/20) | `beam_averaging_blobs/*.py` |
| `moonrunes/beam_tests_lwa35_stockert1420/` | the FWHM-scaling (beam matching) and blob (beam averaging) tests on LWA1 35 MHz and Stockert 1420 MHz, compared with Haslam; write-up in its `README.md` (added 2026-10-28) | `beam_tests_lwa35_stockert1420/*.py` |
| `moonrunes/model_comparison/diagnostic_plots/` | Bayesian model comparison diagnostics and convergence tests | `notebooks/11_model_comparison_diagnostics.ipynb` |
| `moonrunes/model_comparison/sky_models_408/` | sky models at 408 MHz against the power-law and curved fits; the fit and TOD cache | `notebooks/12_sky_models_vs_fit_408.ipynb` |
| `moonrunes/model_comparison/phil_plots/` | plots 1–3 (sky models vs the power-law prediction), plots 4–6 (T–T β, models at the TRIS frequencies) | `notebooks/13`, `notebooks/14` |
| `moonrunes/model_comparison/model_choice/` | the model-choice decision sheet | `notebooks/15_model_choice_decision_sheet.ipynb` |
| `moonrunes/model_comparison/plot3_haslam_vs_fit/` | Haslam and the sky models in the scale–offset plane against the power-law prediction | `notebooks/16_haslam_vs_power_law_408.ipynb` |
| `moonrunes/model_comparison/plot3_haslam_calibration/` | Haslam gain–offset calibration | `notebooks/17_haslam_gain_offset_calibration.ipynb` |
| `moonrunes/model_comparison/archive/<run>/` | frozen copies of a model comparison run and its notebook 11 figures, with a write-up (`README.md`); first: `2026-10-27_beam-matched_4-point_run/` | copied from `results/` and `diagnostic_plots/` |
| `moonrunes/model_comparison/beam_matched_tods/` | beam-matched ring TODs of every survey and sky model, and their change from the old TODs | `notebooks/18_beam_matched_tods.ipynb` |
| `moonrunes/tod_prediction/` | power-law and curved fits per ring sample, and the predicted ring TOD at any frequency (Tests 1–5); the prediction at RHINO frequencies through the RHINO beam (Test 6) | `notebooks/19_power_law_tod_prediction.ipynb` |
| `moonrunes/model_comparison/zero_level_offsets/` | zero-level offsets as shared fit parameters: profiles, example SEDs, residuals per survey; Tables 1–3 | `notebooks/20_zero_level_offsets.ipynb` |
| `moonrunes/model_comparison/robustness_check/` | robustness of the Haslam gain and offset across survey subsets, spectral models, sample subsets and zero-level treatment (added 2026-11-05) | `notebooks/21_robustness_check.ipynb` |

Cell numbers are 0-based, as `nbformat` counts them, and were taken on 2026-10-15. They shift if cells are added above. The section heading is the stable way to find a plot.

All SED results here exclude ARCADE 2 (`sed.exclude_maps`); since 2026-10-27 that list also excludes LWA1 35, Maipu 45 and Juno, so the SED fits have 4 points per ring sample (TRIS 600/820, LWA1 80, Stockert 1420). They use the 23.366° circular beam-matching kernel (`beam_matching.kernel: circular_gaussian`) on the map route.

**Which TODs a plot uses.** Since 2026-10-27 the pipeline itself makes *beam-matched* ring TODs (`sed.tod_beam: beam_matched` in `configs/run_config.yaml`, `sed_prototype.Inputs.ring_tod`): the measured TRIS cut beam with each map's own beam taken out, on a 0.92° grid, the same TODs as `notebooks/18`. Everything built from `sed_prototype.Inputs` uses them when rerun: `run_comparison.py` (the `model_comparison/results/` CSVs of 2026-10-27 17:13 already do, `tod_beam` column), notebooks 07, 08, 10, 11, 12, and notebook 14. Caches in notebooks 11, 12 and 15 carry a `tod_key` and are redone when the TODs change. **Not yet rerun** (still the earlier limTOD TODs and the 6-point SEDs): notebook 12's fit cache and so the saved plots of notebooks 12, 13, 15, 16 and 17; notebook 11's saved figures (fig1–4 show three models, fig5–9 the old TODs) until it is run again. The TOD difference itself is small (T–T β change by ≤ 0.01; Haslam −0.02 K median, LWA1 35 −30 K of ~14 000 K); dropping LWA1 35 and Maipu from the fits matters more.

## Plots by project step

The five steps of the paper (`PAPER_OUTLINE.md`).

| Step | Question | Plots |
| --- | --- | --- |
| 1 | Which sky models, or models + corrections, stretched to 600/820 MHz, agree with the TRIS TODs? | `phil_plots/plot6_tods_600_820_data_vs_models.png` (14, Plot 6); `plot6a_tris_minus_calibrated_models.png`, `plot6b_gsm_correction.png`, `plot6c_predicted_408_vs_haslam.png` (14, Plots 6a–6c); `phil_plots/plot5_beta_matrix_data_vs_models.png` (14, Plot 5) |
| 2 | β from T–T plots, and which maps to put in the SED fit | `phil_plots/plot4_tt_ring_all_pairs.png`, `plot4a_beta_vs_pair_frequency.png`, `plot4c_beta_vs_mean_frequency_span.png`, `plot4b_beta_matrix_ra_split.png`, and without the Maipu / RA splits `plot4d_*`, `plot4e_*`, `plot4f_*` (14); 08 Figs 12, 13, 13b |
| 3 | Model comparison with the selected maps, and the winner | `model_choice/decision_sheet.png`, `dlnZ_vs_RA_with_windows.png`, `stacked_residuals_vs_freq.png` (15); `diagnostic_plots/fig1`–`fig9` (11) |
| 4 | The SED-fit sky model against the existing models at 408 MHz | `phil_plots/plot1_legacy_models_vs_fit.png`, `plot2_corrected_models_vs_fit.png`, `plot3_chi2_vs_power_law.png` (13); `plot3_haslam_vs_fit/fig1`–`fig3` (16); `plot3_haslam_calibration/fig1`–`fig3` (17); `sky_models_408/fig1`–`fig3` (12); the 408 MHz models against the prediction through the TRIS beam: paper `skymodel_models_vs_prediction_408` (10); how strongly TRIS anchors the Haslam gain and offset: 08 Figs 15, 16, paper `skymodel_tris_anchor_gz` (10); its robustness across subsets and spectral models: 21; with zero-level offsets fitted: 20 Figs A–C |
| 5 | Foreground predictions for RHINO | `tod_prediction/test6_rhino_freqs_tris_beam`, `test6_rhino_freqs_rhino_beam` (19 Test 6); paper `rhino_tods_rhino_beam` (10; `rhino_tods_tris_beam` made but not in the paper) |
| — | Beam matching of the TODs | `beam_matched_tods/fig1_beam_matched_minus_old.png` (18); the TRIS beam and its orientation: 08 Figs 14, 14b, paper `method_tris_beam` (10); how much the T-vs-LST curve depends on the beam's roll and E–H interpolation: 03 cells 107, 109 (`outputs/beam_variants/`) |

## Notebooks 11–21 (saved files)

Cell numbers are 0-based, taken on 2026-10-27 (notebooks 19–21 on 2026-11-05). Every figure is a PNG.

**Notebook 11 — model comparison diagnostics** → `model_comparison/diagnostic_plots/`

| File | Shows | Cell |
| --- | --- | --- |
| `fig1_dlogZ_vs_ra` | ΔlnZ (curved − PL, moment order 2 − PL, moment order 3 − PL) per ring sample vs RA, tie band, Cygnus A | 6 |
| `fig2_logZ_vs_ra` | ln Z of each of the four models vs RA | 8 |
| `fig3_dlogZ_histograms` | histograms of ΔlnZ, three panels (curved, moment 2nd, moment 3rd − PL) | 10 |
| `fig4_preference_counts` | samples per evidence category, five comparisons (incl. moment 3rd − PL, moment 3rd − 2nd) | 12 |
| `fig5_sed_fits_per_unit` | SED fits of the four models with 68% bands at six ring samples; Haslam plotted, not fitted | 17 |
| `fig6_maxlogl_vs_logZ` | fit gain (Δ max ln L) vs evidence gain; Occam penalty | 19 |
| `fig7_C_posterior` | posterior of the curvature C with the prior | 21 |
| `fig8_prior_sensitivity_C` | curved − PL ΔlnZ vs C prior width | 25 |
| `fig9_convergence_dlogZ` | ΔlnZ (curved, moment 2nd, moment 3rd − PL) from repeated / re-tuned UltraNest runs vs the exact grid integral; tables `convergence_check1`–`4*.csv` | 37 |
| `fig10_models_at_tris_frequencies` | Part 5: each model's maximum-likelihood fit at every ring sample, at 600.5 / 817.8 MHz, with the TRIS rings (published σ); model − TRIS [K]; tables `part5_model_tods_at_tris.csv`, `part5_model_minus_tris_summary.csv` | 42 |
| `fig11_chi2_per_sample` | Part 5: χ² and χ²/ν of every fit along the ring; tables `part5_chi2_summary.csv`, `part5_chi2_per_sample.csv` (with each point's χ² contribution and the best-fit parameters) | 44 |

**Notebook 12 — sky models vs the fits at 408 MHz** → `model_comparison/sky_models_408/`

| File | Shows | Cell |
| --- | --- | --- |
| `fig1_T408_vs_ra` | PL and curved 408 MHz predictions with ±1σ, seven sky models, TRIS scaled with β = −2.7 (labelled) | 9 |
| `fig2_residual_vs_ra` | sky model − fit [K], against each fit | 11 |
| `fig3_chi2_per_sample` | χ²/N per model and fit; table `chi2_sky_models_vs_fits.csv` | 15 |

**Notebook 13 — sky models vs the power-law prediction** → `model_comparison/phil_plots/`

| File | Shows | Cell |
| --- | --- | --- |
| `plot1_legacy_models_vs_fit` | Haslam, GSM2008, GSM2016, LFSM vs the PL prediction; residual [K] | 6 |
| `plot2_corrected_models_vs_fit` | Carter, Irfan, GINAN (and Haslam) vs the PL prediction; residual [K] | 8 |
| `plot3_chi2_vs_power_law` | χ² contribution per ring sample, and χ²/N per model; table `residual_stats.csv` | 12 |

**Notebook 14 — T–T β along the ring, and models at the TRIS frequencies** → `model_comparison/phil_plots/` (beam-matched TODs)

| File | Shows | Cell |
| --- | --- | --- |
| `plot4_tt_ring_all_pairs` | T–T plots for every survey pair, York fits; table `plot4_tt_betas_data.csv` | 8 |
| `plot4a_beta_vs_pair_frequency` | β of every pair vs √(ν₁ν₂), labelled pairs | 10 |
| `plot4c_beta_vs_mean_frequency_span` | β vs √(ν₁ν₂), horizontal bar from ν₁ to ν₂ | 12 |
| `plot4b_beta_matrix_ra_split` | β matrices for RA < 200° and > 200°; table `plot4b_tt_betas_ra_split.csv` | 14 |
| `plot4d_beta_vs_pair_frequency_all_pairs` | Plot 4a with every pair alike (no Maipu split) | 16 |
| `plot4e_beta_vs_mean_frequency_span_all_pairs` | Plot 4c with every pair alike (no Maipu split) | 17 |
| `plot4f_beta_matrix_whole_ring` | β matrix over the whole ring (no RA = 200° split) | 18 |
| `plot6_tods_600_820_data_vs_models` | TRIS 600/820 MHz TODs with GSM2008, GSM2016, the stretched Carter / Irfan / GINAN corrections and stretched LFSM; model − TRIS [K] | 23 |
| `plot6a_tris_minus_calibrated_models` | TRIS − (a·model + b) after fitting each Plot 6 model to TRIS at each frequency; table `plot6a_models_fitted_to_tris.csv` | 25 |
| `plot6b_gsm_correction` | suggested GSM2008 / GSM2016 correction: one constant offset with a constant gain or a gain × (ν/408)^δ, carried to 408 MHz; table `plot6b_gsm_corrections.csv` | 27 |
| `plot6c_predicted_408_vs_haslam` | 408 MHz TOD from the TRIS-corrected GSMs against Haslam, Carter and Irfan; table `plot6c_haslam_vs_predicted_408.csv` | 29 |
| `plot5_beta_matrix_data_vs_models` | β matrix: data pairs and model-implied β; table `plot5_model_implied_betas.csv` | 33 |

**Notebook 15 — model-choice decision sheet** → `model_comparison/model_choice/`

| File | Shows | Cell |
| --- | --- | --- |
| `decision_sheet.png` | five criteria × three sample sets (all, Cygnus ±15°, ±25° excluded) | 20 |
| `stacked_residuals_vs_freq.png` | stacked (data − model)/σ per survey, three models, error scales | 22 |
| `dlnZ_vs_RA_with_windows.png` | ΔlnZ per sample with the exclusion windows; tables `decision_summary.csv`, `per_sample_values.csv`, `loo_standardised_errors.csv`, `stacked_residuals.csv` | 24 |

**Notebook 16 — Haslam vs the power-law prediction** → `model_comparison/plot3_haslam_vs_fit/`

| File | Shows | Cell |
| --- | --- | --- |
| `fig1_haslam_vs_power_law_tt` | Haslam (uncorrected, Carter, Irfan) vs the PL prediction, York fit | 7 |
| `fig2_scale_offset_plane` | all seven 408 MHz models in the (a, b) plane | 9 |
| `fig3_ratio_along_ring` | model / prediction along the ring; table `haslam_vs_power_law_fits.csv` | 11 |

**Notebook 17 — Haslam gain–offset calibration** → `model_comparison/plot3_haslam_calibration/`

| File | Shows | Cell |
| --- | --- | --- |
| `fig1_haslam_vs_power_law_tt` | Haslam vs the PL prediction with the gain–offset fit | 9 |
| `fig2_calibration_plane` | fitted (s, o) with its ellipse, and the published corrections | 11 |
| `fig3_along_ring` | Haslam versions and the fitted correction along the ring; table `haslam_calibration.csv` | 13 |

**Notebook 18 — beam-matched ring TODs** → `model_comparison/beam_matched_tods/`

| File | Shows | Cell |
| --- | --- | --- |
| `fig1_beam_matched_minus_old` | beam-matched − old TOD [K] for the surveys and the 408 MHz models; tables `beam_matched_minus_old.csv`, `beam_matched_tods.csv`, `beam_matched_tods_meta.csv` | 9 |

**Notebook 19 — power-law TOD prediction** → `tod_prediction/`

| File | Shows | Cell |
| --- | --- | --- |
| `test1_fit_parameters` | Test 1: the fitted parameters along the ring; tables `fit_per_sample_<model>.csv` | 12 |
| `test2_tris_residuals_K` | Test 2: the models against TRIS, residuals in K | 14 |
| `test3_chi2_per_frequency`, `test3_residuals_K_all_frequencies` | Test 3: χ² of the fit at every frequency, and residuals in K; tables `test3_chi2_per_frequency.csv`, `test3_chi2_per_sample.csv` | 16 |
| `test4_haslam_408_heldout` | Test 4: the held-out check at 408 MHz (Haslam) | 18 |
| `test5_predicted_tods` | Test 5: predicted ring TODs at other frequencies; tables `predicted_tods_<model>.csv` (and `_with_cmb_`) | 20 |
| `test6_rhino_freqs_tris_beam` | Test 6: the power-law prediction at 60–80 MHz through the TRIS beam, against GSM2008, GSM2016, LFSM, GSM2008 + Carter, GSM2008 + Irfan and GSM2016 + GINAN (one panel per frequency); table `test6_rhino_freqs_tris_beam.csv` (added 2026-11-05) | 23 |
| `test6_rhino_freqs_rhino_beam` | Test 6: the same, deconvolved from the TRIS beam (Wiener, GSM2016 prior) and observed with the RHINO DryHornGround beam at the Jodrell Bank zenith; table `test6_rhino_freqs_rhino_beam.csv` (added 2026-11-05) | 23 |

**Notebook 20 — zero-level offsets** → `model_comparison/zero_level_offsets/`

| File | Shows | Cell |
| --- | --- | --- |
| `figA_offset_profiles` | ΔF profile of each survey's offset over ±10 prior widths, for configurations (1)–(3), with the TRIS II Table 3 zero-level uncertainties marked; tables `table1_offsets.csv`, `table2_chi2.csv`, `table3a/b/c_*.csv` (cells 11, 13, 15) | 17 |
| `figB_spectra_with_without_offsets` | spectral fits at three samples with and without the configuration (1) offsets; data − model [K] | 19 |
| `figC_residuals_by_survey_K` | data − model [K] at each survey over all unmasked samples, baseline and with offsets; table `figC_residuals_by_survey.csv` | 21 |

**Notebook 21 — robustness check of the Haslam gain and offset** → `model_comparison/robustness_check/`

| File | Shows | Cell |
| --- | --- | --- |
| `robustness_gz_plane` | the g–z plane (Haslam = g P + z), two panels: survey subsets (power law and curved power law), and Cygnus A mask width, RA halves and notebook 20's zero-level rows; the baseline (six-survey power law, 103 samples) error box; table `robustness_table.csv` (cell 8) | 10 |

## Paper figures — `moonrunes_paper/figures/`

Eleven figures, ten of them in `main.tex` (`rhino_tods_tris_beam` was taken out of the paper on 2026-11-06), made by `notebooks/10_paper_plots.ipynb` (rebuilt 2026-10-28), drawn at their printed MNRAS size, PDF + PNG, under the file names `main.tex` includes. `main.tex` uses only the PDFs. File names are `<section>_<content>`, the prefix being the paper section (`data`, `method`, `tt`, `modelcomp`, `skymodel`, `rhino`, as in the `sec:` labels); the old names are in brackets. Captions are in `main.tex`. Each notebook section starts with a note on where the figure is featured.

| File | Shows | Paper section (label) | Notebook 10 cell |
| --- | --- | --- | --- |
| `data_maps_tris_slice` (was `fig1_all_maps_tris_slice`) | the five surveys that reach the ring, native resolution, over the TRIS declination band | 2.2 Surveys (`fig:maps-slice`) | 5 |
| `method_tris_maps` (was `tris_maps`) | TRIS 600 and 820 MHz maps from stage 1 (no 2.5 GHz panel any more) | 3.3 Maps from the TRIS rings (`fig:tris-maps`) | 7 |
| `method_models_vs_tris_tods` (was `plot6_tods_600_820_data_vs_models`) | sky models and stretched corrections at 600.5 / 817.8 MHz vs TRIS; model − TRIS [K] | 3.4.2 The TOD comparison (`fig:tods-600-820`) | 9 |
| `tt_beta_matrix` (was `plot4f_beta_matrix_whole_ring`) | T–T β ± σ and χ²/ν for every survey pair, whole ring | 4.1 T–T fits (`fig:beta-matrix-ring`) | 13 |
| `tt_beta_vs_frequency` (was `plot4c_beta_vs_mean_frequency_span`) | T–T β vs √(ν₁ν₂), bar from ν₁ to ν₂; one column wide | 4.1 T–T fits (`fig:beta-vs-freq`) | 15 |
| `modelcomp_fits_vs_tris` (was `model_comparison_fits_vs_tris`) | maximum-likelihood power law, curved power law and moment expansion at 600.5 / 817.8 MHz vs TRIS; model − TRIS [K] | 5.1 Models, priors and sampler (`fig:modelcomp-tris`) | 17 |
| `method_tris_beam` (was `tris_beam`) | the TRIS beam on the sky (dB, north up, east left) with the E/H axes as dashed crosshairs, and the E- and H-plane cuts; the same figure as 08 Fig 14 at MNRAS size (added 2026-11-02) | 3.2 Beam-matched ring TODs (`fig:tris-beam`) | 19 |
| `skymodel_tris_anchor_gz` (was `tris_anchor_gz_plane`) | how strongly TRIS anchors the Haslam gain and offset: the g–z plane with the 1σ ellipses of TRIS only, without TRIS and all six surveys, and TRIS only against a fixed extrapolation β; one column wide; the same figure as 08 Fig 15 at MNRAS size (added 2026-11-05) | 6.3 The mismatch at 408 MHz (`fig:tris-anchor`) | 22 |
| `skymodel_models_vs_prediction_408` | the six-survey power-law prediction at 408 MHz (±1σ) against the beam-matched 408 MHz TODs of Haslam, Haslam + Carter, Haslam + Irfan, GSM2008, GSM2016, GSM2016 + GINAN and LFSM along the ring, through the TRIS beam; model − prediction [K]; full width; reads notebook 12's fit cache and notebook 18's TODs (added 2026-11-06) | 6.2 Comparison with the existing models at 408 MHz (`fig:models-vs-prediction-408`) | 27 |
| `rhino_tods_tris_beam` | this work's predicted ring TODs (power law ±1σ) at 60, 65, 70, 75, 80 MHz against GSM2008, GSM2016, LFSM, GSM2008 + Carter, GSM2008 + Irfan and GSM2016 + GINAN, all through the TRIS beam; full width; notebook 19 Test 6's figure at MNRAS size (added 2026-11-06) | not in `main.tex` (removed 2026-11-06; its numbers stay in 7.1, The prediction along the TRIS ring) | 24 |
| `rhino_tods_rhino_beam` | the same, deconvolved from the TRIS beam (GSM2016 prior) and observed with the RHINO `DryHornGround` beam at the Jodrell Bank zenith; full width (added 2026-11-06) | 7.2 From the TRIS ring to the RHINO beam (`fig:rhino-rhino-beam`) | 25 |

## Notebook 08 — update plots (`outputs/paper_plots/`)

`notebooks/08_update_plots.ipynb`, target = destriped-only Haslam. Every figure is saved as `<file>.pdf` and `.png`.

| Fig | File | Shows | Cell |
| --- | --- | --- | --- |
| 1 | `fig01_maps_used` | every map in the SED fit, native resolution, with the TRIS stripe | 4 |
| 2 | `fig02_target_stripe` | the destriped Haslam target in equatorial, with the stripe | 6 |
| 3 | `fig03_stripe_slices` | RA/Dec slices of every map along the TRIS stripe | 8 |
| 4 | `fig04_tris_maps_and_tods` | the stage 1 TRIS maps and the ring TODs | 10 |
| 5 | `fig05_sed_fits` | 12 pixel-by-pixel SED fits (map route) with residuals | 12 |
| 5b | `fig05b_sed_fits_tod` | 12 sample-by-sample SED fits (TOD route) with residuals | 14 |
| 6 | `fig06_fit_residual_maps` | best fit − map in K, every survey, over the stripe | 16 |
| 6b | `fig06b_tod_residuals` | best fit − data along the ring, every survey | 18 |
| 7 | `fig07_bestfit_haslam_both_routes` | best-fit Haslam from both routes, and each minus old Haslam | 20 |
| 7b | `fig07b_tod_route_ring` | old Haslam TOD vs the best-fit 408 MHz prediction along the ring | 22 |
| 7c | `fig07c_fractional_difference_hist` | histogram of (best fit − old)/old, both routes | 24 |
| 7d | `fig07d_scaled_tods_vs_lst` | every map's TOD scaled to 408 MHz by a power law fitted without TRIS, vs LST, with GSM2008, GSM2016 and LFSM; dashed lines mark the Haslam survey-zone boundaries along the ring | 26 |
| 8 | `fig08_residuals_by_survey` | pixel-by-pixel residuals (data − fit)/σ, per survey | 29 |
| 9 | `fig09a_chi2_maps`, `fig09b_chi2_tod` | χ²/dof maps (map route) and along the ring (TOD route) | 31 |
| 10 | `fig10a_beta_maps`, `fig10b_beta_tod` | β maps (map route) and along the ring (TOD route) | 33 |
| 11 | `fig11_corrected_full_resolution` | corrected Haslam at full resolution, both routes | 35 |
| 12 | `fig12a_tt_raw`, `fig12b_tt_beam_convolved` | T–T corner plots, raw and beam-convolved | 37 |
| 13 | `fig13_tt_beta` | β from every T–T pair, raw and beam-convolved, with bootstrap errors | 39 |
| 13b | `fig13b_beta_matrix_tt_vs_sed` | β matrix: T–T (below the diagonal) vs map-route SED fit (above), and their difference | 41 |
| 14 | `fig14_tris_beam_tod_route` | the TRIS beam as the beam-matched TOD route uses it (`tris_cuts` kernel, pipeline code), two panels without annotations: the beam on the sky at the zenith in dB (north up, east left) with the E-plane (7°) and H-plane (97°) axes as white dashed crosshairs; the E- and H-plane cuts. The cell prints the FWHM check (19.155° / 23.366°). Added 2026-11-02, trimmed the same day | 46 |
| 14b | `fig14b_tris_beam_limtod_pointing` | the same beam pointed by limTOD (`pointing_beam_in_eq_sys`, azimuth 0°, elevation 90°, roll −7°) on the RA/Dec grid; fitted orientation 6.99° / 96.99°, against the kernel's 7° / 97° (added 2026-11-02) | 48 |
| 15 | `fig15_tris_anchor_gz_plane` | how strongly TRIS anchors the Haslam gain and offset: the g–z plane (Haslam = g P + z, York, 103 samples, Cygnus A masked, errors as now) with the 1σ ellipses for three predictions P of Haslam at 408 MHz: (a) TRIS only, (b) a power law through the four other surveys, (c) all six surveys; and the TRIS-only fit with the extrapolation index fixed from −2.2 to −3.0, with three reference indices marked. CSVs: `fig15_part1_*`, `fig15_part3_*`. Added 2026-11-04 | 58 |
| 16 | `fig16_sigma_lnA408_leave_out` | Fisher σ(ln A at 408 MHz) of the six-survey power law (median over samples), with all six surveys, each survey left out, and both TRIS bands left out. CSVs: `fig16_part2_*` (σ values and each survey's share of the information). Added 2026-11-04 | 60 |

Also in the folder:

- `haslam_zones_tris_ring_check.png`: the TRIS ring traced over Fig. 1 of Nasirudin & Bull (2025), coloured by the Haslam zone it is read as. It is the check behind the zone lines in Figure 7d. Not made by a notebook cell.
- `fig07a_old_vs_bestfit_map_route` and `fig07b_old_vs_bestfit_tod_route` are leftovers from an earlier version of Figure 7. No current cell writes them.

## Notebook 07 — SED calibration

`notebooks/07_sed_calibration.ipynb`, target = the standard Haslam map. Plots stay in the notebook.

| Cell | Section | Shows |
| --- | --- | --- |
| 10 | 3. Calibration stage 1 — the TOD route | Haslam TOD vs predicted 408 MHz; prediction − Haslam; β and χ²/dof vs RA |
| 12 | β vs LST with different points left out | β along the ring with each survey left out, and with groups left out (no TRIS, no surveys below 100 MHz) |
| 18 | 5. The maps | 3 figures: old and TOD-route corrected Haslam; map-route corrected, correction and hybrid maps; map-route β, χ²/dof and points per fit |
| 20 | 6. The two routes against each other | map-route vs TOD-route correction, per stripe pixel |

## Notebook 03 — diagnostics

`notebooks/03_diagnostics.ipynb` has six parts. Plots stay in the notebook, except the beam-variant section (cells 100–114, added 2026-11-02), which also saves its figures and tables to `outputs/beam_variants/`.

| Cell | Part / section | Shows |
| --- | --- | --- |
| 5 | Stage 1 · 2. Does the model describe the ring? | ring data vs model, residuals |
| 7 | Stage 1 · the uncertainties the solver uses | the per-sample error bars |
| 9 | Stage 1 · 3. The zero level | fitted zero level against the template deficit |
| 12 | Stage 1 · 4. Did the data move the map? | prior → posterior σ shrinkage, z-scores |
| 14 | Stage 1 · the band, drawn | the solved band maps |
| 16, 17 | Stage 1 · 5. Is this fit to hand to stage 3? | β between the solved maps; the degrade to nside 16 |
| 24, 26 | Stage 1 · prior-only runs | what the data moved: histograms and sky maps |
| 30, 32, 34 | Stage 1 · high-noise run | original vs high-noise run; per-pixel σ; per-pixel RMS error |
| 37 | Stage 1 · residuals with error bars | residuals with error bars, and the RMS maps |
| 42 | Elliptical beam · 1 | the three kernels on the sky |
| 44 | Elliptical beam · 2 | kernel profiles along E and H, and weight in the shoulders |
| 47 | Elliptical beam · 3 | which kernel each survey gets |
| 49 | Elliptical beam · 4 | a point source through the smoothing code; broad native beams |
| 52 | Elliptical beam · 5 | Haslam with the old vs new kernel |
| 55 | Elliptical beam · 6 | ARCADE 2's mask edge under both kernels |
| 58 | Elliptical beam · 7 | how the new kernel changes every map in the stripe |
| 64 | Harmonic checks · H1 | C_ℓ of every map, raw and after smoothing |
| 67 | Harmonic checks · H2 | the TRIS beam's C_ℓ at nside 16–128 |
| 69 | Harmonic checks · H3 | limTOD's pointed beam vs the exact beam |
| 71 | Harmonic checks · H4 | TODs: harmonic vs exact, at nside 16, 64, 128 |
| 73 | Harmonic checks · H5 | Haslam rotation round trip: ringing at Cas A |
| 75 | Harmonic checks · H6 | harmonic vs real-space Gaussian smoothing |
| 78 | Convolve → deconvolve | the two pointed beams and the operators' singular values |
| 80 | Convolve → deconvolve | Haslam: input, convolved, deconvolved, residual (both beams) |
| 83 | Convolve → deconvolve · the two beams | elliptical − circular, convolved and deconvolved |
| 85 | Convolve → deconvolve · which scales come back | transfer function, correlation by ℓ, error vs smoothing scale |
| 87 | Convolve → deconvolve · every map | rms error per survey, both beams |
| 89 | Convolve → deconvolve · along the TRIS ring | map / input − 1 vs LST, pipeline noise, both beams |
| 93 | Two-step vs one-step | the complicated survey beam: image and profiles |
| 94 | Two-step vs one-step | difference maps: survey beam × kernel |
| 96 | Two-step vs one-step | the difference along the TRIS ring, and a summary |
| 103 | Beam variants · roll check | table: narrow-axis position angle of each variant's limTOD-pointed beam (the roll applied once) |
| 107 | Beam variants · Fig A | `figA_T_vs_RA_beam_variants.png`: T vs RA (= LST) through 7 beam variants (roll +7°/−7°/0°/97°; dB, power, elliptical-FWHM and linear-FWHM interpolations), Haslam 408 and GSM2016 600.5/817.8, V0 in black, mask RA 299.9 ± 25° shaded |
| 109 | Beam variants · Fig B | `figB_difference_from_V0.png`: variant − V0 in K and in %, TRIS ±1σ shaded; tables `table3_difference_from_V0.csv`, `table4_model_vs_TRIS_score.csv` (cells 111, 113) |

## Notebooks 01, 02, 04, 05, 06, 09

Plots stay in the notebooks.

| Notebook | Cell | Section | Shows |
| --- | --- | --- | --- |
| 01 explore data | 5 | the four TRIS products | the 600 / 820 MHz rings with errors |
| 01 | 9 | the temperature convention | β between the rings, with and without the CMB |
| 01 | 11 | the beam | archive E/H cuts vs a Gaussian |
| 01 | 15 | the prior template | GSM2008 through the TRIS beam vs the rings |
| 01 | 21, 24 | ARCADE 2 3.15 / 3.41 GHz | the maps and their masks |
| 01 | 30 | Haslam 408 MHz | the map |
| 01 | 37 | TRIS maps | the stage 1 maps |
| 01 | 47, 49, 51 | the annulus; nside 8 vs 16 | the TRIS footprint and its declination profile |
| 01 | 55, 56 | every map in `res/` | all maps side by side; slices along the TRIS stripe |
| 01 | 61 | the TRIS TOD each map predicts | one TOD panel per map, with coverage |
| 01 | 63, 65 | TODs against each other | normalised by the mean; scaled with β = −2.7 |
| 01 | 67 | spectral index against RA | β of each map vs TRIS 600 along the ring |
| 02 Berkhuijsen prior (superseded) | 4–13 | Haslam; region operators; prior circularity; independent prior | 6 figures |
| 04 coordinate frames | 5 | TRIS ring peak | does the ring peak at the plane crossing? |
| 04 | 14 | sanity check | known sources against literature positions |
| 04 | 17 | ARCADE 2 mask | does the zero mask survive rotation? |
| 05 beam matching (circular kernel) | 9, 10 | before and after | Haslam and ARCADE 2 smoothed to 23.366° |
| 05 | 17 | usable overlap | the four-way overlap with TRIS |
| 06 SED fit (TRIS + ARCADE 2 cross-check) | 10 | sanity checks | the four-point SEDs |
| 06 | 12 | full maps | best-fit Haslam, β and χ² |
| 09 beam FWHM tests | 11 | 4. archive-cut beam FWHM | cut beam vs Gaussian, in dB |
| 09 | 15 | 6. Test 2: round trip | FWHM error vs nside |

## Standalone tests (saved files)

Both are described, with rerun commands, under "Standalone tests" in `README.md`.

**FWHM scaling** — `outputs/fwhm_scaling_test/`, made by `scripts/fwhm_scaling_test.py`

| File | Shows |
| --- | --- |
| `fwhm_scaling_comparison.png` | two-step / one-step − 1 [%]: median, std (and std in the TRIS stripe), max \|Δ\|, for each FWHM scaling, one line per kernel |
| `fwhm_scaling_residual_maps.png` | the residual two-step − one-step [K] on the sky, 3 kernels × 4 scalings, one shared scale |

**Beam averaging with spectral-index blobs** — `beam_averaging_blobs/<map>/` (`HASLAM_DESTRIPED_ONLY/` is the recorded run; `HASLAM/` the earlier one)

| File | Shows | Made by |
| --- | --- | --- |
| `plotA_shape_change_vs_lst.png` | change in shape of T vs LST [%] at 600 and 820 MHz, blobs marked, largest change annotated | `replot_from_csv.py` |
| `plotB_beam_dilution.png` | largest change in shape without vs with the TRIS beam (bars, separate axes) | `replot_from_csv.py` |
| `t_vs_lst_comparison.png` | per frequency: normalised T, T_ν/T₄₀₈ and shape change vs LST, blobs vs uniform-β control | `beam_averaging_blobs.py` |
| `beta_map.png` | the β map in equatorial with the TRIS ring | `beam_averaging_blobs.py` |

## Find a plot by topic

| Question | Where to look |
| --- | --- |
| How low is Haslam against TRIS? | 08 Fig 7b; paper `sed_fit_tod_route`; 07 cell 10 |
| What do the SED fits look like? | paper `sed_fit_tod_route`, `sed_fit_map_route`; 08 Figs 5, 5b |
| Which surveys disagree with the power law? | 08 Figs 6, 6b, 8; 08 Fig 13b |
| How does β vary along the ring / over the sky? | 08 Fig 10; 07 cell 12 (β with surveys left out) |
| Is one power law a good fit? | 08 Fig 9 (χ²); 07 cell 12 |
| How do TRIS, Haslam and the sky models compare? | 08 Fig 7d; at 600/820 MHz: 14 Plot 6; at 408 MHz vs the fit: 13 Plots 1–3, 16, 17 |
| Which sky model or correction agrees best with TRIS? | 14 Plot 6 (step 1) |
| Which spectral model wins (power law, curved, moment)? | 15 decision sheet; 11 Figs 1–8 |
| Is the Bayesian evidence converged? | 11 Fig 9 and `convergence_check*.csv` |
| Haslam's gain and offset against the SED fit | 16 Figs 1–3; 17 Figs 1–3 |
| How strongly does TRIS anchor the Haslam gain and offset? How much does it depend on the extrapolation index? | 08 Figs 15, 16; paper `skymodel_tris_anchor_gz` (10) |
| Is the Haslam gain and offset robust across survey subsets, spectral models, masks and zero-level treatment? | 21 `robustness_gz_plane` |
| Should the zero levels be fitted as offsets, and what changes if they are? | 20 Figs A–C |
| What ring TOD does the fit predict at another frequency? | 19 Test 5 |
| What does the fit predict for RHINO (60–80 MHz), against GSM2008, GSM2016, LFSM and the Carter, Irfan and GINAN corrections? | 19 Test 6 |
| How much does beam matching change the TODs? | 18 Fig 1 |
| Where are the Haslam survey zones along the ring? | 08 Fig 7d (dashed lines); `outputs/paper_plots/haslam_zones_tris_ring_check.png` |
| β from T–T plots, and vs the SED fit | 14 Plots 4, 4a, 4c, 4b, 5 (beam-matched); 08 Figs 12, 13, 13b |
| The corrected Haslam map | 08 Figs 7, 11; 07 cell 18 |
| What the TRIS maps look like | paper `method_tris_maps`; 08 Fig 4; 01 cell 37 |
| What beam matching does to a map | paper `beam_matching_stockert`; 05 cells 9–10 |
| Every survey through the TRIS beam | paper `tods_all_maps`; 01 cells 61–67 |
| Is the beam width / FWHM right? | 09 cells 11, 15 |
| How is the TRIS beam oriented on the sky (E/H planes, the 7° roll)? | 08 Figs 14, 14b |
| Does the ring TOD depend on the beam's roll or on how E and H are interpolated? | 03 cells 107, 109 (`outputs/beam_variants/`) |
| Elliptical vs circular kernel | 03 cells 42–58 |
| Harmonic-space ringing, C_ℓ | 03 cells 64–75 |
| Convolve then deconvolve | 03 cells 78–89 |
| Two-step vs one-step beam matching | 03 cells 93–96 |
| Which FWHM to give the pipeline for an elliptical survey beam | `outputs/fwhm_scaling_test/` |
| Does beam averaging change the shape of T vs LST with varying β? | `beam_averaging_blobs/HASLAM_DESTRIPED_ONLY/plotA_*`, `plotB_*` |
| Stage 1 fit quality and residuals | 03 cells 5–37 |
| Coordinate frames | 04 cells 5–17 |
