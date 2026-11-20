# Handover

What someone taking this project over needs: how to reproduce the paper's results, which code makes which figure,
table and number, the environment, the data that may not be shared, constants defined in more than one place, the
state of the stored notebook outputs, and the open items. The method itself is in `docs/PIPELINE.md`; the datasets in
`docs/DATA_SOURCES.md`; every plot in `docs/PLOT_INDEX.md`; the decisions in force in `CHANGELOG.md`.

The paper lives in a separate repository (`moonrunes_paper`, `main.tex`). Paper items are referred to by their
LaTeX labels, which do not move when the text is edited.

---

## 1. Reproducing the results, in order

1. **The pipeline:** `python run_pipeline.py --step all` (one to two hours). Step 0 writes the TRIS maps to
   `outputs/stage1/`; steps 1–7 write `paper_outputs/pipeline/step*/`; the report is
   `paper_outputs/pipeline/PIPELINE_REPORT.md`. Every step writes a `manifest.json` with the settings it used.
   The `.npz` products are git-ignored, so a fresh clone must run this first.
2. **`notebooks/tests_analysis/`**, one notebook per Part (`README.md` there gives the order: A, then B before D–I and
   K2, H before I). They read the pipeline's products. Parts C–I and L make paper numbers; A, B, J and K are checks.
   Part D runs about 840 UltraNest fits (tens of minutes); K2 takes about 35 minutes.
3. **`notebooks/10_paper_plots.ipynb`** writes the figures straight into the paper repository's `figures/` folder
   (the path is set in cell 1, `OUT`; change it on another machine). Back up the old figures first. The last line of
   cell 25 writes `rhino_tods_tris_beam`, which the paper does not use; skip it unless that file is wanted.
4. `notebooks/diagnostics.ipynb` is independent: checks of the input files, the coordinate frames and the TRIS maps.

`scripts/model_comparison/run_comparison.py` is a standalone version of the model comparison (the same configs, plus
the map route); the pipeline does not use it. Its last output is in `paper_outputs/model_comparison/results/`.

## 2. Environment

`environment.yml` describes the conda environment `limtod` (Python 3.14, healpy 1.20.0, matplotlib 3.11.1, pygdsm
1.7.1, ultranest 4.6.3, joblib 1.6.0). limTOD is imported from the checkout `../limTOD` (`paths.limtod`;
github.com/zzhang0123/limTOD at commit `dbe720b`), which has changes to `limTOD.tris` that the PyPI release does not;
do not `pip install limTOD`. The TRIS archive files are read from `../limTOD/downloads/TRIS`
(`paths.tris_archive_dir`); the RHINO beam and horizon mask from `../rhino-beams/` (step 7 only). The `tests_analysis`
Parts find the repository root themselves (from any folder inside it); `10_paper_plots` and `diagnostics` must be run
from `notebooks/` or the repository root. Headless: `jupyter execute <notebook>`.

## 3. Paper figures

All made by `notebooks/10_paper_plots.ipynb` (cell numbers 0-based), as PDF + PNG. Ring = 120 samples, RA 0–357° in 3°
steps, Dec +42.433°; beam = `tris_cuts`, 19.155° × 23.366° (E plane 7° east of the meridian); samples kept where the
observed beam weight is ≥ 0.5. "Step 5's data sets" = every map `sed.exclude_maps` leaves in with data on the ring
(LWA1 35, 38, 40, 45, 50, 60, 70, 74, 80; Maipu 45; Stockert 1420) plus TRIS 600.5 and 817.8: 13 points.

| File | Label | Cell(s) | Inputs and configuration |
|---|---|---|---|
| `fig1_data_maps_tris_slice` | `fig:maps-slice` | 5 | the eleven survey maps of step 5's data sets and Haslam, native resolution, no beam, two columns |
| `fig2_method_tris_beam` | `fig:tris-beam` | 20 | `beam_spec_from_config(kernel="tris_cuts")`; asserts 19.155° / 23.366° |
| `fig3_method_tris_maps` | `fig:tris-maps` | 7 | step 0's TRIS maps; the fitted zero level is not added |
| `fig4_method_models_vs_tris_tods` | `fig:tods-600-820` | 10 | TRIS ring data (labelled on the plot); model TODs from step 2 (`beam_matched_tods.npz`); Carter and GINAN carried with the base model's spectrum, Irfan with Irfan & Puglisi's β_s and c_s |
| `fig5_tt_beta_matrix` | `fig:beta-matrix-ring` | 12 (York fits), 14 | 14 ring TODs (step 5's data sets + Haslam), 90 pairs; York, σ × √(χ²/ν); pairs need ≥ 20 samples; LWA1 45 vs Maipu 45 skipped (equal frequencies). Drawn full width |
| `fig6_tt_beta_vs_frequency` | `fig:beta-vs-freq` | 16 | as Figure 5 |
| `fig7_modelcomp_fits_vs_tris` | `fig:modelcomp-tris` | 18 | its own maximum-likelihood fits (24 starts, in parallel) on step 5's data sets; priors from `configs/model_comp_config.yaml` |
| `fig8_skymodel_models_vs_prediction_408` | `fig:models-vs-prediction-408` | 28 | our model at 408 MHz (step 5: maximum-likelihood A, posterior σ) against step 2's 408 MHz model TODs |
| `fig9_skymodel_tris_anchor_gz` | `fig:tris-anchor` | 22–23 | Haslam = g·P + z (York), 103 samples (Cygnus A ± 25° out), step 5's data and errors; the β grid and reference points are read from `paper_outputs/nb08_paper_plots/fig15_part3_g_z_vs_beta.csv` and the TRIS-only fits checked against `fig15_part1_haslam_gain_offset.csv` (open item 6) |
| `fig10_rhino_tods_rhino_beam` | `fig:rhino-rhino-beam` | 25 (function), 26 | step 7's `rhino_tods_rhino_beam.csv`: our model, Wiener-deconvolved with a GSM2016 prior, observed with the RHINO beam at the Jodrell Bank zenith |

## 4. Paper tables and quoted numbers

The tables are typed into `main.tex` by hand.

| Paper item | Made by | Configuration to state |
|---|---|---|
| `tab:surveys` | by hand, from `docs/DATA_SOURCES.md` and `configs/error_bars.csv` | survey uncertainties as published (`CHANGELOG.md`) |
| `tab:step1` | Part E, cell 26 → `paper_outputs/phil_plots/plot6a_models_fitted_to_tris.csv` | **TRIS = a·model + b**, weighted least squares, errors on TRIS only, covariance × χ²/ν |
| `tab:408-residuals` | step 6 → `paper_outputs/pipeline/step6_predict_408/residuals_408.csv` (printed by `10_paper_plots` cell 28) | χ² with σ = the prediction's posterior σ only |
| `tab:scale-offset` | step 6 → `scale_offset_408.csv` (plotted by Part F) | **model = a·P + b**, York, errors on both axes, σ inflated; every model given Haslam's σ (10 % ⊕ 3 K); Cygnus A ± 25° out → 103 samples |
| Section 5 counts (wins, ties, borderline samples) | Part C, cells 13 and 49 (`part6_counts_masked.csv`), from step 3 | step 3's data sets: 11 surveys, **TRIS excluded by design** |
| Stacked residuals by survey | Part D, cell 17 → `paper_outputs/model_choice/stacked_residuals.csv` | criterion 1 uses step 3 (no TRIS); criteria 2–5 use fits with TRIS |
| 408 MHz prediction (median and σ) | step 6 → `prediction_408.csv`; also `10_paper_plots` cell 28 | our model (step 5) |
| Zero levels as fitted offsets | Part H, cells 12 and 14 → `paper_outputs/zero_level_offsets/table1_offsets.csv`, `table2_chi2.csv` | one offset per survey, shared by all samples, priors from the published systematics |
| g, z robustness | Part I, cell 9 → `paper_outputs/robustness_check/robustness_table.csv` | baseline: step 5's data sets, 103 samples |
| Haslam against the prediction (ratios) | Part F (reads step 6) | York, as `tab:scale-offset` |
| Section 7 (RHINO) numbers | step 7 → `rhino_summary.csv`, `rhino_checks.csv`; `10_paper_plots` cells 25–26 print the medians | our model; the closure, the uncertainty limits and the Irfan split are in `rhino_checks.csv` (`CHANGELOG.md`) |
| Irfan spread along the ring at 600 / 820 MHz | step 2 → `irfan_split_tris.csv` (`range (max − min) total (K)`) | Irfan-corrected model − base model, through the TRIS beam |
| Fisher shares, σ(ln A) leave-out, g at the median index | Part L, cells 9, 11, 13 → `paper_outputs/fisher_share/` | Fisher (ln A, β) of T = A (ν/408)^β at step 5's fit, 103 samples |
| Beam tests | Parts K3 (`scripts/beam_averaging_blobs/`), K4 (`scripts/fwhm_scaling_test.py`), K5 (`scripts/beam_tests_lwa35_stockert1420/`), shown from saved results | — |
| TRIS zero levels +1.74 / +0.52 K; 23 / 7 of 1272 pixels narrowed > 5 %; posterior/prior width 0.998 / 0.999 | step 0 → `outputs/stage1/stage1_manifest.json` (`frequencies.<f>.zero_level_k`, `pixels_shrunk_gt_5pct`, `sigma_shrinkage_median`); `10_paper_plots` cell 7 prints them | nside 16, GSM2008 prior; the zero level is fitted, not added to the maps |
| TRIS stripe Dec +7.8° to +77.2° | `10_paper_plots` cell 1 (`DEC_LO`, `DEC_HI`) | the 1272 pixels where the beam response exceeds 0.01 of its peak |

## 5. Data that may not be shared

Every dataset, its source and its uncertainty: `docs/DATA_SOURCES.md`. Do **not** copy or redistribute:

- `res/Juno_600MHz.fits`, `res/Juno_1248MHz.fits` (git-ignored): no provenance (the header has only NSIDE and
  ORDERING; no unit, frame, beam or citation). Excluded from every fit.
- `external/irfan_puglisi_2025/*.npy` (git-ignored): Irfan & Puglisi 2025's β_s and c_s maps; the licence is not
  cleared with the authors. Used by the Irfan correction off 408 MHz (steps 2 and 7) and by Part J's β comparison.
- Anything derived from `../rhino-beams/` (the RHINO beam and the Jodrell Bank horizon mask): used by step 7 and so
  by Figure 10.
- Logs or outputs that contain local absolute paths.

`res/CHIPASS_1400MHz.fits` is git-ignored (too large). The other maps in `res/` are LAMBDA products.

## 6. Constants defined in more than one place

Not merged; listed so that a change is made everywhere. "PP" = `notebooks/10_paper_plots.ipynb`; "A c2:15" = Part A,
cell 2, line 15.

| Constant | Where | Agree? |
|---|---|---|
| Cygnus A RA 299.9° | `run_pipeline.yaml` `pipeline.predict_408.cygnus_ra_deg`; A c2:15, B c2:52, C c3:56, D c2:24, F c2:18, G c2:20, H c2:34, I c3:48, L c1:20; PP c1:43 | yes |
| Cygnus half-window | 25° in `pipeline.predict_408.cygnus_half_width_deg`, C c47:2, F c2:18, G c2:20, I c3:48, L c1:20, PP c22:13; **15°** as an alternative in C c47:3 and D c2:25 | two windows, by design |
| Carter 1.029, 0.91 K | `run_pipeline.yaml` `pipeline.sky_models.carter`; B c2:36 | yes |
| Irfan monopole 8.9 K | `pipeline.sky_models.irfan_monopole_k`; B c2:33 | yes |
| Irfan ν₀ = 45 MHz and the β_s / c_s files | `pipeline.sky_models.irfan_*`; Part J cell 11 | yes |
| GINAN anchor 350 MHz | `pipeline.sky_models.ginan_freq_mhz`; B c2:49 | yes |
| Coverage floor 0.5 | `sed_prototype.py` `MIN_COVERAGE`; `beam_matching.weight_floor`; E c2:24, F c2:17, G c2:18, K2 c64:12 (`COV_MIN`) | yes |
| σ × √(χ²/ν) when χ²/ν > 1 (`INFLATE`) | `pipeline.york`; E c2:27, G c2:19, H c16:1, I c3:50; PP c1:41 | yes |
| ν₀ = 408 MHz | `sed_prototype.py`; `model_comp_config.yaml` `models.nu0_mhz`; I c3:47, L c1:19, PP c22:13 | yes |
| T–T pair cut `MIN_RATIO` 1.0 | E c2:25, PP c1:42 | yes |
| T_CMB | 2.72548 K in `sky_maps.py`; 2.725 K inside pygdsm, used for GINAN | **no** (≤ 22 mK on the GINAN map) |
| Survey uncertainties | `sed_prototype.py` `PUBLISHED_ERRORS` and `configs/error_bars.csv` | yes (two copies) |
| York fit | `pipeline.york`; copies in Parts E, G, H, I and PP cells 12, 22 | yes |

## 7. Stored notebook outputs

The notebooks are committed with their outputs. Last runs: Parts B–G, I–L on 2026-11-17 or 2026-11-18; Part K2 on
2026-11-18. Older: Part A (several cells have no stored output), Part H (its outputs are from the six-survey subset),
and `diagnostics.ipynb` (Part A shows the TRIS maps of an nside 8 run; Part C runs on earlier grids). Each Part's first
code cell records the start time and its last cell prints the run time.

## 8. Open items

**Method**

1. **Correlated errors.** Survey scale and zero-level terms, and TRIS's zero level, are the same for every point of a
   survey but are added per point as independent noise. Consequences: χ²/dof is not a goodness-of-fit test; the
   summed ln Z of step 3 counts each independent beam about six times (it ranks the models, but overstates the
   preference; thinning to one sample per beam gives Δln Z ≈ 5–7 instead of 36 for the curved power law); and step 7's
   RHINO-beam σ is a lower bound (0.4 % quoted, 0.8–0.9 % with the ring errors correlated, 3.6–4.0 % with the prior
   correlated: `rhino_checks.csv`). Part H fits one offset per survey as a partial test.
2. **Stage 1 at nside 16** (`tris.nside_hires: 16`): reduced χ² 16.1 / 5.9 against 3.09 / 3.22 with 64. Kept for
   speed and so that every product shares the TRIS grid.
3. **The stage 1 convergence cross-check is off** (`tris.solver.cross_check: false`). Switch it on before using the
   TRIS maps quantitatively; the TOD route does not use the maps.
4. **The 820 MHz zero level is asymmetric** (+0.43 / −0.30 K). Stage 1 uses one symmetric 0.5 K prior; the ring
   errors in the fits use the larger side.
5. **Beam matching of broad-beam surveys** (`scaled_cuts`, natives above ~1.2°: LWA1, Maipu): the matched beam comes
   out slightly wider than the target (LWA1 35 MHz +2.4 % / +1.4 %). `elliptical_gaussian` hits the widths but has
   Gaussian shoulders.
6. **Figure 9 reads two tables nothing in the repo writes**: `paper_outputs/nb08_paper_plots/fig15_part1_*.csv` and
   `fig15_part3_*.csv` (the β grid and reference points, and a check of the TRIS-only fits). The calculation should
   move into `10_paper_plots` cell 22 so the files are produced here. The rest of `nb08_paper_plots/` is not read by
   anything.
7. **Maipu 45 MHz:** Monsalve et al. 2021's corrections (×1.076, −160 K) are not applied; only their uncertainties are
   used. Maipu carries much of the low-frequency weight.
8. **Judgement calls in the survey uncertainties:** LWA1 20 % is the middle of 15–25 % (Dowell et al. 2017 also publish
   per-pixel uncertainty maps, not used); CHIPASS borrows Stockert's scale and zero-level errors; Rhodes's are upper
   limits used as 1σ.
9. **The Irfan correction:** the licence of the β_s / c_s maps; their full model also adds a free-free template, not
   added here; the maps' frame and ordering (Galactic, RING) were determined from the data, not stated in the paper.
10. **GINAN anchor frequency (to be confirmed):** the GINAN-corrected GSM2016 is 1.534 × GSM2016 plus an offset;
    anchoring at 300 MHz (F = 1.40) instead of 350 MHz would reduce the bumps where the ring crosses the Galactic plane
    by about a quarter (multiplicative term only; the offset term was not recomputed).
11. **TRIS beam width:** the archive header says the 3 dB beam is 18° wide in the E plane; the code uses 19.155° from
    the archive's beam cuts (checked in Part K1).
12. **The TRIS effective frequencies** (600.5, 817.8 MHz) and bandwidths come from limTOD's table
    (`../limTOD/limTOD/tris/archive.py`), attributed there to Zannoni et al. 2008 Sect. 2 and Tab. 1; not checked against
    the paper here. A 2 MHz error at 600 MHz shifts a synchrotron index by about 0.03.
13. **The 2.5 GHz TRIS channel** (six points, not a ring) is not used.

**Housekeeping**

14. **Regenerate the tracked step 2 and step 7 tables** before quoting them: `python run_pipeline.py --step ring_tods
    --no-report` (~20 min), then `--step rhino --no-report`, then `--step report`. The committed
    `rhino_tods_*.csv`, `rhino_summary.csv`, `rhino_checks.csv`, `ring_tods.csv`, `beam_matched_tods*.csv` and step 2's
    `manifest.json` predate the code (the Irfan method, the split columns and step 7's checks); `PIPELINE_REPORT.md`
    predates its tables. Then rerun `10_paper_plots`.
15. `paper_outputs/model_comparison/results/` is the standalone script's output from 2026-10-29 (4 data sets); rerun
    `run_comparison.py` or ignore it.
16. `T_CMB` is defined twice (§6).
17. A colour-vision check of the figures (protanopia, deuteranopia, tritanopia, greyscale) has been done with
    scripts outside the repo; there is no checker in the repo.
18. Part A and Part H have not been rerun since the pipeline's 2026-11-18 step 2 run (their inputs did not change).
19. The end-to-end tests need the TRIS archive, the FITS files and `pygdsm`; without them they skip, so CI on a bare
    machine tests only the algebra.

## 9. Where outputs go

`paper_outputs/` = what a paper figure, table or number is made from; `test_outputs/` = tests and checks the paper
does not use; `outputs/` = step 0 / stage products; the paper repository's `figures/` = the figures.

| Folder | What | Written by |
|---|---|---|
| `paper_outputs/pipeline/step1_beam_matched_maps/` | every allowed map beam-matched to the TRIS beam (nside 16) | step 1 |
| `paper_outputs/pipeline/step2_ring_tods/` | ring TODs and σ of every data set (`ring_tods.npz`); every survey and sky model through the TRIS beam (`beam_matched_tods.npz`, `sky_tods.npz`); the Irfan split (`irfan_split_tris.csv`) | step 2 |
| `paper_outputs/pipeline/step3_model_comparison/` | ln Z per ring sample and model, TRIS excluded | step 3 |
| `paper_outputs/pipeline/step4_model_selection/` | `selection.json`: the selected model | step 4 |
| `paper_outputs/pipeline/step5_fit/` | our spectral model: the selected model's UltraNest fit at every sample, TRIS included | step 5 |
| `paper_outputs/pipeline/step6_predict_408/` | the 408 MHz prediction against Haslam and the 408 MHz models | step 6 |
| `paper_outputs/pipeline/step7_rhino/` | the RHINO prediction through both beams, its checks and the Irfan split | step 7 |
| `paper_outputs/pipeline/PIPELINE_REPORT.md` | the last run's results, step by step | `--step report` |
| `paper_outputs/model_comparison/diagnostic_plots/` | diagnostics of the step 3 comparison | Part C |
| `paper_outputs/model_choice/` | the model-choice decision sheet | Part D |
| `paper_outputs/phil_plots/` | T–T β; sky models against TRIS (`plot6a_models_fitted_to_tris.csv` = `tab:step1`) | Part E |
| `paper_outputs/zero_level_offsets/` | zero-level offsets as fit parameters | Part H |
| `paper_outputs/robustness_check/` | robustness of the Haslam gain and offset | Part I |
| `paper_outputs/fisher_share/` | Fisher information on ln A(408) | Part L |
| `paper_outputs/nb08_paper_plots/` | tables read by Figure 9 (open item 6) | not written by the repo |
| `paper_outputs/model_comparison/results/` | the standalone comparison script's output | `scripts/model_comparison/run_comparison.py` |
| `test_outputs/beam_matched_vs_old/` | step 2's TODs against native-resolution TODs | Part A |
| `test_outputs/sky_models_408/` | curved power-law fits; sky models against both fits | Part B |
| `test_outputs/plot3_haslam_vs_fit/` | figures of step 6's scale–offset fits | Part F |
| `test_outputs/plot3_haslam_calibration/` | Haslam with one gain and offset | Part G |
| `test_outputs/tod_prediction/` | 13-point fits and predicted TODs (Tests 1–5) | Part J |
| `test_outputs/beam_variants/`, `fwhm_scaling_test/`, `beam_averaging_blobs/`, `lwa35_stockert/` | beam tests | Part K2; the beam-test scripts |
| `outputs/stage1*/` | TRIS maps (`outputs/stage1/` is step 0's; the others are variant runs) | step 0 / `--stage 1` |
| `outputs/stage2/` | Haslam and ARCADE 2 beam-matched on their own | `--stage 2` |
