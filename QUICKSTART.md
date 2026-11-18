# Handover

What someone taking this project over needs: which notebook cell makes which paper figure, table and number; the
environment; the data and what may not be shared; constants defined in more than one place; notebook previews that
are out of date; and the open items. As of **2026-11-16**. Line numbers in `main.tex` refer to the paper repo
(`Papers_and_reports/moonrunes_paper/main.tex`) as of that date.

Notebook overview: `NOTEBOOK_INDEX.md`. Every plot: `PLOT_INDEX.md`. Every dataset: `DATA_SOURCES.md`. Pipeline:
`PIPELINE.md`. Rules for working in this repo: `CLAUDE.md`.

---

## 1. Reproducing the paper, in order

1. **The pipeline:** `python run_pipeline.py --step all` (steps 0–7 and the report; `PIPELINE.md` Part I). Step 0
   writes the TRIS maps to `outputs/stage1/`; steps 1–7 write `paper_outputs/pipeline/step*/`; the report is
   `paper_outputs/pipeline/PIPELINE_REPORT.md`. It ran end to end on 2026-11-16 (18:17–19:59), with the LWA1
   38–74 MHz maps included (11 points per sample in step 3, 13 in step 5). Settings: `configs/run_pipeline.yaml` and
   `configs/model_comp_config.yaml` (`PIPELINE.md` P4; no setting is in both).
2. **`notebooks/tests_analysis/`**, Parts C–I, K and L, one notebook per Part (they read the pipeline's products; A, B and J are checks). `notebooks/test_results.ipynb` is the merged notebook they were split from (record only).
   Parts C–I make the paper's remaining numbers (`tab:step1`, Section 5's counts and residuals, zero levels,
   robustness).
3. **`notebooks/10_paper_plots.ipynb`** → the paper repo's `figures/` folder. Do not run the last line of cell 25
   (`rhino_tods_tris_beam`, not in the paper) unless that file is wanted.
4. `notebooks/diagnostics.ipynb` is independent (checks of the inputs and of stage 1).

`scripts/model_comparison/run_comparison.py` is not part of the chain: the comparison is pipeline step 3. It is
a standalone tool (the same two configs and data sets, plus the map route) writing to
`paper_outputs/model_comparison/results/`.

## 2. Environment

`environment.yml` describes the conda environment **`limtod`** (Python 3.14, healpy 1.20.0, matplotlib 3.11.1,
pygdsm 1.7.1, ultranest 4.6.3, joblib 1.6.0 for the parallel fits of `10_paper_plots` Figure 7). limTOD is imported from the local checkout `../limTOD` (`configs/run_pipeline.yaml`
`paths.limtod`), which carries local changes; do not `pip install limTOD`. The TRIS archive files are read from
`../limTOD/downloads/TRIS` (`paths.tris_archive_dir`). Notebooks run from `notebooks/`.

## 3. Paper figures

All in `notebooks/10_paper_plots.ipynb` (cell numbers 0-based), written to the paper repo's `figures/` as PDF + PNG.
Six-point set = TRIS 600, TRIS 820, LWA1 35, Maipu 45, LWA1 80, Stockert 1420. Ring = 120 samples, RA 0–357° in 3°
steps, Dec +42.433°; beam = `tris_cuts` 19.155° × 23.366° (E plane 7° east of the meridian), coverage ≥ 0.5.

| File | Label | Cell(s) | Inputs and configuration |
|---|---|---|---|
| `fig1_data_maps_tris_slice` | `fig:maps-slice` | 5 | five maps (`SURVEY_FILES`, cell 1 lines 37–38: LWA1 35, Maipu 45, LWA1 80, Haslam, Stockert), native resolution, no beam. Not the LWA1 38–74 MHz maps used by Figure 10 |
| `fig2_method_tris_beam` | `fig:tris-beam` | 20 | `beam_spec_from_config(kernel="tris_cuts")`; asserts 19.155/23.366 |
| `fig3_method_tris_maps` | `fig:tris-maps` | 7 | stage-1 maps (map route, circular 23.366° kernel in stage 2); stage-1 zero level **not added** |
| `fig4_method_models_vs_tris_tods` | `fig:tods-600-820` | 10 | TRIS ring data; model TODs from pipeline step 2 (`beam_matched_tods.npz`: Carter, Irfan at 408 MHz and GINAN at 350 MHz; Carter and GINAN carried with the base model's spectrum, Irfan with Irfan & Puglisi's β_s and c_s maps) |
| `fig5_tt_beta_matrix` | `fig:beta-matrix-ring` | 12 (York), 14 | 7 ring TODs (six-point set + Haslam); York, σ × √(χ²/ν); pairs need ≥ 20 samples |
| `fig6_tt_beta_vs_frequency` | `fig:beta-vs-freq` | 16 | as fig5 |
| `fig7_modelcomp_fits_vs_tris` | `fig:modelcomp-tris` | 18 | its own maximum-likelihood refits (24 starts), six-point set, priors from `configs/model_comp_config.yaml` |
| `fig8_skymodel_models_vs_prediction_408` | `fig:models-vs-prediction-408` | 28 | prediction = pipeline step 5 (our model: UltraNest maximum-likelihood A and posterior σ; 13 points with TRIS since the 2026-11-16 run, six in the paper's version); models from step 2 |
| `fig9_skymodel_tris_anchor_gz` | `fig:tris-anchor` | 22–23 | Haslam = g·P + z (York), 103 samples (Cygnus A ± 25° out); data sets, data and errors from step 5's fit files (every map `sed.exclude_maps` leaves in with data on the ring, plus TRIS: 13 points); cell 22 checks its TRIS-only fits against `paper_outputs/nb08_paper_plots/fig15_*.csv` (open item 7) |
| `fig10_rhino_tods_rhino_beam` | `fig:rhino-rhino-beam` | 25 (function), 26 | `paper_outputs/pipeline/step7_rhino/rhino_tods_rhino_beam.csv`: **our model (step 5, UltraNest; 13 points with TRIS in the 2026-11-16 run)**, Wiener deconvolution with a GSM2016 prior, RHINO beam from `../rhino-beams/`; GSM2008 + Irfan carried with Irfan & Puglisi's β_s and c_s (stale in the 2026-11-16 products until step 7 is rerun). Step 7 has run; **the figure in the paper is the earlier 13-point weighted-least-squares version and needs redrawing** |

## 4. Paper tables and quoted numbers

The tables are typed into `main.tex` by hand (the paper's `tables/` folder is empty). Cell numbers below are those of
`notebooks/test_results.ipynb` unless stated; in `notebooks/tests_analysis/`, subtract the Part's offset (table in its `README.md`).

| Paper item | Made by | Configuration to state |
|---|---|---|
| `tab:surveys` (main.tex:140) | by hand; matches `DATA_SOURCES.md:173–181` and `src/moonrunes/sed_prototype.py:47–122` | lists LWA1 35 and 80 only (Figure 10 uses nine LWA1 maps) |
| `tab:step1` (main.tex:363) | Part E, cell 148 → `paper_outputs/phil_plots/plot6a_models_fitted_to_tris.csv` | **TRIS = a·model + b**, weighted least squares, errors on TRIS only, covariance × χ²/ν; TRIS data from Part B's cache, models from Part A |
| `tab:408-residuals` (main.tex:541) | pipeline step 6 → `paper_outputs/pipeline/step6_predict_408/residuals_408.csv` (also printed by `10_paper_plots` cell 28; Part B's check: `test_outputs/sky_models_408/chi2_sky_models_vs_fits.csv`) | χ² with σ = prediction posterior σ only (no model error). **The % columns are in no file in the repo** |
| `tab:scale-offset` (main.tex:570) | pipeline step 6 → `paper_outputs/pipeline/step6_predict_408/scale_offset_408.csv` (the same York fit as Part F; its table: `test_outputs/plot3_haslam_vs_fit/`) | **model = a·P + b**, York, errors on both axes, σ inflated; every model given Haslam's σ (10% ⊕ 3 K); Cygnus A ± 25° → 103 samples |
| Section 5 verdict "0/114/6", "68 borderline" (main.tex:501) | Part C, cells 44 and 80, from pipeline step 3 (`paper_outputs/pipeline/step3_model_comparison/`; **the 2026-11-16 run has 11 points, the paper's counts came from the 4-point 2026-10-29 run: rerun Part C**) | **four points, TRIS excluded by design** (open item 1) |
| Stacked residuals "TRIS 817.8 −1.46σ, LWA1 80 +0.80σ" (main.tex:503) | Part D, cell 110 | six-point UltraNest fits from Part B (open item 2) |
| 408 MHz prediction "median 32.13 K, σ 0.17 K" (main.tex:523) | pipeline step 6 (`prediction_408.csv`); also `10_paper_plots` cell 28 | our model (step 5). **The 2026-11-16 run (13 points) gives 32.15 K, σ 0.168 K**; the paper's value is from six points |
| Zero levels as fitted offsets (main.tex:394) | Part H, cells 206–208 (Tables 1–2) | six-point set |
| g, z robustness (main.tex:396; by keyword also 490, 507–508, 607) | Part I, cell 228 → `paper_outputs/robustness_check/robustness_table.csv` | six-point baseline, 103 samples |
| Ratio "0.81 … 0.97" (main.tex:647) | Part F, Summary (cell 175) | the notebook says **0.98** (open item 5) |
| Section 7 numbers (main.tex:601–617, 628, 671, 689) | pipeline step 7 → `paper_outputs/pipeline/step7_rhino/rhino_summary.csv`; `10_paper_plots` cells 25–26 print them | **our model (step 5); the paper's numbers come from the earlier 13-point weighted-least-squares fit and are stale** |
| Beam tests (main.tex:271–275) | `notebooks/tests_analysis/` Parts K3 (`scripts/beam_averaging_blobs/`) and K4 (`scripts/fwhm_scaling_test.py`), shown from saved results | — |
| Zero levels +1.74 and +0.52 K (main.tex:303; Figure 3 caption) | stage 1 (pipeline step 0) → `outputs/stage1/stage1_manifest.json`, `frequencies.600.zero_level_k` = 1.7369 K and `frequencies.820.zero_level_k` = 0.5158 K; also printed by `10_paper_plots` cell 7 | the 2026-11-16 product: nside 16, GSM2008 prior; fitted, not added to the maps; `outputs/stage1/` is git-ignored, so step 0 rebuilds it |
| 23 and 7 of 1272 pixels narrowed by more than 5 per cent (main.tex:303) | `stage1_manifest.json`, `frequencies.<600, 820>.pixels_shrunk_gt_5pct` = 23, 7; `n_sky_pixels` = 1272 | posterior against prior σ per pixel; same product |
| Median posterior-to-prior width 0.998 and 0.999 (main.tex:303) | `stage1_manifest.json`, `frequencies.<600, 820>.sigma_shrinkage_median` = 0.9985, 0.9987 | same product |
| Stripe Dec +7.8° to +77.2° (Figure 3 caption; Figure 1's stripe) | `10_paper_plots` cell 1 (`DEC_LO`, `DEC_HI`: the stage 1 band mask's pixel-centre declinations, +9.59° to +75.34°, ± half an nside 16 pixel); its stored output prints it, and it was recomputed from `outputs/stage1/tris_maps_tris_haslam_v0.npz` | band = the 1272 pixels where the beam response exceeds 0.01 of its peak |
| Fisher shares and σ(ln A) leave-out (main.tex:625: 88, 79–95 and 11 per cent; 0.006, 0.065, 0.25; < 2e-5) and g at the median index (main.tex:625–626) | `notebooks/tests_analysis/` Part L, cells 11 (each quoted number with its definition, six-survey and 13-point values) and 13 (g) → `paper_outputs/fisher_share/` | Fisher (ln A, β) of T = A (ν/408)^β at step 5's 13-point fit (UltraNest maximum-likelihood point), errors as now, 103 samples (Cygnus A ± 25° masked); definitions of the old notebook 08 cell 54; the quoted values are the six-survey ones, which Part L reproduces |

## 5. Data, and what may not be shared

Every dataset, its source and its uncertainty: `DATA_SOURCES.md`. Do **not** copy or expose:

- `res/Juno_600MHz.fits`, `res/Juno_1248MHz.fits`: no provenance in the repo (header has only NSIDE/ORDERING; no unit,
  frame, beam or citation). Excluded from every fit (`sed.exclude_maps`).
- `external/irfan_puglisi_2025/*.npy`: third-party maps (Irfan & Puglisi 2025), git-ignored; licence not in the repo
  (to be cleared with the authors). Used by the pipeline's Irfan correction off 408 MHz (steps 2 and 7,
  `pipeline.sky_models.irfan_beta_map`, `irfan_cs_map`) and by Part J's β comparison.
- Anything derived from `../rhino-beams/` (RHINO beam and Jodrell Bank horizon mask): outside the repo; used by Part J
  Test 6 and so by Figure 10.
- `log_28thoct_allmaps_3rdmomentadded.txt`: contains local absolute paths.

`res/CHIPASS_1400MHz.fits` is git-ignored (`.gitignore:249`). The other maps in `res/` are LAMBDA products.

## 6. Constants defined in more than one place

Not merged; listed so that a change is made everywhere. "TR" = `notebooks/test_results.ipynb` cell numbering (in `notebooks/tests_analysis/`, subtract the Part's offset: `README.md` there), "PP" =
`notebooks/10_paper_plots.ipynb`, "c:L" = cell:line.

| Constant | Where | Agree? |
|---|---|---|
| Cygnus A RA 299.9° | TR A c3:32, B c15:51, C c34:55, D c95:22, F c162:17, G c178:18, H c196:33; PP c1:44 | yes |
| Cygnus half-window | 25° in TR F c162:17, G c178:18, PP c22:14; Parts C and D also use ±15° | **no** (two windows) |
| Carter 1.029, 0.91 K | TR A c3:27, J c235:49, B c15 (`HASLAM_CORRECTIONS`) | yes |
| Irfan monopole 8.9 K | `run_pipeline.yaml` `pipeline.sky_models.irfan_monopole_k`; TR A c3:28, B c15:32, J c235:50 | yes |
| Irfan ν₀ = 45 MHz and the β_s / c_s files | `run_pipeline.yaml` `pipeline.sky_models.irfan_*`; TR J c243:2–3 (`IRFAN_NU0`, `IRFAN_FILES`) | yes (both read the files as Galactic, RING) |
| GINAN anchor 350 MHz | TR A c3:29, B c15:48 | yes |
| Coverage floor 0.5 | `sed_prototype.py:123` `MIN_COVERAGE`; `run_pipeline.yaml` `beam_matching.weight_floor`; TR E c124:22, F c162:15, G c178:16 (`COV_MIN`); diagnostics A c62:9 | yes |
| `INFLATE` (σ × √(χ²/ν)) | TR E c124:25, F c162:16, G c178:17, H c210:1, I c222:47; PP c1:42 | yes |
| ν0 = 408 MHz | `sed_prototype.py:36`, `configs/model_comp_config.yaml` `models.nu0_mhz`, TR I c222:44, PP c22:14; TR H and J read the config | yes |
| TRIS frequencies 600.5 / 817.8 | TR A c3:31 (hard-coded); computed from the ring files in TR E and PP c3 | yes |
| T–T pair cut `MIN_RATIO` | 1.0 in TR E c124:23 and PP c1:43 (the removed notebook 08 used 1.3) | yes now |
| T_CMB | 2.72548 K in `sky_maps.py:316`; 2.725 K in pygdsm, used for GINAN (TR A, B) | **no** (≤ 22 mK on the GINAN map) |
| Survey error bars | `sed_prototype.py:47–122` (`PUBLISHED_ERRORS`) and `configs/error_bars.csv` | yes (two copies) |

## 7. Stale notebook previews

Stored outputs that may not match the current code or files. "Stale" means the stored text or image predates a
change; it is not proof that a number is different.

- **`10_paper_plots.ipynb`:** run on 2026-11-16 before the save names became `figN_…`, so its "saved …" lines print the
  old file names. The figures in the paper repo are the `figN_` files.
- **`tests_analysis/` Parts and `diagnostics.ipynb`:** stored outputs are from the original notebooks, not
  re-run after the merge (Part L excepted). Part D cell 6 holds the error of the failed 2026-11-17 run (fixed since). Part K2's first code cell (a copy of the old setup cell) has no output.
- **Part B Summary:** χ²/N 134.0 / 134.2 (run of 2026-10-23) against 130.72 / 137.11 in its own CSV and in the paper.
- **Part J first cell:** says "4 survey points … and both TRIS frequencies"; `DATASETS` and its outputs have 13.
- **Diagnostics Part C and `tests_analysis` Part K2:** partly on the earlier grid (TODO N-1).
- **`notebooks/tests_analysis/` Parts A and B:** the cells that read the pipeline's TODs and 408 MHz sky-model TODs have no
  stored output. Part J has no Test 6 cells (RHINO is step 7). Other cells' stored outputs print the old paths.
- **`10_paper_plots.ipynb`:** its stored outputs predate the switch to the pipeline's products (paths in printed
  lines are the old ones); Figure 10's cells need step 7 to have run. Cells 18 (parallel fits), 22–23 (Figure 9 on the
  pipeline's data sets), 25 (RHINO line widths) and the run-time cells have changed since the stored outputs.

## 8. Open items

1. **Section 5's comparison excludes TRIS by design** (pipeline step 3, `comparison.exclude_tris: true` in
   `configs/model_comp_config.yaml`): 11 points per sample (LWA1 35, the seven LWA1 38–74 MHz maps, Maipu 45,
   LWA1 80, Stockert 1420; 4 in the 2026-10-29 run). main.tex:501 describes it as made on the six ring surveys
   including TRIS; the text should change, not the code.
2. **Part D combines two data sets** by design: criterion 1 is the step 3 evidence (no TRIS); criteria 2–5 use the
   six-point fits (TRIS included).
3. **Survey sets:** `sed.exclude_maps` does not list the LWA1 38–74 MHz maps, so the pipeline uses
   11 points (step 3) and 13 (step 5, and everything built on it: steps 6–7, and the `notebooks/tests_analysis/` Parts and
   `10_paper_plots` cells that read step 5). **The paper's Sections 5–7, `tab:408-residuals`, `tab:scale-offset`,
   Figure 8 and Figure 10 were made from the 4- and 6-point sets and are stale until regenerated** (`PIPELINE.md`
   P10 gives the 2026-11-16 run's headline numbers). `tab:surveys` and Figure 1 show the six-point surveys. Figure 7
   (`10_paper_plots` cell 18) and `notebooks/tests_analysis/` Part J choose their data sets themselves (six and 13 points).
4. **"a" means opposite things in `tab:step1` and `tab:scale-offset`,** and the two tables use different fitters and
   error models (see §4).
5. **main.tex:647 says 0.97; Part F says 0.98.** The "2–3 K after the best single gain and offset" in the same sentence
   may come from Part G; not confirmed.
6. **`tab:408-residuals` % columns** are not in any file in the repo.
7. **Figure 9's check reads `paper_outputs/nb08_paper_plots/fig15_*.csv`**, written by notebook 08 (archived in
   `../moonrunes_archive_2026-11-17/notebooks/08_update_plots.ipynb`; also in git at `a5c827c`). Nothing in the repo
   writes these files. Its Fisher calculation (cell 54) is `tests_analysis` Part L, which reproduces its printed
   values on the six surveys.
8. **GINAN as an outlier (to be confirmed):** the GINAN-corrected GSM2016 is 1.534 × GSM2016 plus an offset carried
   with GSM2016's spectrum; the multiplicative term spans 8.5 K along the ring at 600 MHz and puts the bumps where the
   ring crosses the Galactic plane (RA ≈ 69° and 300–321°). **To be confirmed:** anchoring at a lower frequency inside
   the correction's valid range (e.g. 300 MHz, F = 1.40) would reduce the bumps by about a quarter (multiplicative term
   only; the offset term was not recomputed).
9. **`CLAUDE.md` names a colour-vision checker `check_figures_cvd.ipynb`; no such notebook exists.** The checks so far
   were done with throwaway scripts outside the repo.
10. **Old notebook numbers left in runtime strings** (not changed, because they end up in printed output or figure
    text): `src/moonrunes/sed_prototype.py:117, 121` ("notebook 03's noise test" → `notebooks/tests_analysis/` Part K2),
    `scripts/fwhm_scaling_test.py:195, 276` and `scripts/beam_tests_lwa35_stockert1420/fwhm_scaling_test_maps.py:187`
    ("notebooks/03" → `notebooks/tests_analysis/` Part K2).
11. **Juno maps:** provenance, units and frame unknown (see §5).
12. **TRIS beam width:** the archive header (`res/TRIS_absolute_600.txt`, line 6) says the 3 dB beam is 18° wide in the E
    plane; the code uses 19.155° (checked in `notebooks/tests_analysis/` Part K1).
13. **Irfan correction off 408 MHz** follows Irfan & Puglisi 2025 (their eqs. 1, 2, 6, 11): β_s and c_s from their maps,
    ν₀ = 45 MHz (`pipeline.irfan_scaling`). Open: (a) the licence of `external/irfan_puglisi_2025/*.npy`; (b) their full
    sky model also adds a free-free template (§3.2), not added here; (c) the maps' frame and ordering are not stated in
    the paper (Galactic, RING determined from the data, see `DATA_SOURCES.md`); (d) the 2026-11-16 step 2 and step 7
    Irfan products, Figure 4, Figures 9–10, the step 7 numbers and `notebooks/tests_analysis/` Part E's Irfan plots (6a–6c) predate it; (e) `CLAUDE.md` line 25 says every
    correction is carried with the base model's spectrum.
14. **Folders moved on 2026-11-17** (by hand): docs to `docs/`; `model_comparison/`, `beam_averaging_blobs/`,
    `beam_tests_lwa35_stockert1420/` and `tests/` to `scripts/`. Code, `pyproject.toml` (`testpaths`) and the docs
    point to the new places; stored notebook outputs and archived files still print the old paths.

## 9. What is in the archive folder

`../moonrunes_archive_2026-11-16/` (outside the repo) holds what was taken out on 2026-11-16, at its original relative
path:

- `notebooks/`: the 14 notebooks merged into `diagnostics.ipynb` and `test_results.ipynb` (see `NOTEBOOK_INDEX.md`).

`../moonrunes_archive_2026-11-17/` holds the configs before they moved to `configs/` and `notebooks/08_update_plots.ipynb`
(its Fisher cells are ported to `tests_analysis` Part L).
- `src/moonrunes/stage3_assemble.py`, `stage4_calibration.py`, `stage5_validate.py`: empty placeholders.
- `src/moonrunes/sed_prototype.py.before_2026-11-16_cleanup` and `stage2_beam_matching.py.before_2026-11-16_cleanup`:
  the two modules before `apply_full_resolution` and `plot_beam_matched` were removed.

Notebooks 02, 05, 06, 07, 08 and 13 are not there; they are in git at `a5c827c`.

## 10. Where outputs go

`paper_outputs/` = what the paper is made from, `test_outputs/` = tests and checks the paper does not use, `outputs/` = pipeline stage products, paper repo `figures/` = the figures. Stored notebook outputs print the paths from before this layout (2026-11-16).

| Folder | What | Written by |
|---|---|---|
| `paper_outputs/pipeline/step1_beam_matched_maps/` | every allowed map beam-matched to the TRIS beam (nside 16) | pipeline step 1 |
| `paper_outputs/pipeline/step2_ring_tods/` | ring TODs and σ of every data set (`ring_tods.npz`); every survey and sky model through the TRIS beam (`beam_matched_tods.npz`, `sky_tods.npz`) | step 2 |
| `paper_outputs/pipeline/step3_model_comparison/` | ln Z per ring sample and model, TRIS excluded | step 3 |
| `paper_outputs/pipeline/step4_model_selection/` | `selection.json`: the best model | step 4 |
| `paper_outputs/pipeline/step5_fit/` | our spectral model: the best model's UltraNest fit at every ring sample, TRIS included | step 5 |
| `paper_outputs/pipeline/step6_predict_408/` | the 408 MHz prediction against Haslam and the 408 MHz models (`tab:408-residuals`, `tab:scale-offset`) | step 6 |
| `paper_outputs/pipeline/step7_rhino/` | the RHINO prediction through the TRIS and RHINO beams (Figure 10) | step 7 |
| `paper_outputs/pipeline/PIPELINE_REPORT.md` | what the last run found, step by step | `--step report` |
| `paper_outputs/model_comparison/diagnostic_plots/` | diagnostics of the step 3 comparison | `notebooks/tests_analysis/` Part C |
| `paper_outputs/model_choice/` | the model-choice decision sheet | Part D |
| `paper_outputs/phil_plots/` | T–T β, sky models against TRIS (`plot6a_models_fitted_to_tris.csv` = `tab:step1`) | Part E |
| `paper_outputs/zero_level_offsets/` | zero-level offsets as fit parameters | Part H |
| `paper_outputs/robustness_check/` | robustness of the Haslam gain and offset | Part I |
| `paper_outputs/model_comparison/results/` | the standalone `run_comparison.py` output (the 2026-10-29 run on 4 points) | `scripts/model_comparison/run_comparison.py` |
| `paper_outputs/nb08_paper_plots/` | the removed notebook 08's figures and tables; `fig15_*.csv` is what Figure 9's cell checks against | notebook 08 (removed) |
| `test_outputs/beam_matched_vs_old/` | step 2's TODs against the old (native-beam) TODs | Part A |
| `test_outputs/sky_models_408/` | curved power-law fits, sky models against both fits (`chi2_sky_models_vs_fits.csv`) | Part B |
| `test_outputs/plot3_haslam_vs_fit/` | the scale-offset fits and figures (a check of step 6) | Part F |
| `test_outputs/plot3_haslam_calibration/` | Haslam with one gain and offset | Part G |
| `test_outputs/tod_prediction/` | 13-point fits and predicted TODs (Tests 1–5); Test 6's last outputs (2026-11-05) | Part J |
| `test_outputs/beam_variants/`, `fwhm_scaling_test/`, `scripts/beam_averaging_blobs/<map>/`, `lwa35_stockert/` | beam tests | Part K2; the beam-test scripts (Parts K3–K5) |
| `outputs/stage1*/` | TRIS maps (the pipeline's step 0 writes `outputs/stage1/`) | step 0 / `--stage 1` |
| `outputs/stage2/` | the older Haslam + ARCADE 2 beam matching (nothing reads it) | `--stage 2` |
| paper repo `figures/` | the paper figures `fig1_…`–`fig10_…` | `10_paper_plots.ipynb` |
