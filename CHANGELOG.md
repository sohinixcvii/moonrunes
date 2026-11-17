# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Stages refer to `tris_haslam_pipeline.md`, and the blockers to its "three blockers to
resolve before running anything".

## [Unreleased]

### Changed

#### `external/bayesian_skymap` no longer used (2026-11-17)

The submodule (a Gibbs recalibration, a preliminary design that was not adopted; kept in an older branch) is no
longer imported or referred to. No map, TOD, fit or number changes.
- `src/moonrunes/stage1_tris_maps.py`: `import_bayesian_func` and `_check_beam_against_bayesian_skymap` removed. The new
  `_beam_record` keeps the check that `beam.beam_deg` is the E/H mean and records `beam_deg`, `nside_new_pinned` and
  `pixels_per_fwhm`. The stage 1 manifest's `blocker_1_beam_reconciliation` loses `matches_bayesian_skymap`,
  `nside_new_from_bayesian_func` and `bayesian_func`. Stage 1 no longer needs the submodule to run.
- `scripts/tests/test_stage_io.py`: the nside test and the preconditioner-diagonal test, which compared against
  `bayesian_func`, removed (73 tests).
- `configs/run_pipeline.yaml`: `paths.bayesian_skymap` removed; the Gibbs-design comments removed (the superseded
  `haslam`, `assemble` and `gibbs` keys stay; only `gibbs.use_woodbury` is read, by a test).
- `pyproject.toml`, `environment.yml`: `emcee` removed (nothing imports it); scipy kept (used throughout), comments
  updated.
- `notebooks/diagnostics.ipynb`: the beam check in cell 13 and the mentions in cells 4, 12, 18, 93, 95, 110 and 111
  removed (sources only).
- Docs: the live references removed (README, HANDOVER, `tris_haslam_pipeline.md`, `docs/TODO.md` including the
  "Stage 4 prerequisites (patches to bayesian_skymap)" section, `docs/DATA_SOURCES.md`, `docs/TRIS_MAP_PROVENANCE.md`).
  HANDOVER's provenance table gains the stage 1 numbers of main.tex:303 and the Figure 3 caption, checked against the
  stage 1 manifest and band mask.

#### `test_results` split into Part notebooks, debugged for 13 points; Part L (Fisher) added; notebook 08 archived (2026-11-17)

- **Split:** `notebooks/test_results.ipynb` → `notebooks/tests_analysis/Part_<letter>_<heading>.ipynb`, one per Part
  (A–K5), cells and stored outputs copied; `tests_analysis/README.md` has the overview, the Part K banner and each
  Part's cell offset (old cell − offset = new cell). `test_results.ipynb` stays as a record, not maintained. Each Part
  finds the repo root from `configs/run_pipeline.yaml`, records its start time in the first code cell and prints the
  total run time in a last cell.
- **The failure** (2026-11-17 run, Part D cell 99): `assert sorted(keys) == sorted(SURVEYS)` — step 5's fits have 13
  data sets since the 2026-11-16 run, Parts D, G, H and I hard-coded the six ring surveys (G, H and I would have failed
  the same way). They now take step 5's data sets from `moonrunes.pipeline.fit_data_sets` (new, with `data_label`).
  Part D: `N_PTS` from the data; Part G: the assert keeps "Haslam is never a fit point"; Part H: data sets without an
  entry get the "published" zero-level rule (the LWA1 38–74 MHz maps quote none, so no offset); Part H's per-survey
  figure grid is sized to the data sets (a 2 × 3 grid dropped seven surveys); Part I: subsets "all", one left out per
  data set, "without the LWA1 maps", the best-measured three, and notebook 08's six surveys, whose fit is now the
  check against notebook 08's g and z (it was the baseline). Part H/I row labels "six-survey …" → "all-survey …".
  Part E keeps its own T–T set.
- **Stale caches that would have been reused:** Part B's curved fits (six-point, 2026-11-03) were accepted for any
  data set; Part D's UltraNest cache checked only the TOD key. Both now require the data sets to match (Part D saves
  them with each run).
- **Efficiency, same results:** the UltraNest runs of Parts B, C and D, and the maximum-likelihood fits of Part C, run
  in parallel (joblib). Seeded runs seed themselves and random starting points are drawn in the serial order first, so
  the results are those of a serial run (checked bit for bit on the real models; Part C's `run_nested` was and is
  unseeded).
- **Part L — Fisher information in the 408 MHz prediction** (new): notebook 08 cell 54's Fisher matrix and shares,
  reusing `pipeline.load_ring_tods`, `pipeline.our_model` (step 5's 13-point fit), `sed_prototype.fit_power_law` and
  `pipeline.york`. It reproduces notebook 08's printed values on the six surveys, then gives the 13-point values of
  every number main.tex:625 quotes, and g at the median index. Writes `paper_outputs/fisher_share/`.
- **Notebook 08** (`notebooks/08_update_plots.ipynb`, put back as a reference) moved to
  `../moonrunes_archive_2026-11-17/notebooks/` after Part L's check passed.
- **Paths after the folders moved to `docs/` and `scripts/`:** `pipeline._mc_modules`, `run_comparison.py` (repo root),
  the four tests and the beam-test scripts (repo root one level up), `pyproject.toml` `testpaths`, the notebooks'
  `MC_DIR`, and the paths in the docs.

#### Irfan correction carried with Irfan & Puglisi's own spectrum; Figure 7 fits in parallel (2026-11-17)

- **Irfan correction off 408 MHz** (`pipeline.irfan_scaling`, `_sky_model_maps_tris`, step 7): the corrected 408 MHz
  map, (T_tot − 8.9 K) as their eq. 6 (unchanged), is now carried to other frequencies with Irfan & Puglisi 2025's
  synchrotron model, T(p, ν) = T(p, ν₀) (ν/ν₀)^(β_s(p) + c_s(p) ln(ν/ν₀)), ν₀ = 45 MHz (their eqs. 1, 2, 11), using their
  `external/irfan_puglisi_2025/final_beta_1deg.npy` and `final_cs_1deg.npy` (nside 256, RING, Galactic; frame checked
  against Haslam's small-scale structure, +0.23 Galactic vs +0.02 equatorial; brought to the model's nside with
  `hp.ud_grade`). Until now it used GSM2008's / GSM2016's own per-pixel spectrum. New config keys
  `pipeline.sky_models.irfan_nu0_mhz`, `irfan_beta_map`, `irfan_cs_map`; step 7's sky-model cache key includes them.
  **Numbers change** at the TRIS frequencies (step 2 `modeltris::*Irfan*`, Figure 4) and at 60–80 MHz (step 7,
  Figures 9–10) once steps 2 and 7 are rerun; 408 MHz values do not change. Licence of the two maps not yet cleared.
- **`10_paper_plots` cell 18 (Figure 7):** the maximum-likelihood fits run in parallel (joblib). The starting points are
  drawn in the serial order first, so the fits are bit-identical (checked on the real models).

#### Pipeline report tables; run times (2026-11-17)

- `PIPELINE_REPORT.md` gains plain-text tables for steps 3 (evidence of every model), 4 (selection criteria and each
  model against them), 6 (top 3 sky models against the 408 MHz prediction, by χ²) and 7 (top 3 per beam, by χ²/N over
  every frequency), read from each step's products (`pipeline.REPORT_TABLES`; `write_report` takes an optional `path`).
  No step's computation changed.
- `run_pipeline.py` prints the total run time at the end. `10_paper_plots.ipynb` (cells 1 and 29) and
  `diagnostics.ipynb` (cell 3 and a new last cell, 132) record the start time in the first code cell and print the
  total notebook run time in the last cell.

#### Notes, comments and docs: change history taken out of the running text (2026-11-17)

- Notebook Markdown, code comments, docstrings, configs and docs no longer narrate what was moved, merged, removed or
  renamed and when; that history is in this file. Banner labels in the notebooks: **Origin:**, **Stored outputs**,
  **Pipeline:**. Facts, paths and numbers are unchanged; no code changed. Fixed on the way: `model_comparison/README.md`
  gave `max_likelihood` as the default statistic (the code's default and the config's choice are `evidence`);
  `10_paper_plots` cells 27–28 and `PLOT_INDEX.md` called the 408 MHz prediction six-survey and its source notebook
  12's cache (it is pipeline step 5); `configs/run_pipeline.yaml` said the 2026-11-16 products were made without the
  LWA1 38–74 MHz maps (they include them).

#### `10_paper_plots`: the TRIS-anchor figure on the pipeline's data sets; RHINO figures less crowded (2026-11-17)

- **`fig9_skymodel_tris_anchor_gz` (cells 0, 21–23):** fits (b) without TRIS and (c) all surveys now use every map
  `sed.exclude_maps` leaves in with data on the ring, plus TRIS for (c): 11 maps + TRIS 600/820, read with their data and
  errors from step 5's fit files (the same values as step 2's TODs and `sigma_for`, checked). It used notebook 08's six
  fixed surveys and stopped with `KeyError: 'LWA1_38.fits'` once step 5 had 13. CHIPASS and Rhodes drop out (no data on
  the ring; asserted). The open marker "six-survey median β" is now (c)'s median β (−2.487). The TRIS-only fits (a) and
  the fixed-β line do not depend on the survey set and are still checked against notebook 08's CSVs. **Numbers change:**
  (b) g 1.3748 → 1.4374, z −9.875 → −15.015 K; (c) g 1.0548 → 1.0553, z −4.771 → −4.802 K; (a) unchanged.
- **`rhino_tods_tris_beam`, `fig10_rhino_tods_rhino_beam` (cell 25):** model lines 1.3 → 0.8 pt, markers 4 → 2.6,
  figure height 4.6 → 6.4 in. No data or colour changed.

#### Both configs in `configs/`, no setting in both (2026-11-17)

- **Moved:** `configs/run_config.yaml` → `configs/run_pipeline.yaml` (the run config, `load_config()`'s default);
  `model_comparison/config.yaml` → `configs/model_comp_config.yaml`; `model_comparison/error_bars.csv` →
  `configs/error_bars.csv` (`data.error_bars_csv` is now relative to `configs/`). Originals in
  `../moonrunes_archive_2026-11-17/`. New in `moonrunes.config`: `load_model_comp_config`, `model_comp_config_path`,
  `error_bars_path`, `DEFAULT_MODEL_COMP_CONFIG_PATH`.
- **Duplicates removed:** the Haslam target (`data.target_file` removed; `pipeline.target_file` is read by the
  pipeline, `run_comparison.py` and `test_results` Parts B, C); the comparison statistic (`pipeline.compare.statistic`
  removed; `comparison.statistic` is read by the pipeline too, so `test_results` Part C can no longer disagree with
  step 3); TRIS's exclusion (`pipeline.compare.exclude_tris` → `comparison.exclude_tris`, now honoured by
  `run_comparison.py` and Part C); `data.datasets` removed, so `run_comparison.py` and Part C use the maps
  `sed.exclude_maps` leaves in, like the pipeline (`run_comparison.select_datasets` gained `exclude_tris`; an
  explicit list still works for Part J). No value changed.
- **Consequence:** `run_comparison.py` and `test_results` Part C, when rerun, use the pipeline's 11 points (they used
  `data.datasets`' 4). `run_comparison.py` gained `--run-config`.
- **Docs:** the config names everywhere outside archives and the changelog's history; `PIPELINE.md` P4 (the two
  files), P6 (data sets: 11/13 points with the LWA1 38–74 MHz maps, removed from `sed.exclude_maps` on 2026-11-16),
  P10–P11 rewritten for the first full run (2026-11-16); HANDOVER §1, §3, §4, §8, §10; NOTEBOOK_INDEX; PLOT_INDEX;
  `model_comparison/README.md`; CLAUDE.md (paths only). The `sed.exclude_maps` comment in `run_pipeline.yaml` now
  says the LWA1 maps are in.

#### The end-to-end pipeline: `run_pipeline.py --step all` (2026-11-16)

- **New `src/moonrunes/pipeline.py`, `run_pipeline.py --step <name> | all | report` (`--from-step`):** 0 `tris_maps`
  (stage 1), 1 `beam_match` (every allowed map → TRIS beam, `tris_cuts`, nside 16), 2 `ring_tods` (120-sample TODs and
  σ of every allowed data set and of Haslam; the sky models at 408 MHz and the TRIS frequencies), 3 `compare` (every
  model of `model_comparison/config.yaml`, TRIS excluded, UltraNest evidence), 4 `select` (summed ln Z, ties within
  `tie_lnz` to the simpler model), 5 `fit` (the best model with TRIS, UltraNest posterior: our spectral model),
  6 `predict_408` (against Haslam and the 408 MHz models: residuals, York scale-offset), 7 `rhino` (TRIS beam →
  Wiener deconvolution → RHINO beam), and a report (`paper_outputs/pipeline/PIPELINE_REPORT.md`) from each step's
  `manifest.json`. `--stage 1|2` still work. New `pipeline:` section in `configs/run_config.yaml`.
- **Code moved out of `notebooks/test_results.ipynb`:** Part A's TOD build (now step 2), Part B's power-law fits and
  408 MHz sky-model TODs (steps 5, 2), Part J's Test 6 (step 7). Those Parts now read the pipeline's products; Parts
  B (curved power law), F and J (Tests 1–5) are checks. Parts C–I, K2 and `10_paper_plots.ipynb` (cells 1, 22, 25,
  28) read `paper_outputs/pipeline/`. The notebooks as they were: `../moonrunes_archive_2026-11-16/notebooks/*.before_pipeline`.
- **Outputs moved:** `paper_outputs/beam_matched_tods/` → `pipeline/step2_ring_tods/` (TODs) and
  `test_outputs/beam_matched_vs_old/` (the comparison); `paper_outputs/sky_models_408/` → power-law fits to
  `pipeline/step5_fit/fits/`, `sky_tods.npz` to `step2_ring_tods/`, the rest to `test_outputs/sky_models_408/`;
  `paper_outputs/tod_prediction/` and `plot3_haslam_vs_fit/` → `test_outputs/`; the run_comparison.py results copied to
  `pipeline/step3_model_comparison/`; Part J's RHINO sky-model cache copied to `pipeline/step7_rhino/cache/`.
- **Not run.** Steps 2–5 are **seeded** with the existing products they reproduce (their `manifest.json` says
  `"seeded": true` and from where; step 4's `selection.json` applies the rule to the seeded step 3 numbers: the power
  law). Steps 1, 6, 7 have no products yet. **Running step 7 will change Figure 10 and the Section 7 numbers** (our
  six-point UltraNest model instead of Part J's 13-point weighted least squares).
- **Docs:** `PIPELINE.md` rewritten (Part I: the pipeline, with a flowchart; Part II: the earlier stage reference and
  Appendices A–F, numbering unchanged; old §8, §10, §11 dropped); README (Running, outputs table, Layout), HANDOVER (§1,
  §3, §4, §7, §8, §10), NOTEBOOK_INDEX, PLOT_INDEX (including the locations table, which still had pre-split paths),
  CLAUDE.md, `model_comparison/README.md`.

#### Outputs split into `paper_outputs/` and `test_outputs/` (2026-11-16)

- **Rule:** `paper_outputs/` = anything a paper figure, table or quoted number is made from; `test_outputs/` = tests and
  checks the paper does not use; `outputs/` keeps the pipeline stage products (`run.outputs_dir`, unchanged); figures
  are still written to the paper repo's `figures/`.
- **Moved (whole folders, caches included, nothing re-run):** `model_comparison/{beam_matched_tods, sky_models_408,
  model_choice, phil_plots, plot3_haslam_vs_fit, zero_level_offsets, robustness_check}` → `paper_outputs/`;
  `model_comparison/{results, diagnostic_plots}` → `paper_outputs/model_comparison/`; `tod_prediction/` →
  `paper_outputs/tod_prediction/`; `outputs/paper_plots/` → `paper_outputs/nb08_paper_plots/`;
  `model_comparison/plot3_haslam_calibration/` → `test_outputs/`; `outputs/{beam_variants, fwhm_scaling_test}` →
  `test_outputs/`; `beam_averaging_blobs/{HASLAM, HASLAM_DESTRIPED_ONLY}` → `test_outputs/beam_averaging_blobs/`;
  `beam_tests_lwa35_stockert1420/{fwhm_scaling, beam_averaging_blobs}` → `test_outputs/lwa35_stockert/`.
  1418 + 53 files; `model_comparison/` now holds only code, config, `error_bars.csv`, its README and `archive/`.
- **Paths updated:** path constants in `notebooks/test_results.ipynb` (Parts A–J, K2) and `notebooks/10_paper_plots.ipynb`
  (cells 1, 22, 25, 28); Markdown path mentions in both; `model_comparison/run_comparison.py` (output folder),
  `scripts/fwhm_scaling_test.py` (`--out` default), `beam_averaging_blobs/beam_averaging_blobs.py`, `replot_from_csv.py`,
  `beam_tests_lwa35_stockert1420/*.py` (output and input folders, and the provenance line they write into their CSVs);
  `mkdir(parents=True, ...)` where the target is now nested. No other code, value or number changed.
- **Docs:** README (Layout, new "Where outputs go"), `NOTEBOOK_INDEX.md`, `PLOT_INDEX.md`, `HANDOVER.md` (new §10),
  `PIPELINE.md`, `CLAUDE.md` (one rule), `model_comparison/README.md` and `config.yaml` comments, both beam READMEs
  (and their copies in `test_results.ipynb` Parts K3 and K5), `.gitignore` comment. Dated records (older entries here,
  `TODO.md`, `inputs_audit/inputs_audit.md`, `model_comparison/archive/`) keep the old paths.

#### Python files: dead code removed (2026-11-16)

- **Checked:** every `.py` file outside `external/` and the frozen `model_comparison/archive/` (24 files). Each top-level
  function, class and constant was looked up by name in every other `.py` file and in the code cells of the three
  notebooks; imports were checked with an AST pass.
- **Removed (unused anywhere):** `sed_prototype.apply_full_resolution` (its only caller was notebook 08, removed earlier
  the same day) and `stage2_beam_matching.plot_beam_matched` (never called); the unused `List` import in
  `sed_prototype.py`. The files as they were are in `../moonrunes_archive_2026-11-16/src/moonrunes/*.before_2026-11-16_cleanup`.
- **Moved to `../moonrunes_archive_2026-11-16/src/moonrunes/`:** the empty placeholders `stage3_assemble.py`,
  `stage4_calibration.py`, `stage5_validate.py` (0 lines, imported nowhere; `run_pipeline.py --stage 3/4/5` still says
  "not implemented yet").
- **Not consolidated:** no function is duplicated between files; the LWA1 35 / Stockert test scripts already import the
  Haslam scripts (`scripts/fwhm_scaling_test.py`, `beam_averaging_blobs/beam_averaging_blobs.py`) instead of copying
  them. Merging `beam_averaging_blobs/replot_from_csv.py` into `beam_averaging_blobs.py` would change its command line
  and could only be checked by regenerating its figures, so it was left.
- **Checks:** `pytest tests` 75 passed before and after; the edited modules import; no numbers or outputs changed.
- **Docs:** `PIPELINE.md` (stage table rows 3–5, line 779, the H-a row), `README.md` (Layout), `TODO.md` (N-11).

#### Notebooks merged into three: `diagnostics`, `test_results`, `10_paper_plots` (2026-11-16)

- **`notebooks/diagnostics.ipynb` (new):** Part A = `01_explore_data`, Part B = `04_coordinate_frames`, Part C =
  `03_diagnostics` cells 0–38 (stage-1 diagnostics and prior-only runs).
- **`notebooks/test_results.ipynb` (new), in dependency order:** A = `18_beam_matched_tods`, B = `12_sky_models_vs_fit_408`,
  C = `11_model_comparison_diagnostics`, D = `15_model_choice_decision_sheet`, E = `14_ring_tt_betas`,
  F = `16_haslam_vs_power_law_408`, G = `17_haslam_gain_offset_calibration`, H = `20_zero_level_offsets`,
  I = `21_robustness_check`, J = `19_power_law_tod_prediction`; K = beam tests: K1 = `09_beam_fwhm_tests`,
  K2 = `03_diagnostics` cells 39–114 (with a copy of its setup cell, cell 1, so the Part runs on its own), K3–K5 = the
  saved results of `beam_averaging_blobs/`, `scripts/fwhm_scaling_test.py` and `beam_tests_lwa35_stockert1420/`
  (their READMEs and figures; the scripts stay where they are).
- **`notebooks/10_paper_plots.ipynb`:** unchanged.
- **How:** every code cell copied with its source, stored outputs and execution count unchanged (checked: all 228
  original code cells present once); Markdown headings one level lower; one overview cell per notebook and one banner
  per Part (what it was, what it reads and writes, which paper items it feeds). **Nothing was re-run.** Both notebooks
  declare the `limtod` kernel (03 had `21cmfast`).
- **Moved, not deleted:** the 14 merged originals are in `../moonrunes_archive_2026-11-16/notebooks/`. Notebooks 02, 05,
  06, 07, 08 and 13 were removed earlier on 2026-11-16 (in git at `a5c827c`).
- **Docs:** `NOTEBOOK_INDEX.md` rewritten (with an old → new table); new `HANDOVER.md`; `README.md` (Notebooks, Standalone
  tests, Layout), `PLOT_INDEX.md` (made-by column, section headings, `figN_` paper file names, and the exclude-list
  sentence, which wrongly said LWA1 35 and Maipu 45 were excluded), `PIPELINE.md` (§10 and references),
  `DATA_SOURCES.md`, `TRIS_MAP_PROVENANCE.md`, `tris_haslam_pipeline.md`, `model_comparison/README.md`,
  `beam_tests_lwa35_stockert1420/README.md`; a pointer note at the top of `TODO.md` and `inputs_audit/inputs_audit.md`
  (which keep the old numbers as dated records). Comments and docstrings that named notebooks: `configs/run_config.yaml`,
  `environment.yml`, `.gitignore`, `src/moonrunes/{frames,elliptical_beam,sed_prototype,sky_maps,stage2_beam_matching}.py`,
  `tests/test_beam_fwhm.py`, `tests/test_elliptical_beam.py`, `scripts/fwhm_scaling_test.py`,
  `beam_tests_lwa35_stockert1420/fwhm_scaling_test_maps.py`. No code, value or number changed.
- **Not changed:** strings that end up in output (`sed_prototype.py:117, 121`; `scripts/fwhm_scaling_test.py:195, 276`;
  `beam_tests_lwa35_stockert1420/fwhm_scaling_test_maps.py:187`) still name notebook 03; older entries of this file and
  the frozen `model_comparison/archive/` keep the old numbers.

#### Paper figure: the 408 MHz models against the power-law prediction, notebook 10 Figure 11 (2026-11-06)

- **`notebooks/10_paper_plots.ipynb`, Figure 11 (cells 26–27), `figures/skymodel_models_vs_prediction_408.pdf/.png`.**
  The six-survey power-law prediction at 408 MHz along the ring (notebook 12's fit cache: maximum-likelihood A, σ the
  standard deviation of the posterior samples) against the beam-matched 408 MHz TODs of Haslam, its Carter and Irfan
  corrections, GSM2008, GSM2016, GSM2016 + GINAN and LFSM (notebook 18), through the TRIS beam; model − prediction
  below. In `main.tex` §6.2 (`fig:models-vs-prediction-408`). The medians reproduce the paper's 408 MHz residual table
  to within 0.05 K. Run once on its own (cell 1's settings only); the notebook was not rerun.
- **`rhino_tods_tris_beam` taken out of the paper** (`main.tex`); notebook 10 still makes it. `PLOT_INDEX.md`,
  `NOTEBOOK_INDEX.md` and `README.md` updated.

#### Paper figures renamed, and every figure placed in `main.tex` (2026-11-06)

- **New names, `<section>_<content>`,** the prefix being the paper section the figure sits in (`data`, `method`, `tt`,
  `modelcomp`, `skymodel`, `rhino`, as in the `sec:` labels): `fig1_all_maps_tris_slice` → `data_maps_tris_slice`,
  `tris_beam` → `method_tris_beam`, `tris_maps` → `method_tris_maps`, `plot6_tods_600_820_data_vs_models` →
  `method_models_vs_tris_tods`, `plot4f_beta_matrix_whole_ring` → `tt_beta_matrix`,
  `plot4c_beta_vs_mean_frequency_span` → `tt_beta_vs_frequency`, `model_comparison_fits_vs_tris` →
  `modelcomp_fits_vs_tris`, `tris_anchor_gz_plane` → `skymodel_tris_anchor_gz`. The two RHINO figures keep their
  names. Notebook 08's and 14's own working files (`outputs/paper_plots/`, `phil_plots/`) are not renamed.
- **`notebooks/10_paper_plots.ipynb`:** the `save()` calls and the "Featured in" notes use the new names and the paper
  sections the figures are now in (3.2, 6.3, 7.1, 7.2 for the four that were not yet in `main.tex`). `save()` still
  writes PDF + PNG; `main.tex` includes only the PDFs, and the old PNGs were removed from `figures/`. Not rerun.
- `PLOT_INDEX.md`, `NOTEBOOK_INDEX.md` and `README.md` updated to match.

#### RHINO figures in the paper plots, notebook 10 Figures 9–10 (2026-11-06)

- **`notebooks/10_paper_plots.ipynb`, Figures 9–10 (cells 23–25), `figures/rhino_tods_tris_beam.pdf/.png` and
  `figures/rhino_tods_rhino_beam.pdf/.png`.** Notebook 19 Test 6's figures at MNRAS full width, read from its CSVs in
  `tod_prediction/`, with Figure 3's model styles. Not yet in `main.tex`. Only cells 1 and 24–25 were run.

#### RHINO predictions, notebook 19 Test 6 (2026-11-05)

- **`notebooks/19_power_law_tod_prediction.ipynb`, Test 6 (cells 21–23; settings added to cell 2).** The power-law
  prediction at 60, 65, 70, 75 and 80 MHz through the TRIS beam, then deconvolved from it (Wiener solve on the nside 16
  grid, stage 1's prior width, GSM2016 as the prior) and observed with the RHINO
  `DryHornGround` CST beam at the Jodrell Bank zenith (`JodrellBank.fits` horizon, ground term left out). GSM2008,
  GSM2016, LFSM, GSM2008 + Carter and GSM2008 + Irfan (applied at 408 MHz, carried with GSM2008's per-pixel spectrum) and GSM2016 + GINAN
  (evaluated directly, valid 60–350 MHz) go through both beams. Writes `tod_prediction/test6_rhino_freqs_{tris,rhino}_beam.png/.csv`; the
  sky models are cached in `tod_prediction/cache/rhino_sky_models.npz`. The beam's azimuth orientation is limTOD's
  identity mapping, not checked against the as-built horn.

#### Robustness check of the Haslam gain and offset, notebook 21 (2026-11-05)

- **`notebooks/21_robustness_check.ipynb`, new** (a collaborator's comment 2). Haslam = g P + z (York, as notebook 08
  Figs 15–16) on a fixed grid: survey subsets (each survey left out, both LWA1 bands out, Maipu 45 + TRIS 600 + 820,
  TRIS only), power law and curved power law (≥ 4 surveys), Cygnus A mask ±15° / ±35°, the RA halves split at 200°,
  and three rows copied from notebook 20's Table 3(a). The baseline (six-survey power law, 103 samples) reproduces
  g = 1.0548 ± 0.0536, z = −4.77 ± 1.47 K. Writes `model_comparison/robustness_check/robustness_table.csv` and
  `robustness_gz_plane.png/.pdf`. No spectral fit failed the convergence check.

#### Paper figure: the TRIS anchor of the Haslam gain and offset (2026-11-05)

- **`notebooks/10_paper_plots.ipynb`, Figure 8 (cells 20–22), `figures/tris_anchor_gz_plane.pdf/.png`.** Notebook 08
  Fig 15 at MNRAS column width. The cell redoes notebook 08's fits (the ellipses need the g–z covariance, which 08 does
  not save) and asserts they match 08's CSVs in `outputs/paper_plots/`. Not yet in `main.tex`.

#### Fix: notebook 08 Figure 14b `NameError: PA` (2026-11-05)

- The 2026-11-02 trim of Figure 14 (cell 46) had removed the `PA` / `FW` arrays that Figure 14b (cell 48) uses, so a
  fresh run stopped there. The two lines are back in cell 46; Figure 14 is unchanged. Notebook 08 runs from the top.

#### Notebook 08 Figures 15–16: how strongly TRIS anchors the Haslam gain and offset (2026-11-04)

- **Cells 49–61** (a collaborator's comment 1): Haslam = g P + z for P from TRIS only, from the four other surveys
  and from all six; σ(ln A at 408 MHz) from the Fisher matrix with each survey left out; the TRIS-only fit against a
  fixed extrapolation β. Figures `fig15_tris_anchor_gz_plane`, `fig16_sigma_lnA408_leave_out`, tables
  `fig15_part*`, `fig16_part*` in `outputs/paper_plots/`.

#### Notebook 20: TRIS II Table 3 zero levels filled in (2026-11-03)

- `ZERO_LEVEL_SIGMA_K`: TRIS 600 = 0.066 K, TRIS 820 = 0.430 K (the larger side of +0.430 / −0.300 K; footnote a's
  0.659 K recorded, not used). They equal the zero-level terms already inside `ring_sigma`, so every number is
  unchanged; the text no longer says they are missing, and Figure A marks them.

#### Documentation (2026-11-03 to 2026-11-05)

- **`NOTEBOOK_INDEX.md`, new:** what each notebook does, what it reads and writes, and its status.
- `PLOT_INDEX.md`: notebooks 19–21, notebook 08 Figs 15–16 and the paper figure `tris_anchor_gz_plane` added.
- `README.md`: the notebook table covers 01–21 (notebook 08's name corrected to `08_update_plots`).

#### Not recorded here at the time

- Notebook 08 Figures 14 / 14b and notebook 10's `tris_beam` (2026-11-02), notebook 03's beam-variant section
  (`outputs/beam_variants/`), and notebooks 19 and 20 themselves were added without a changelog entry; see
  `NOTEBOOK_INDEX.md` and `PLOT_INDEX.md` for what they contain.

#### Model comparison: choose the data sets (2026-10-28)

- **`model_comparison/config.yaml`: `data.datasets`,** the list of data sets that go into every SED (keys as in
  `error_bars.csv`; TRIS bands included). It overrides `sed.exclude_maps` for the comparison only;
  `run_config.yaml` is not edited. `null` keeps the old behaviour. Set to the current 6: TRIS 600, TRIS 820,
  LWA1 35, Maipu 45, LWA1 80, Stockert 1420.
- **`run_comparison.py`:** new `select_datasets()`; `tod_units()` / `map_units()` take the TRIS bands to include
  (default both, so notebooks calling them are unchanged); summaries get a `datasets` column.

#### Model comparison without the Occam penalty (2026-10-28)

- **`model_comparison/config.yaml`: new `comparison` block, `statistic: max_likelihood`.** The models are now
  compared on each one's best fit (minimum χ² inside the prior box, ΔlnL_max = −Δχ²/2 against the power law),
  so no model is penalised for having more parameters. `statistic: evidence` gives the UltraNest ln Z
  comparison as before. `ml_starts: 24` optimiser starts per fit.
- **`model_comparison/run_comparison.py`:** new `max_likelihood()` (L-BFGS-B in the prior's unit cube, the same
  fitter as notebook 11 Part 5) and `write_maxlike()`; writes `results/maxlike_per_unit_<route>.csv` (χ², ln L_max,
  degrees of freedom, best-fit parameters per unit and model) and `results/maxlike_summary_<route>.csv` (medians
  and counts, no sums). The `evidence_*` files are not touched by a max-likelihood run. `models.py` and the
  priors are unchanged. Not yet run on the data.
- Why: the evidence's Occam cost for C (prior ±0.5, posterior about ±0.07) is about 1.7 in ln Z, which cancelled
  the curved power law's typical χ² gain of about 3 and gave ties on most of the ring.

#### Notebook 11 diagnoses the order-3 moment expansion (2026-10-28)

- **`notebooks/11_model_comparison_diagnostics.ipynb`: four models.** The order-3 moment expansion
  (`models.moment_orders: [2, 3]`) is now in every part except Part 3 (the C prior test): Figs 1–3 get an
  order-3 panel, Fig 4 two more comparisons (order 3 − power law, order 3 − order 2), Figs 5–6 and the Part 4
  convergence checks and Fig 9 include it. The exact grid integral uses 25 / 41 points per axis for the
  4-parameter model (41 / 81 for the others); `convergence_check1_exact.csv` now has the columns
  `grid points per axis`, `grid ln Z (coarse)`, `grid ln Z (fine)`.
- **Fix: `ImportError: cannot import name 'tod_method_key'` in cell 15.** A kernel that had imported
  `sed_prototype` before 2026-10-27 kept the old module. Cell 15 now drops the already-loaded `moonrunes`,
  `models` and `run_comparison` modules and imports them again.
- Its saved figures and the "Part 4 in plain numbers" / Summary text are from earlier runs until it is rerun.

#### More maps left out of the SED fits (2026-10-27)

- **`sed.exclude_maps`** in `configs/run_config.yaml` now also lists `LWA_35MHz.fits`, `MAIPU_45mhZ.fits`,
  `Juno_600MHz.fits` and `Juno_1248MHz.fits` (besides ARCADE 2). Each TOD-route ring sample now has 4 points
  (TRIS 600/820, LWA1 80, Stockert 1420), so the order-3 moment expansion has no degrees of freedom left.
  `model_comparison/results/` (2026-10-27 17:13) was made this way. The reason is not yet recorded.

#### Beam-matched ring TODs everywhere (2026-10-27)

- **`sed_prototype.Inputs.ring_tod`, new; `sed.tod_beam: beam_matched` in `run_config.yaml`.** Every ring TOD
  is now stage 2's beam matching (`tris_cuts` kernel, the map's own beam taken out) evaluated at the 120 ring
  pointings, as notebook 18 built them. Checked equal to notebook 18's cache to 0 K for every survey, Haslam
  and the 408 MHz sky models. The older route (native map observed with limTOD) stays as
  `sed.tod_beam: native` / `Inputs._native_tod`.
- **Users:** `run_comparison.py` (prints the TOD method, and writes it to the summary CSV as `tod_beam`),
  notebooks 07 (its own copy of the TOD code), 08 (sky-model TODs), 10, 11 and 12 (sky-model TODs, with the
  pygdsm models' native beams). Notebooks 13 and 15–17 read notebook 12's cache.
- **Caches** in notebooks 11, 12 and 15 now store a `tod_key` (`sed_prototype.tod_method_key`), and entries made
  from other TODs are redone. Notebook 15 stops if notebook 12's cache is stale.
- **Notebook 18** compares against a snapshot of the old TODs, `model_comparison/beam_matched_tods/old_native_tods.npz`.
- **Not changed:** the map route and stage 2 (`beam_matching.kernel: circular_gaussian`), and the beam studies
  (notebooks 01, 03, 09 and `beam_averaging_blobs/`).

### Added

#### Beam tests on LWA1 35 MHz and Stockert 1420 MHz (2026-10-28)

- 🧪 **`beam_tests_lwa35_stockert1420/`, new.** The FWHM-scaling test (two-step vs one-step beam matching,
  `scripts/fwhm_scaling_test.py`) and the β-blob beam-averaging test (`beam_averaging_blobs/`) rerun on the lowest- and
  highest-frequency maps that reach the TRIS ring, importing the original code unchanged. Adds masked-map handling, a
  "survey ⊕ own beam" width variant and a CMB-removed variant. Write-up and numbers in its `README.md`.
- **Findings:** the Haslam ranking holds at both ends (TRIS-cut kernel and major axis best); LWA1 35's own 4.65° beam
  must be included in the width (TRIS-stripe spread 0.48% → 0.20%). The blob shape change is ≈ 0.6–0.7% per unit
  ln(ν/ν_ref) on all three maps: 1.9–2.2% from LWA1 35 to TRIS, 0.3–0.5% from Stockert.

#### Beam averaging with spectral-index blobs (2026-10-19/20)

- 🧪 **`beam_averaging_blobs/`, new.** Tests whether beam averaging changes the *shape* of
  T vs LST across frequencies when regions inside the TRIS beam have different spectral
  indices.
  - Haslam 408 MHz is extrapolated with a β map: β₀ = −2.7 plus four 4° blobs on the TRIS
    stripe, with β from −2.5 to −2.9.
  - The maps are observed along the ring at 408/600/820 MHz with the TOD route's own code
    (`sed_prototype.Inputs._native_tod`), and compared with a uniform-β control.
  - Rerun: `python beam_averaging_blobs/beam_averaging_blobs.py`. Details, parameters and
    caveats are in `beam_averaging_blobs/README.md`.
- **Findings** (`HASLAM_DESTRIPED_ONLY.fits`):
  - The control is shape-invariant to 7×10⁻¹⁶.
  - With blobs, the shape changes by up to 0.27% / 0.49% (600 / 820 MHz), and the
    normalisation by less than 0.1%.
  - Without the beam the shape changes by 8% / 15%, so the beam dilutes the effect strongly.
  - The largest effect is from the blob next to Cygnus A, which the destriped map keeps:
    0.49% here against 0.23% on `HASLAM.fits`.

#### Two-step vs one-step beam matching with four FWHM scalings (2026-10-19)

- 🧪 **`scripts/fwhm_scaling_test.py`, new.** Reruns notebook 03's two-step vs one-step
  test (`PIPELINE.md` App. E) with the complicated survey beam's FWHM given to the pipeline
  as the major axis (5.010°, as in the notebook), the minor axis (3.007°), the arithmetic
  mean (4.009°) or the geometric mean (3.881°). All three kernels.
  - Rerun: `python scripts/fwhm_scaling_test.py` (~7.5 min). `--replot-maps` redraws the
    residual maps from the saved `.npz`.
  - Outputs are in `outputs/fwhm_scaling_test/`: a CSV with one row per test, a comparison
    plot, and a residual map [K] for every test.
- **Findings:**
  - The major-axis rows reproduce the notebook's table exactly.
  - The major axis gives the smallest discrepancy for every kernel and the minor axis the
    largest; the two means fall in between.
  - Every scaling still under-corrects (a positive median), because the sidelobes and error
    beam smooth beyond the main lobe.
- 📝 README "Standalone tests", `PIPELINE.md` App. E and App. F, `PLOT_INDEX.md`.

### Changed

#### Survey error bars: published values replace the 10% placeholder (2026-10-21)

- 🔧 **`sed_prototype.sigma_for` now uses each survey's own published uncertainty.**
  The new table `sed_prototype.PUBLISHED_ERRORS` holds, per file, the temperature-scale
  error (fraction of T), the zero-level error (K) and the rms noise (K), each with its
  source. They are added in quadrature for every point: σ = scale·|T| ⊕ zero ⊕ noise.

  | Survey | scale | zero level | noise | source |
  |---|---|---|---|---|
  | LWA1 35, 80 MHz | 20% | – | – | Dowell et al. 2017 (LAMBDA: 15–25%; middle used) |
  | Maipu + MU 45 MHz | 1.7% | 39 K | – | Monsalve et al. 2021 (2σ: ±0.034, ±78 K) |
  | CHIPASS 1400 MHz | 5% | 0.5 K | 30 mK | Calabretta et al. 2014; scale and zero from Stockert |
  | Stockert + Villa-Elisa 1420 MHz | 5% | 0.5 K | 17 mK | Testori et al. 2001 / LAMBDA |
  | Rhodes/HartRAO 2326 MHz | 5% | 80 mK | 30 mK | Jonas et al. 1998 (upper limits) |
  | Haslam 408 MHz (target; notebook 03 only) | 10% | 3 K | – | Haslam et al. 1982 |

  - TRIS (ring errors, or stage 1 posterior ⊕ ẑ) and ARCADE 2 (Fixsen et al. 2011) are
    unchanged. `ASSUMED_FRACTION` (10%) stays only as a fallback for a map with no entry.
  - Judgement calls, recorded in each entry's `choice`: LWA's 20% is the middle of a range;
    Maipu's are the errors of corrections (×1.076, −160 K) that are **not** applied to the
    map in `res/`; CHIPASS quotes no scale or zero error of its own; Rhodes's are limits.
  - The errors are still diagonal: a scale or zero-level error is really common to every
    point of a survey, but each point is treated as independent.
- 📝 Notebook 07 (its own copy of `sigma_for`) and notebook 03's noise test now use the
  same table; their saved outputs were made with the 10% and were not rerun.
  `model_comparison/error_bars.csv` has the same values (new kind `published`, with
  `scale_frac`, `zero_level_k`, `noise_k` columns), so no point is a placeholder any more.
- Effect on the SED fit: `PIPELINE.md` §6.

#### ARCADE 2 excluded from the beam-convolved SED studies (2026-10-15)

- ✂️ **New config key `sed.exclude_maps: [ARCADE2_315.fits, ARCADE2_341.fits]`**, applied
  by `sed_prototype.Inputs` and by notebook 07's own map list. It covers the TOD route, the
  map route, and every figure built on them (notebooks 07, 08 and 10).
  - Why: ARCADE 2 covers ~7% of the sky as a thin ring and reaches only 16–17 of the 120
    TRIS ring samples, and its 6–7.5 mK errors let those few points pull the fits.
  - Not affected: stage 2 (still beam-matches it), notebook 06 (the TRIS + ARCADE 2
    four-point cross-check) and notebook 03's method diagnostics.
- 🔁 Notebooks 07, 08 and 10 re-run. The SED-fit figures no longer reserve panels for
  ARCADE 2; notebook 07's stage 2 check and "TRIS + ARCADE 2 only" group handle its
  absence.
- **Results:**

  | | with ARCADE 2 | without |
  |---|---|---|
  | TOD route, standard Haslam (07): prediction − Haslam | +4.2 K (+17.8%) | **+4.0 K (+17.6%)** |
  | TOD route, destriped Haslam (08, 10) | +3.8 K (+12.6%) | **+3.7 K (+12.5%)** |
  | TOD route: median β; χ²/dof | −2.55; up to ~50 at ARCADE 2 | **−2.53; median 3.9** |
  | map route (07): prediction − Haslam | +0.3 K (+0.9%) | +0.4 K (+1.5%) |
  | Fig 7d: TRIS 600 / 820 ÷ non-TRIS fit | 1.24 / 1.13 | 1.23 / 1.12 |

  Figure 7d's steps at LST 255–280° and 320–345° are gone.
- 📝 Readings in notebooks 07 and 08, notes in 10, the paper captions in
  `moonrunes_paper/figures/README.md`, and `PIPELINE.md` (§1 table, §6) updated.

#### Back to the circular beam-matching kernel (2026-10-14)

- ↩️ **`beam_matching.kernel: circular_gaussian`.** The config goes back to the 23.366°
  circular Gaussian, the kernel the current stage 2 product was made with. The elliptical
  kernels (`tris_cuts`, `elliptical_gaussian`; 2026-10-13) stay implemented and tested,
  but are not in use.
- 🖼️ The paper figures in `moonrunes_paper/figures/` were regenerated
  (`notebooks/10_paper_plots.ipynb`), so the beam-matching example and the map-route SED
  fits both use the circular kernel. Their captions were updated.
- 📝 README, `PIPELINE.md` (stage map, §3, §5, §5.4, config table, A.8, A.9, C.3),
  `tris_haslam_pipeline.md`, `TODO.md` N-8, the stage 2 docstring and the notes in
  notebooks 03, 05, 07 and 08 now say which kernel is in use.

#### Stage 1's convergence cross-check switched off (2026-10-14)

- 🔕 **`tris.solver.cross_check: false`.** The exact reference solve (`woodbury_reference`)
  is no longer computed. With the flag off, stage 1 now skips it entirely instead of only
  skipping the error. That robustness is not needed yet.
  - Products record `cross_check: false` and `null` for the comparison.
  - The progress line says "exact cross-check off".
  - The manifest's solver note says the maps are limTOD's answer, unverified.
  - The function and its unit tests stay. Setting the key to `true` restores the check
    (`TODO.md` S1-8).
- 🧪 `tests/test_stage_io.py`: the end-to-end test expects `null` when the check is off.
- 📝 Notebook 03 §1 and its verdict, README, `PIPELINE.md` (§4.4 and the config table),
  `TRIS_MAP_PROVENANCE.md`, the config comments and `TODO.md` (S1-8, S1-4) all say the
  check is off.

#### Two-step against one-step beam matching (2026-10-14)

- 🧪 **`notebooks/03_diagnostics.ipynb`, new last section.** Using the pipeline's own
  `match_to_tris_beam`, Haslam goes either through a survey beam first and then the
  pipeline, or straight to the target. The survey beam is a Gaussian control or a
  complicated beam (elliptical main lobe, sidelobe ring, error beam); all three kernels.
  Plots: the beam, a 2 × 3 grid of difference maps, the difference along the TRIS ring,
  and a summary.
- **Findings:**
  - With a Gaussian survey beam the two paths agree to 0.02% (Gaussian kernels) and 0.14%
    (the default measured-beam kernel).
  - With the complicated beam they differ by 0.3–0.5% rms, up to 1–2% on the Galactic
    plane: the pipeline only knows the main lobe's FWHM.
- 📝 `PIPELINE.md` Appendix E.

#### Convolve → deconvolve with limTOD (2026-10-14)

- 🔁 **`notebooks/03_diagnostics.ipynb`, new last section.** Every SED map is convolved with
  limTOD's circular (23.366°) and elliptical (19.155° × 23.366°) Gaussian beams, pointed at
  every nside 16 pixel with TRIS's geometry. It is then deconvolved with limTOD's
  `wiener_filter_map` at three noise levels. Plots: the pointed beams and the operators'
  singular values; input / convolved / deconvolved / residual maps; elliptical − circular
  for both; transfer function, correlation by ℓ and error by smoothing scale; every map.
- **Findings:**
  - Convolution matches `generate_TOD_sky` exactly.
  - Deconvolution does not return maps pixel by pixel: 14–39% rms, worse than the convolved
    map at realistic noise.
  - It returns them on large scales: correlation ≥ 0.9 to ℓ ≈ 21 at 10⁻⁴ noise, and
    2.8–3.8% at 10° on Haslam.
  - The elliptical beam passes more modes (377 against 320 above 10⁻²) and recovers better
    for every map.
  - Elliptical and circular convolved maps differ by 1.9% rms.
  - Implication for stage 1's TRIS maps and the map route: `TODO.md` N-10.
- 📈 **The same round trip along the TRIS ring (LST)**, with each survey's pipeline noise
  (`sed_prototype.sigma_for`):
  - Convolution changes the maps mainly at the Cygnus crossing: −45 to −80% at LST ≈ 305°,
    +15 to +55% beside it, within ±20% elsewhere.
  - The circular and elliptical beams are nearly identical along the ring, which runs along
    the H plane.
  - Deconvolution scatters by 10–25% rms at every LST.
- 📝 `PIPELINE.md` Appendix D.

#### Harmonic-space checks (2026-10-14)

- 🔭 **`notebooks/03_diagnostics.ipynb`, new section "Harmonic-space checks"**, in response
  to a collaborator's concern about rotating and convolving in harmonic space. It finds the
  four harmonic operations in the code (Haslam's `rotate_map_alms`, `hp.smoothing`, and
  limTOD's beam and horizon pointing via `map2alm`/`rotate_alm`). It tests each with its
  C_ℓ and against an exact pixel-space calculation:
  - **H1:** the C_ℓ of every map, raw and smoothed;
  - **H2:** the TRIS beam's C_ℓ at nside 16–128;
  - **H3:** limTOD's pointed beam against the exact beam;
  - **H4:** the TODs;
  - **H5:** Haslam's rotation ringing;
  - **H6:** harmonic against real-space smoothing.
- **Findings:**
  - Every effect on the nside 16 products is ≤ 0.05% (harmonic smoothing on Haslam ≤ 0.5%,
    but not ringing).
  - The smoothed maps' a_ℓm reach 10⁻¹¹–10⁻¹⁶ of the peak by ℓ ≈ 35.
  - Haslam's full-resolution rotation rings at bright sources: 33 K at Cas A, which reaches
    notebook 08's Figure 11 (`TODO.md` N-11).
  - In the TODs, nside 16 sky pixelisation (≤ 0.53%) outweighs every harmonic step
    (≤ 0.043%).
  - limTOD's beam orientation confirms the elliptical kernel's +7° E plane.
- 📝 `PIPELINE.md` Appendix C; `TODO.md` N-11, N-12, S1-7 and N-9 updated.

#### The map route smooths to the measured, elliptical TRIS beam (2026-10-13)

> **Run status (2026-10-14):**
> - **Done:** `tests/test_elliptical_beam.py` passes (15/15), and notebook 03's elliptical
>   section ran. The default `scaled_cuts` kernel for broad-beam surveys gives ARCADE 2 an
>   effective beam 8% / 5% too wide (`TODO.md` N-9).
> - **Not re-run yet:** stage 2 and notebooks 06–08 (`TODO.md` N-8).

- 🎯 **`beam_matching.kernel`** chooses the shape of the TRIS beam that stage 2 and the map
  route smooth every survey to. The default changes from the old 23.366° circular Gaussian
  to **`tris_cuts`**:
  - it is the measured beam: the archive E/H cuts, blended with cos²φ / sin²φ, the same
    beam the TOD route and stage 1 use;
  - it is 19.155° (E, north–south) × 23.366° (H, east–west), with the E plane 7° east of
    north, as TRIS saw the sky from the zenith;
  - it is zero beyond 90° (the horizon).

  `elliptical_gaussian` (same widths and orientation, Gaussian shape) and
  `circular_gaussian` (the original) remain available.
- 🧮 **New module `src/moonrunes/elliptical_beam.py`.**
  - `hp.smoothing` cannot apply a kernel that differs by direction, so the elliptical
    kernels are a **direct sum on the sphere**: normalised convolution over a working grid
    (`nside_work` 64), evaluated at 2×2 points per nside 16 pixel.
  - The survey's own beam is allowed for per axis:
    - beams that widen the E core by ≤ 0.2% are neglected (Haslam, Stockert, Rhodes,
      CHIPASS);
    - broader ones get the measured beam compressed to the quadrature widths (LWA1 35/80,
      Maipu, ARCADE 2: `broad_native_kernel: scaled_cuts`), or the elliptical Gaussian.
  - Derivation in the module docstring and `PIPELINE.md` Appendix A.9.
- 🔌 **`stage2_beam_matching.match_to_tris_beam()`**: smooth + mask + regrid for one map,
  with whichever kernel is configured.
  - `beam_match()`, `sed_prototype.Inputs` (the map route) and `notebooks/07` all call it,
    so they cannot drift apart.
  - The stage 2 npz gains a `kernel` field.
  - The manifest records, per map, the kernel's widths, its compression, the expected
    effective beam and why that kernel was chosen.
- ⚙️ **Config:** `beam_matching.kernel` plus `beam_matching.elliptical.*`, each with its
  reasoning. `target_fwhm_deg` now applies only to `circular_gaussian`.
- 🧪 **Tests:**
  - new `tests/test_elliptical_beam.py`: geometry, exact agreement with limTOD's cut
    beam, half-power widths, compressed widths, constant-sky and mask-edge behaviour, a
    point source 19.155° north–south and 23.366° east–west, and the keep/floor rules;
  - `tests/test_stage_io.py` updated for the kernel choice.
- 📊 **`notebooks/03_diagnostics.ipynb`, new last section, "The elliptical beam for the map
  route"**, 7 plots:
  1. the three kernels on the sky;
  2. profiles against the archive cuts and the weight in the shoulders;
  3. the kernel per survey;
  4. a point source through the real code, for narrow and broad native beams;
  5. Haslam, old vs new: maps, difference, along the ring;
  6. ARCADE 2's mask edge;
  7. the change in every map-route survey.
- 📝 Notes added to notebooks 05 (circular kernel only), 07 (section 2 now uses
  `match_to_tris_beam`) and 08 (which kernel the figures use). Updated: README,
  `PIPELINE.md` (§3, §5, §9, §10, §11, A.8, new A.9), `TODO.md` (N-8 to N-10) and
  `tris_haslam_pipeline.md` (Step 0.1 marked superseded).

#### Beam-width tests and SED diagnostics (2026-10-12/13)

- 🧪 **`notebooks/09_beam_fwhm_tests.ipynb` and `tests/test_beam_fwhm.py`.** These test the
  beam widths: the σ ↔ FWHM constant; the FWHM function on Gaussian and archive cuts; deg/rad
  and σ/FWHM slips; the cut beam's principal planes (φ = 0 → E, 90° → H); smoothing width
  and quadrature; and a round trip (Gaussian beam map → cuts → FWHM) at nside 16–1024. It
  passes at 0.1% from nside 128. Findings:
  - linear-in-dB interpolation on the archive's 2° grid reads widths ~0.08° narrow;
  - `hp.smoothing` at 23.366° has a real-space FWHM of 23.425°;
  - a beam map at nside 16 reproduces the Gaussian only to ~1%.
- 📊 `notebooks/07`: β vs LST with surveys left out (fits with ≥ 3 points).
- 📊 `notebooks/08`: Figure 13, β from every T–T pair, with bootstrap error bars.
- 📊 `notebooks/01`: a table of the nside of every file the notebook loads.

#### Notebook 06 on nside 16, TRIS + ARCADE 2 only (2026-10-07)

- 📉 **`notebooks/06_sed_fit.ipynb` ported to nside 16.** It reads the stage 1 and stage 2
  products instead of rebuilding maps on the old grid; nside, beams and floor come from
  the config (ARCADE 2 11.6°). It stays the four-trusted-point analysis. New: full-sky maps
  of best-fit Haslam, old Haslam, their difference, β, χ²/dof and σ(β).
- **Fitted at every TRIS-stripe pixel:** 4 points where both ARCADE 2 bands exist (153), 3
  where one does (14), and the 2 TRIS bands elsewhere (1105; exact, no χ²). With ARCADE 2:
  median β −2.79, best fit − Haslam **+3.4 K (+9.6%)**. TRIS only: σ(β) 0.5–0.9, +31%,
  driven by ẑ (+1% without it).

#### Paper plots (2026-10-07)

- 🖼️ **`notebooks/08_paper_plots.ipynb`** and **`moonrunes.sed_prototype`** (notebook 07's
  calculation as a module). Target: `HASLAM_DESTRIPED_ONLY.fits`. Figures: maps used, target
  with the TRIS stripe, stripe slices, TRIS maps and TODs, pixel-by-pixel (5) and
  sample-by-sample (5b) SED fits with residuals, fit-residual maps, old vs best-fit Haslam
  for both routes, residuals by survey, χ² and β maps, full-resolution corrected Haslam for
  both routes, and T–T corner plots (raw and beam-convolved). PNG + PDF in
  `outputs/paper_plots/`.
- With the destriped-only target, the TOD route finds Haslam **+12.6%** low along the ring
  (+3.8 K). The map route finds −0.4%, still circular.

#### Stage 3 prototype with every map: TOD route and map route (2026-10-07)

- 📈 **`notebooks/07_sed_calibration.ipynb`**. Every usable `res/` map plus TRIS, fitted with a
  power law two ways. The **TOD route** uses limTOD TODs of the native maps, fits each ring
  sample against the real TRIS data, and maps the correction TOD with stage 1's Wiener
  filter and a zero-mean prior. The **map route** fits each TRIS-stripe pixel on the
  nside 16 grid, with stage 1 + ẑ. Plots: old Haslam (nside 16), the corrected maps for
  both routes, and a hybrid (corrected in the stripe, old elsewhere).
- 🔎 **The routes disagree.** The TOD route finds Haslam ~18% low along the ring; the map
  route finds +0.9%. The map route is circular through stage 1's GSM2008 prior.
- 📦 **`moonrunes.sky_maps`**: the `res/` reader (four layouts, metadata from keywords,
  comments, filename or `EXTERNAL_METADATA`), now with a numeric native beam for the
  smoothing kernel. The WMAP K / Ka files are excluded: their headers say V band, and their
  monopole is removed.

#### Stage 2 — complete: matched to TRIS, written to disk (2026-10-07)

- 🛰️ **Stage 2 now rotates, smooths, masks, regrids and writes.** Per map: the declared
  frame is checked, the map is rotated into equatorial (`moonrunes.frames`: harmonic for
  Haslam, mask-aware pixel rotation for ARCADE 2), smoothed to the TRIS beam (23.366°),
  cut to observed pixels whose smoothed weight reached 0.5, and regridded onto
  `beam.nside_new` = `tris.nside` = 16 with the new `regrid_to_nside` (a coarse pixel
  needs every sub-pixel observed). TRIS is not touched: everything else is matched to it.
- 💾 **`run_stage2()` / `load_stage2()`** write and read
  `outputs/stage2/beam_matched_<run>.npz` + `stage2_manifest.json`, with stage 1's
  overwrite guard. `run_pipeline.py --stage 2` calls it. ARCADE 2 keeps 170 / 192 pixels
  at nside 16, after removing 4 / 5 extrapolated ones.
- 📏 **`beam_matching.arcade2_native_fwhm_deg` 12.0 → 11.6°**, as both headers state. New
  config key `beam_matching.coarse_pixel_rule: all_children`.
- 🧪 Six new tests: the regrid rule, no resolution invented, the ARCADE 2 beam against its
  headers, the TRIS grid and no extrapolation end to end, and the product round trip.
- 🗑️ The stale `outputs/stage2/` product from the deleted `stage2_haslam_prep` was removed.


#### Everything on nside 16 (2026-10-02)

- 🔢 **`tris.nside`, `tris.nside_hires` and `beam.nside_new` all 16.** nside 16 is the grid
  the TRIS stage 1 maps are made on, so every other product is put on the TRIS grid. Stage 1
  at the current config: 1272 pixels, reduced χ² **16.1 / 5.9**, zero level +1.737 / +0.516 K, and
  23 / 7 pixels tightened by >5%. With `nside_hires: 64` the same run gives 3.09 / 3.22 and
  388 / 201; `nside_hires` is kept at 16 by decision (`TODO.md` S1-7).
- 🔓 **The `nside_new` ↔ `bayesian_func.nside_for_beam` check is no longer fatal.** It
  only guarded the superseded Gibbs path (the heuristic derives a coarser grid from the
  21.26° beam). A mismatch is now
  recorded in the manifest (`matches_bayesian_skymap`) and warned about. The `beam_deg`
  = E/H mean check stays fatal.
- 🧪 `test_stage1_refuses_to_clobber` pins `run.overwrite: false` itself, instead of
  depending on the config default.
- 📝 ARCADE 2's beam: both headers carry `BEAMSZ = 11.6°`. The config, stage 2's dataset
  table and the README said there was no beam keyword. The 12.0° in use is unchanged
  (`TODO.md` N-4).


#### Stage 1 — solver, resolution, and pixel selection

- 🔁 **The solve is now limTOD's own** — `TRISMapMakingInputs.solve`, wrapping
  `limTOD.wiener_filter_map` (dense, direct). The in-repo matrix-free `krylov_map_solve`
  stays as the test suite's reference implementation but no longer makes the maps. The
  exact-reference cross-check is retained and repointed at limTOD's answer: a dense factorisation
  with `εI` on the diagonal can return a biased result without raising. New key
  `tris.solver.regularization`; the Krylov keys and stage 5's `require_solver_info_zero`
  are marked unused.
- 📉 **`tris.nside` 64 → 16.** This is supposed to be a coarse map: the beam is 19–23°, so
  nside 16 is already ~6× finer than the resolution behind it, and `beam.nside_new` is 16
  as well, so every downstream product sits on the same grid as the TRIS maps. It also puts the dense solve back in reach — at nside 64
  the 25,705-parameter normal matrix needs ~25 GB and did not complete on 16 GB of RAM.
- 🎯 **Pixel selection is now on beam response, not a declination band.** A pixel is
  retained if the beam response exceeds `tris.beam_response_threshold` (0.01) of the peak
  for at least one observation — the collaborator's criterion, and the physical statement,
  where the old `dec_half_width_deg` cut was only a geometric proxy. 1272 pixels, 41.4% of
  sky, 99.60% of beam power retained.
- 📐 **Threshold convergence checked at 0.003 / 0.01 / 0.03.** The maps agree to within
  **0.40 σ** (max, on the common pixel set) — immaterial, as expected. The **fitted zero
  level does not converge**: it moves +0.112 K at 600 MHz and +0.065 K at 820 MHz across the
  decade, monotonically, which is 2.6 σ / 2.4 σ of its own fitted error though only 22% /
  13% of the 0.5 K prior width. Recorded as a ~0.1 K pixel-selection systematic on the zero
  level; the maps are unaffected.
- 📊 Products regenerated. Reduced χ² falls from 33.6 / 13.3 to **3.09 / 3.22**, and the
  posterior σ now beats the prior by >5% for **388 / 201** of 1272 pixels rather than 2 of
  25,704 — coarsening the grid is what bought that.

### Added

#### Coordinate frames — everything into Equatorial

- 🧭 **`src/moonrunes/frames.py` + `notebooks/04_coordinate_frames.ipynb`.** Frames read
  from the files, then checked against the sky: Haslam declares `COORDSYS='GALACTIC'`;
  ARCADE 2 declares `SKYCOORD='Galactic'` (not `COORDSYS`, which is why it looked
  unstated); the stage 1 TRIS products declare nothing and are equatorial by construction.
  Haslam is rotated G→C in harmonic space (full sky, band-limited at nside 512); ARCADE 2
  in pixel space with the mask interpolated separately, so the exact-zero fill cannot leak
  into edge pixels (a plain `rotate_map_pixel` contaminates ~157 per band). Sgr A*, Cas A
  and Tau A land within 0.08–0.15° of their literature positions after rotation.
  Beam matching does not use it yet.

#### Notebooks and analysis

- 📈 **`notebooks/03`**: residuals with error bars and the RMS maps, with the derivation:
  r = d − A x̂, Cov(r) = N − AΣAᵀ, and a χ² normalised by n − Σhᵢ.
- 🛰️ **`notebooks/05_beam_matching.ipynb`**: Haslam and ARCADE 2 smoothed in the equatorial
  frame. UNSEEN counts, the final ARCADE 2 mask (observed **and** w ≥ 0.5, removing 4 / 5
  extrapolated pixels), and the four-way overlap with TRIS (33 pixels on the earlier, coarser grid with the
  "all 4 children" rule; to be re-run at nside 16).
- 📉 **`notebooks/06_sed_fit.ipynb`**: stage 3 prototype. The 5-point SED and the 4-point
  power-law fit at the 33 overlap pixels, ARCADE 2 in thermodynamic → RJ units with
  calibration errors from Fixsen et al. 2011, and the sensitivity to TRIS beam matching and
  the zero level (`TODO.md` N-2, N-3).
- 🗺️ **`notebooks/01` sections 05–06**: every map in `res/`, discovered automatically, with
  the TRIS strip, Dec +30–55° strip slices and a frequency / frame / nside / beam table.
  Then the TRIS TOD each map predicts, from limTOD's `generate_TOD_sky` (`TODO.md` N-5).
- 📦 **`res/`**: Maipu + MU 45 MHz and Stockert + Villa-Elisa 1420 MHz added as reference
  maps (`DATA_SOURCES.md`).

#### Documentation

- 🔬 **`TRIS_MAP_PROVENANCE.md`** — a collaborator-facing record of exactly what produced
  `outputs/stage1/tris_maps_tris_haslam_v0.npz`: the limTOD commit (`dbe720b`, with the dirty
  files enumerated and shown to be outside the map-making path), the five `limTOD.tris` entry
  points in call order, and every setting tagged as archive-stated, our decision, or a limTOD
  default we did not pass. Records the convergence and per-frequency results, and three gaps
  found while writing it: the manifest pins limTOD's path but not its commit,
  `beam.selfrot_deg` is never read, and the E/H HPBWs are never checked against the cut table
  limTOD measures them from.
- 📓 **`notebooks/01_explore_tris_data.ipynb` → `01_explore_data.ipynb`**, now three
  sections rather than one: the TRIS products stage 1 consumes (unchanged), a first look
  at ARCADE 2 at 3.15/3.41 GHz, and the Haslam 408 MHz download. The ARCADE 2 section
  establishes the three properties Step 0 of the SED design turns on — `NESTED` at
  nside 16, no `COORDSYS` in the header, and the unobserved ~93% of the sky stored as
  exact zero instead of `hp.UNSEEN`, which is what the mask-before-smoothing procedure
  exists to handle. Re-executed end to end; every stage 1 number in it is byte-identical
  to the previous run.
- 📡 **`DATA_SOURCES.md`** — the download list. Every external dataset, its LAMBDA archive
  page, the config key that points at it, and the stage that reads it: the TRIS product
  table (stage 1), the reprocessed Haslam 408 MHz map (stage 2), and ARCADE 2, whose two
  lowest bands are trusted points 3 and 4 of the SED design. Also states
  where the two non-LAMBDA inputs come from — Berkhuijsen (1972), still missing, and
  GSM2008, which `pygdsm` generates rather than downloads.

#### Stage 2 — Haslam prep

- 🌌 **Stage 2 — Haslam, the gain operator, and the prior covariance**
  (`moonrunes.stage2_haslam_prep`). Loads the destriped/desourced reprocessed Haslam
  408 MHz map (Remazeilles et al. 2015) from pygdsm's astropy cache, degrades it to
  `bayesian_skymap`'s fixed nside 128 grid, and removes the CMB monopole — necessary, not
  cosmetic, because that repository models `sky_gain` as a pure power law with no monopole
  term, so 2.7 K left in it is a forward-model error everywhere. Runs in about 20 s.
- 🧭 **Galactic as the common frame** (`haslam.frame`). Haslam is native there and is never
  rotated; the gain and β regions become galactic latitude bands, which is what
  `bayesian_skymap`'s own `region_operator` assumes and what is physical for synchrotron β.
  Stage 3 will rotate stage 1's equatorial TRIS band in, and the ragged result is carried
  by `pixel_mask`.
- 🧩 **Region operators pinned against `bayesian_skymap`'s own** — `region_operator` builds
  the `(npix, nregions)` 0/1 indicator matrices, and the test suite asserts byte equality
  with `tests/make_test_data.region_operator` rather than trusting a lookalike. An
  `equal_area` geometry is available as an alternative: equal-spaced colatitude bands leave
  a 7.6× pixel-count imbalance between the polar and central gain regions.
- 📐 **Three prior sources, switched by `haslam.prior.source`** — `gsm2008` (default),
  `tris_stage1` (stage 1's TRIS maps extrapolated 600 → 408 MHz; genuinely
  Haslam-independent across the 52% of sky the ring covered, deliberately loose outside),
  and `berkhuijsen_1972` (the real answer, which raises with a clear message until
  `paths.berkhuijsen_map` points at a file). Swapping Berkhuijsen in when it arrives is one
  config line and nothing else in the pipeline changes.
- 🚨 **`prior_is_haslam_derived`, and the measurement behind it.** The default prior is
  *not* independent of Haslam, and `notebooks/02` quantifies it rather than asserting it:
  GSM2008 is fitted to eleven input surveys, one of which is Haslam 408 MHz, and at 408 MHz
  its reconstruction tracks Haslam at **log-space r = 0.993, median ratio 1.005, 80% of
  pixels within 5%**. Stage 2 sets the flag in its manifest and prints a warning, so a gain
  recovered against this prior cannot be mistaken for a result.
- 🔧 **Substitutes for two more truth-informed quantities**, both recorded rather than
  silently supplied: `covG`, which `gibbs_*.py` builds from the *true* gain as
  `(0.15 × (1 + g_true))²` (stage 2 emits the same width at `g = 0`), and `sky_haslam`,
  which those scripts read as the true sky (stage 2 emits `s0_init_k` from
  `haslam.s0_init`).
- 📓 **`notebooks/02_build_berkhijsen_prior.ipynb`** — committed with outputs. The Haslam
  map, the region operators with the polar-imbalance comparison, section 3's circularity
  measurement, the `tris_stage1` prior against Haslam (median ratio 1.027, 16–84%
  [0.982, 1.084] in band — a real absolute cross-check), and a table of which `.npz` keys
  stage 2 has supplied and which stage 3 still owes.
- 📦 **`load_stage2()` / `Stage2Products`**, matching stage 1's loader contract, and
  `run_pipeline.py --stage 2`.

#### Stage 1 — TRIS maps

- 🗺️ **Stage 1 — calibrated TRIS maps** (`moonrunes.stage1_tris_maps`). Reads the public
  TRIS rings at 600 and 820 MHz, builds the operator, noise and prior through
  `limTOD.tris.build_tris_mapmaking_inputs`, and solves the prior-regularized map with
  the zero level carried as an explicit nuisance parameter. The forward model is the one
  verified in limTOD's own TRIS walkthrough — cut-built beam, parked-zenith geometry with
  `selfrot = −7`, horizon mask, Rayleigh–Jeans plus the CMB monopole — and nothing here
  re-derives it. The whole stage runs in about 12 s for both frequencies.
- 🧮 **The direct Krylov solve, per Blocker 3** — `krylov_map_solve` wraps the MAP
  normal-equation matvec in a `scipy.sparse.linalg.LinearOperator` with a Jacobi
  preconditioner and solves it with GMRES, mirroring `bayesian_skymap`'s amplitude step.
  No dense normal matrix is ever formed (it would be 25705² at the working nside), and
  the low-rank fast-path acceleration that is wrong by eight orders of magnitude at two channels is
  never on that path.
- 🔍 **A convergence check that is not "did it finish"** — an exact reference solve computes
  the MAP mean and the full posterior variance exactly via
  `S = P − P Aᵀ (N + A P Aᵀ)⁻¹ A P`, whose inner matrix is only `n_samples × n_samples`
  (120 for a TRIS ring). Stage 1 compares its Krylov answer against it **in posterior σ**
  and refuses to write a product that disagrees. The residual norm alone is a poor proxy
  here: the prior width spans an order of magnitude across the band, so a 1e-9 residual
  leaves the map ~1e-3 σ off. Measured agreement is 1.6e-07 σ at 600 MHz and 1.8e-09 σ at
  820 MHz. This identity is also what supplies the posterior σ, which a Krylov solve
  alone cannot give — and it is exact at a single frequency, not the multi-frequency
  approximation Blocker 3 rules out.
- 🧭 **Blocker 1 is re-derived at run time, not trusted.** Stage 1 calls
  `bayesian_func.nside_for_beam` on the configured E/H mean beam width (21.2605°) and
  refused to start if the pinned `beam.nside_new` had drifted from it, so the grid stage 3
  hands to `bayesian_skymap` could not move under the pipeline unnoticed. Relaxed to a
  warning on 2026-10-02, when `nside_new` moved to 16 to match the TRIS maps.
- 📖 **Config loader** (`moonrunes.config`) with dotted lookup. `require()` treats a
  `null` as *"still an open decision"* and raises, rather than letting a stage invent a
  default; `section()` supplies the blocks each stage records in its manifest.
- 📦 **Stage products and their loader** — `outputs/stage1/tris_maps_<run>.npz` (maps,
  posterior σ, prior, reference template, zero levels, residuals, beam coverage, stacked
  over frequency) and `stage1_manifest.json` (the config block used, per-frequency
  diagnostics, and each blocker's resolution). `load_stage1()` / `Stage1Products` is the
  single contract the notebooks and stages 3 and 5 read them through.
- 🚀 **`run_pipeline.py`** — the single `--stage N` entry point. Stages 2–5 say they are
  unimplemented rather than half-running.
- 📓 **`notebooks/01_explore_tris_data.ipynb`** — the inputs, committed with outputs. The
  four archive products, the zero-level and uncertainty-floor decisions traced to the data
  that forces them, the Rayleigh-Jeans/CMB-monopole discriminator, the cut beam against
  the Gaussian, and the GSM2008 template forward-modelled against the absolute rings:
  **1.92 K cold at 600 MHz, 0.56 K cold at 820 MHz**.
- 📊 **`notebooks/03_diagnostics.ipynb`** — the stage 1 results, committed with outputs.
  Convergence, reduced χ² and residuals, the fitted zero level against the template
  deficit, prior→posterior σ shrinkage and z-scores, the band maps, the spectral index
  between the two solved maps (median β = −2.701, 92% inside the −3.0…−2.5 range stage 5
  expects), and a preview of the degrade to `nside_new` on the earlier, coarser grid (432 of 768
  pixels touched, 364 fully covered). At nside 16 the TRIS maps are already on the
  `nside_new` grid.
- ✅ **Tests** (`tests/test_stage_io.py`, 14 of them). The algebra is pinned against dense
  references — the Krylov solve for all three methods, the exact reference posterior including
  the regime where pixels outnumber samples, and `lhs_diagonal` against
  `bayesian_func.estimate_diag_precond` — and runs anywhere; the end-to-end tests skip
  themselves without the archive or `pygdsm`.

### Changed

- 🔢 **`haslam.n_gain_regions` is 12, not 6.** The config's comment claimed 6 "matches the
  six regions bayesian_skymap initialises" — but those six are the **spectral index**
  regions (`BETA_REGIONS`); that repository's own generator uses `ngain = 12`. The gain is
  meant to have more freedom than β. `n_beta_regions: 6` is now explicit and is validated
  against the length of `assemble.beta_init`.
- 📍 **`paths.haslam_map: null` no longer means "open decision".** It means "use the
  reprocessed Remazeilles map that pygdsm already caches", which is the map the pipeline
  document asks for. Set it to override with a local file.
- ⚙️ `configs/run_config.yaml` gained `haslam.frame`, `haslam.remove_cmb_monopole`,
  `haslam.region_geometry`, `haslam.n_beta_regions`, and the extra
  `haslam.prior` keys (`out_of_band_frac_sigma`, `gain_frac_sigma`), with the prior
  source's circularity written out in full next to the value that causes it.
- ⚙️ `configs/run_config.yaml` gained the blocks stage 1 consumes: `tris.prior`
  (GSM2008 template, 10% relative and 0.3 K absolute width), `tris.solver`
  (method, tolerances, cross-check), `tris.nside_hires`, and `paths.limtod`.
  `paths.tris_archive_dir` now points at a local archive copy.
- 🎯 `tris.solver.rtol` is **1e-12**, not the 1e-9 used elsewhere in the pipeline, for the
  reason recorded above: on this operator a 1e-9 residual is not a converged map.
- 🔢 `tris.uncertainty_floor_k: null` now means *"use this ring's smallest strictly
  positive statistical error"* rather than blocking. Each published ring contains exactly
  one row with a 0.000 K error, so a floor is mandatory and the data supplies the only
  defensible value; the choice is recorded per frequency in the manifest.
- 🙈 `outputs/` is gitignored (except `.gitkeep`), matching what the config already
  claimed.

### Known limitations

- 🚧 **Blocker 2 is unresolved.** There is no true sky for real data and Berkhuijsen (1972)
  is not in hand, so GSM2008 currently stands in as *both* the prior template and the
  reference truth. That is circular by construction: the manifest records
  `prior_equals_truth: true`, and agreement between a stage 1 map and `truth_k` is not
  evidence of anything. The Berkhuijsen-based prior is stage 2 work
  (`notebooks/02_build_berkhijsen_prior.ipynb`, still empty).
- 📉 **The stage 1 map is very largely the prior, redisplayed.** 120 samples carry roughly
  15 independent numbers, so at any usable nside the pixels outnumber the data: the median
  posterior/prior σ ratio is 0.99998. The information lives in the beam-convolved
  directions and in the fitted zero level, which lands on the template's own error
  (+1.7361 ± 0.0108 K at 600 MHz against a measured 1.92 K deficit) with a monopole
  degeneracy of ≈1. Read that offset as a model-mismatch amplitude, not as a measurement
  of the archive's zero point. This is the honest limit of one ring plus a template that
  is kelvins off, and it is what stage 2's prior work exists to fix.
- ⚠️ **The shipped prior is circular, by choice.** With `haslam.prior.source: gsm2008`,
  stage 5's comparison of a recovered gain against the ~3% / ~0.91 K literature benchmark
  tests the plumbing, not the sky. This is deliberate — it lets stages 3–5 be built and
  exercised end to end — but a number that comes out of it is not a result. Use
  `tris_stage1` for a genuinely independent prior over the band TRIS observed, and
  `berkhuijsen_1972` once that survey is in hand.
- 🚧 **Stage 2 implements `op_alm = 0` only.** The spherical-harmonic gain parametrisation
  (`op_alm = 1`) builds no operator matrix — `bayesian_skymap` synthesises the gain from
  alms internally — so there is nothing for this stage to produce; it raises rather than
  writing a misleading product.
- 🧱 **Stages 3–5 are not implemented.** `src/moonrunes/stage3_assemble.py` through
  `stage5_validate.py` are empty.
