# moonrunes — notebook index

As of 2026-11-16. There are three notebooks in `notebooks/`. For the individual plots see `PLOT_INDEX.md`; for which
notebook cell makes which paper figure and table, and the open items, see `HANDOVER.md`.

| Notebook | What it is for | Run order |
|---|---|---|
| `diagnostics.ipynb` | Checks and plots of the **inputs and the stage-1 maps**: what every input file contains, whether every map is in the frame it declares, and whether the stage-1 TRIS map-maker converged | independent (needs the stage-1 products in `outputs/stage1*/`) |
| `tests_analysis/` (one notebook per Part, A–L; `README.md` there) | **The analyses and checks around the pipeline, with their results**: beam-matched TODs, spectral fits, model comparison, the 408 MHz sky models, T–T indices, the Haslam scale/offset fits, zero-level and survey-subset robustness, the prediction at other frequencies (RHINO), and the beam tests | after the pipeline; before `10_paper_plots` (Parts C–I make paper numbers) |
| `10_paper_plots.ipynb` | **The paper's figures** (`fig1_…`–`fig10_…` and `rhino_tods_tris_beam`), written to the paper repo's `figures/` folder | last |

**The pipeline comes first:** `python run_pipeline.py --step all` (`PIPELINE.md`) makes the TODs, the model comparison, the model selection, our spectral model, the 408 MHz and the RHINO predictions in `paper_outputs/pipeline/`. `tests_analysis/` and `10_paper_plots.ipynb` read those products; the Parts do not build the TODs, the power-law fits or the RHINO prediction themselves. Every notebook prints its total run time in the last cell (the first code cell records the start). `10_paper_plots` Figure 7 (cell 18) and the UltraNest and maximum-likelihood fits of Parts B, C and D run in parallel with `joblib`. `test_results.ipynb` is the merged notebook the Parts were split from (2026-11-17), kept as a record.

All three are committed **with their outputs**. The outputs in `diagnostics.ipynb` and the `tests_analysis/` Parts are the
ones stored in the notebooks they were built from; **not re-run after the merge**. Each
Part of those two notebooks begins with a banner (what it was, what it reads and writes, which paper items it feeds)
and with its own setup cell, so a Part can be run on its own.

## Where the old notebooks went

The originals of the merged notebooks are in `../moonrunes_archive_2026-11-16/notebooks/`. The other six
(02, 05, 06, 07, 08, 13) are in git at commit `a5c827c`.

| Old notebook | Now |
|---|---|
| `01_explore_data` | `diagnostics.ipynb` Part A |
| `02_build_berkhijsen_prior` | removed 2026-11-16 (superseded; it imported a deleted module) |
| `03_diagnostics`, cells 0–38 (stage-1 diagnostics, prior-only runs) | `diagnostics.ipynb` Part C |
| `03_diagnostics`, cells 39–114 (beam, harmonic and beam-matching tests) | `notebooks/tests_analysis/` Part K2 |
| `04_coordinate_frames` | `diagnostics.ipynb` Part B |
| `05_beam_matching` | removed 2026-11-16 (superseded by the `tris_cuts` kernel) |
| `06_sed_fit` | removed 2026-11-16 (TRIS + ARCADE 2 only; ARCADE 2 is excluded from every paper analysis) |
| `07_sed_calibration` | removed 2026-11-16 (prototype of `sed_prototype`) |
| `08_update_plots` | archived. **Its cells 50 onward wrote `paper_outputs/nb08_paper_plots/fig15_*.csv`, which `10_paper_plots.ipynb` cell 22 checks Figure 9's TRIS-only fits against; nothing in the repo writes those files** (see `HANDOVER.md`) |
| `09_beam_fwhm_tests` | `notebooks/tests_analysis/` Part K1 |
| `10_paper_plots` | unchanged |
| `11_model_comparison_diagnostics` | `notebooks/tests_analysis/` Part C |
| `12_sky_models_vs_fit_408` | `notebooks/tests_analysis/` Part B |
| `13_sky_model_residuals_408` | removed 2026-11-16 (duplicated Part B's numbers) |
| `14_ring_tt_betas` | `notebooks/tests_analysis/` Part E |
| `15_model_choice_decision_sheet` | `notebooks/tests_analysis/` Part D |
| `16_haslam_vs_power_law_408` | `notebooks/tests_analysis/` Part F |
| `17_haslam_gain_offset_calibration` | `notebooks/tests_analysis/` Part G |
| `18_beam_matched_tods` | `notebooks/tests_analysis/` Part A |
| `19_power_law_tod_prediction` | `notebooks/tests_analysis/` Part J |
| `20_zero_level_offsets` | `notebooks/tests_analysis/` Part H |
| `21_robustness_check` | `notebooks/tests_analysis/` Part I |
| `scripts/fwhm_scaling_test.py` (results) | shown in `notebooks/tests_analysis/` Part K4 (the script stays in `scripts/`) |
| `scripts/beam_averaging_blobs/` (results) | shown in `notebooks/tests_analysis/` Part K3 (the scripts stay there) |
| `scripts/beam_tests_lwa35_stockert1420/` (results) | shown in `notebooks/tests_analysis/` Part K5 (the scripts stay there) |

Text inside the merged cells (and in `CHANGELOG.md`, `TODO.md`, `inputs_audit/inputs_audit.md` and the frozen
`scripts/model_comparison/archive/`) still uses the old numbers; this table translates them.

---

## `diagnostics.ipynb`

| Part | Was | Reads | Writes |
|---|---|---|---|
| A — The data the pipeline consumes | `01_explore_data` | `res/`, the TRIS rings, `outputs/stage1_nside8/` | inline figures |
| B — Coordinate frames | `04_coordinate_frames` | Haslam and ARCADE 2 in `res/`, the TRIS rings, `outputs/stage1_nside8/` | a JSON summary |
| C — Stage-1 diagnostics and prior-only runs | `03_diagnostics` cells 0–38 | `outputs/stage1/`, `outputs/stage1_nside16_high_noise/`, `outputs/stage1_nside16_prior_only/` | inline figures |

### Part A — The data the pipeline consumes
- **§01 TRIS:** the four archive products, where `tris.zero_level_sigma_k = 0.5` comes from, the temperature convention,
  the beam, the GSM2008 prior template, and what 120 samples can support.
- **§02–03:** ARCADE 2 at 3.15 and 3.41 GHz, including how its unobserved pixels are marked; Haslam 408 MHz.
- **§04 TRIS maps:** plotting checks, how unseen sky is marked, whether the Haslam and TRIS masks match, the annulus,
  and nside 8 against nside 16.
- **§05:** every map in `res/` side by side, with frequencies, grids and beams (the documented check of
  `moonrunes.sky_maps.describe`).
- **§06:** the TRIS TOD each map predicts, against TRIS, scaled by the mean or by a power law, and β against RA (older
  native TODs, not the beam-matched ones).
- **§07:** the nside of every file loaded.

### Part B — Coordinate frames
- Checks what each input's header declares against the sky; pixel ordering; the rotation method per map (with a round
  trip); known sources after rotation; whether ARCADE 2's mask survives the rotation. The code is in
  `src/moonrunes/frames.py`.

### Part C — Stage-1 diagnostics
- **Stage 1:** has the solve converged, does the model describe the ring, the zero level and its degeneracy,
  prior → posterior shrinkage, and whether the product is fit for stage 3.
- **Prior-only runs:** what the data moved, where on the sky, high-noise runs, per-pixel σ and RMS maps, and
  residuals with error bars.
- The beam sections the introduction cell lists are in `notebooks/tests_analysis/` Part K2. Predates nside 16
  in places (TODO N-1). The original ran on the `21cmfast` kernel.

---

## `notebooks/tests_analysis/`

Parts read the pipeline's products (step 2's TODs, step 3's comparison, step 5's fits); H writes what I reads. Parts A, B, F, G and J are checks whose outputs go to `test_outputs/`.

| Part | Was | Reads | Writes | Feeds the paper |
|---|---|---|---|---|
| A — Beam-matched ring TODs, against the old TODs | `18_beam_matched_tods` | `paper_outputs/pipeline/step2_ring_tods/` (made by pipeline step 2) | `test_outputs/beam_matched_vs_old/` | — |
| B — Curved power-law fits and the 408 MHz sky models (a check) | `12_sky_models_vs_fit_408` | `paper_outputs/pipeline/step5_fit/` (power law), `paper_outputs/pipeline/step2_ring_tods/``sky_tods.npz`, `configs/model_comp_config.yaml` | `test_outputs/sky_models_408/` | — (the paper's numbers are steps 5–6) |
| C — Model comparison diagnostics | `11_model_comparison_diagnostics` | `paper_outputs/pipeline/step3_model_comparison/` (pipeline step 3) | `paper_outputs/model_comparison/diagnostic_plots/` | Section 5 counts (main.tex:501) |
| D — Model-choice decision sheet | `15_model_choice_decision_sheet` | `paper_outputs/pipeline/step3_model_comparison/`, `paper_outputs/pipeline/step5_fit/`, Part B's curved fits | `paper_outputs/model_choice/` | stacked residuals (main.tex:503) |
| E — T–T indices; sky models against TRIS | `14_ring_tt_betas` | `paper_outputs/pipeline/step2_ring_tods/`, `paper_outputs/pipeline/step5_fit/` | `paper_outputs/phil_plots/` | `tab:step1` (Plot 6a) |
| F — Haslam and the sky models against the power-law prediction (a check of step 6) | `16_haslam_vs_power_law_408` | `paper_outputs/pipeline/step2_ring_tods/`, `paper_outputs/pipeline/step5_fit/` | `test_outputs/plot3_haslam_vs_fit/` | main.tex:647 ratio (`tab:scale-offset` is step 6's) |
| G — Haslam with one gain and offset | `17_haslam_gain_offset_calibration` | `paper_outputs/pipeline/step2_ring_tods/`, `paper_outputs/pipeline/step5_fit/` | `test_outputs/plot3_haslam_calibration/` | possibly main.tex:647 (unconfirmed) |
| H — Zero-level offsets as fit parameters | `20_zero_level_offsets` | `paper_outputs/pipeline/step2_ring_tods/`, `paper_outputs/pipeline/step5_fit/` | `paper_outputs/zero_level_offsets/` | main.tex:394; input of I |
| I — Robustness of the Haslam gain and offset | `21_robustness_check` | `paper_outputs/pipeline/step2_ring_tods/`, `paper_outputs/pipeline/step5_fit/`, H's outputs | `paper_outputs/robustness_check/` | main.tex:396 and others |
| J — 13-point fits and the TOD prediction at any frequency (Tests 1–5; RHINO is pipeline step 7) | `19_power_law_tod_prediction` | `res/`, configs, `external/irfan_puglisi_2025/` | `test_outputs/tod_prediction/` | — |
| K1 — Beam widths and smoothing | `09_beam_fwhm_tests` | config, limTOD's beam | — | the 19.155° × 23.366° beam |
| K2 — Beam matching, harmonic checks, beam-specification variants | `03_diagnostics` cells 39–114 | stage-1 outputs, `res/`, B's cache | `test_outputs/beam_variants/` | — |
| K3 — Beam averaging with spectral-index blobs | `scripts/beam_averaging_blobs/` | (saved results) | — | main.tex:275 |
| K4 — Two-step vs one-step beam matching, four FWHM scalings | `scripts/fwhm_scaling_test.py` | (saved results in `test_outputs/fwhm_scaling_test/`) | — | main.tex:273 |
| K5 — K3 and K4 at LWA1 35 and Stockert 1420 MHz | `scripts/beam_tests_lwa35_stockert1420/` | (saved results) | — | — |

**Survey sets are not the same in every Part** (details in `HANDOVER.md`). The Parts that read the pipeline's
step 5 fits (B, D's fits, E, F, G, H, I) take its data sets: 13 points since the 2026-11-16 run (TRIS 600, TRIS
820, LWA1 35, LWA1 38–74 ×7, Maipu 45, LWA1 80, Stockert 1420), six before; their stored outputs are from six. C and
D's evidence read step 3 (11 points without TRIS since that run, four before). J chooses its own 13 points.

### Part A — Beam-matched ring TODs
- Reads the beam-matched TODs that pipeline step 2 makes (every survey and sky model: the `tris_cuts` kernel with the
  native beam taken out, at the ring pointings; the corrections stretched to the TRIS frequencies; LFSM extended with
  its own 300–408 MHz index) and compares them with the old native-beam TODs (Fig 1). The build moved to
  `moonrunes.pipeline.step_ring_tods` on 2026-11-16.

### Part B — Spectral fits and the 408 MHz sky models
- Reads the power-law fits (pipeline step 5) and the 408 MHz sky-model TODs (step 2), fits the curved power law per ring
  sample with UltraNest itself, and ranks every pair of sky model and fit by χ². A check of steps 5–6, which make the
  paper's numbers.
- Its Summary cell quotes χ²/N from a run of 2026-10-23 (GSM2016 134.0, Haslam + Carter 134.2) that differ from its own
  CSV and the paper (130.72, 137.11).

### Part C — Model comparison diagnostics
- Follows `comparison.statistic` and `comparison.exclude_tris` in `configs/model_comp_config.yaml` and the maps
  `sed.exclude_maps` leaves in (`configs/run_pipeline.yaml`), as pipeline step 3 does. Its own Parts 1–6:
  per-sample score differences; posteriors, best fits and the Occam split; the C prior; convergence; fits at
  600/820 MHz and χ²; robustness (Cygnus mask, Occam share, AIC/BIC, held-out prediction, stacked residuals,
  error-bar scaling).
- **It reads step 3's CSVs: since the 2026-11-16 pipeline run, 11 points without TRIS. Its stored outputs are from
  the 2026-10-29 run on four points.**

### Part D — Model-choice decision sheet
- The five criteria: ΔlnZ, χ²/AIC/BIC, leave-one-out prediction, stacked residuals and error-bar sensitivity, with
  and without the Cygnus windows. It does not choose a model.
- **Criterion 1 reads the four-point evidence from `paper_outputs/model_comparison/results/`; criteria 2–5 use B's six-point fits.**

### Part E — T–T indices and the sky models against TRIS
- **Plots 4–4f:** York fits for every survey pair, β against pair frequency or mean frequency, RA-split and
  whole-ring β matrices.
- **Plots 6–6d:** sky models and stretched corrections against TRIS at 600/820 MHz; models calibrated to TRIS
  (`TRIS = a·model + b`, weighted least squares, errors on TRIS only — `tab:step1`); the suggested GSM correction; the
  corrected GSMs against Haslam at 408 MHz.
- **Plot 5:** data β against model-implied β.

### Part F — Haslam against the power law at 408 MHz
- York scale–offset fit (`model = a·P + b`, errors on both axes) for Haslam and six sky models, with and without
  Cygnus (`tab:scale-offset`); the T–T plot, the scale–offset plane and the ratio along the ring. Every model is
  given Haslam's error (10% ⊕ 3 K).

### Part G — Haslam gain–offset calibration
- One gain g and offset z for Haslam against the power-law prediction, compared with the published corrections;
  the calibration plane and the ring view.

### Part H — Zero-level offsets
- One additive offset per survey shared by all samples, with Gaussian priors from the published zero-level
  systematics; Tables 1–3 (offsets, χ² per configuration, effect on the science); Figures A–C.

### Part I — Robustness of the Haslam gain and offset
- Haslam = g P + z (York), against the baseline (six-survey power law, Cygnus A mask ±25°, 103 samples): each survey
  left out, both LWA1 bands out, the best-measured three, TRIS only; power law and curved power law; Cygnus mask
  ±15°/±35° and the RA halves; three rows copied from Part H's Table 3(a). One table and the g–z figure.

### Part J — 13-point fits and the prediction at any frequency
- Power-law and curved fits per ring sample (`FIT_MODELS`, `DATASETS`), compared with the Irfan & Puglisi β.
- **Tests 1–5:** parameters along the ring, models against TRIS, χ² per frequency, held-out Haslam 408 MHz, predicted
  TODs at other frequencies; a prediction function for any frequency.
- **Test 6 (RHINO)** moved to the pipeline on 2026-11-16 (step 7, fed by our spectral model of step 5); a Markdown cell
  says so, and Test 6's last outputs (2026-11-05, 13-point fit) are in `test_outputs/tod_prediction/test6_*`.
- Its first cell still says "4 survey points … and both TRIS frequencies"; `DATASETS` has 13.

### Part K — Beam tests
- **K1:** the σ ↔ FWHM constant, an elliptical Gaussian through the FWHM function, degree/radian and σ/FWHM slips, the
  cut beam's numerical FWHM and principal planes, and a round trip on a known beam. `scripts/tests/test_beam_fwhm.py` repeats
  these checks.
- **K2:** the map-route elliptical beam (kernels, E/H profiles, a point source through the smoothing code, Haslam and
  ARCADE 2 under both kernels); harmonic-space checks H1–H6; convolve → deconvolve with limTOD (including the noise
  test along the ring); two-step against one-step beam matching (`PIPELINE.md` App. E); the beam-specification
  variants (roll, interpolation; Figs A–B, Tables 3–4, writes `test_outputs/beam_variants/`). Starts with a copy of the
  old notebook's setup cell, which has no stored output.
- **K3–K5:** the script-based tests, shown from their saved results and READMEs; rerun commands are in each Part and
  under "Standalone tests" in `README.md`.
