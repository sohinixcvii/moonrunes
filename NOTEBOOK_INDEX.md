# moonrunes — notebook index

As of 2026-11-05. What each notebook in `notebooks/` does, what it reads and where it writes. For the individual plots
see `PLOT_INDEX.md`; for the paper steps see `../PAPER_OUTLINE.md`.

**Status** column:
- **current**: runs on the pipeline as it is now.
- **stale**: its saved outputs predate a later change, so rerun it before using its numbers.
- **superseded**: kept for reading only.

"Notebook 12's cache" is `model_comparison/sky_models_408/cache/`. Its fits (`fits/unit###__*.npz`, 23 Oct) and
`sky_tods.npz` predate the beam-matched TODs. Every notebook that reads this cache inherits that.

## At a glance

| # | Notebook | What it is for | Writes to | Status |
|---|---|---|---|---|
| 01 | `01_explore_data` | First look at every input: TRIS rings, ARCADE 2, Haslam, TRIS maps, every map in `res/`, and the TRIS TOD each map predicts | nothing (inline plots) | current (§06 uses the older TODs) |
| 02 | `02_build_berkhijsen_prior` | The old Gibbs-sampler prior work: GSM2008/Haslam circularity, Berkhuijsen prior | nothing | superseded: imports a deleted module, do not run |
| 03 | `03_diagnostics` | Stage 1 diagnostics, then beam, harmonic and beam-matching tests | inline; `outputs/beam_variants/` (beam-variant section) | current; the "Along the TRIS ring" noise test still shows the 10% errors (N-14) |
| 04 | `04_coordinate_frames` | Puts Haslam, ARCADE 2 and TRIS into one (equatorial) frame and checks it | nothing | current |
| 05 | `05_beam_matching` | Stage 2 beam matching with the circular 23.366° kernel; overlap at nside 8 | nothing | superseded by the elliptical `tris_cuts` kernel (03, 18) |
| 06 | `06_sed_fit` | Per-pixel power-law SED from TRIS + ARCADE 2 only, compared with Haslam at 408 MHz | nothing | current for what it does |
| 07 | `07_sed_calibration` | Haslam calibration prototype: TOD route against map route | nothing | stale: 10% errors (N-14) |
| 08 | `08_update_plots` | The figure set for the SED calibration of Haslam (Figs 1–16) | `outputs/paper_plots/` | Figs 1–13b stale (N-14); Figs 14–16 current |
| 09 | `09_beam_fwhm_tests` | Checks the beam FWHM numbers (19.155° / 23.366°) and the code that makes and uses them | nothing | current |
| 10 | `10_paper_plots` | The paper's figures, to the paper's conventions | `Papers_and_reports/moonrunes_paper/figures/` | current inputs; needs a full rerun (N-14); Figure 8 run 2026-11-05; Figures 9–10 run 2026-11-06 (after notebook 19) |
| 11 | `11_model_comparison_diagnostics` | Diagnostics of the spectral-model comparison (Parts 1–6, Figs 1–17) | `model_comparison/diagnostic_plots/` | stale: saved figures predate the current `data.datasets` |
| 12 | `12_sky_models_vs_fit_408` | Power-law and curved fits per ring sample (UltraNest), and which 408 MHz sky model matches the prediction | `model_comparison/sky_models_408/` (+ the shared cache) | stale: fit cache from 23 Oct |
| 13 | `13_sky_model_residuals_408` | Plots 1–3 for Phil: sky models against the power-law 408 MHz prediction | `model_comparison/phil_plots/` | stale (reads 12's cache) |
| 14 | `14_ring_tt_betas` | T–T β along the ring for every survey pair; sky models against TRIS at 600/820 MHz (Plots 4–6d) | `model_comparison/phil_plots/` | current (beam-matched) |
| 15 | `15_model_choice_decision_sheet` | The model-choice decision sheet: 5 criteria, with and without Cygnus; the reader applies the rule | `model_comparison/model_choice/` | stale (reads 12's cache) |
| 16 | `16_haslam_vs_power_law_408` | Haslam and the other sky models against the power-law prediction: scale–offset fits | `model_comparison/plot3_haslam_vs_fit/` | stale (reads 12's cache) |
| 17 | `17_haslam_gain_offset_calibration` | Phil's plot 3: Haslam calibrated with one gain and offset along the ring | `model_comparison/plot3_haslam_calibration/` | stale (reads 12's cache) |
| 18 | `18_beam_matched_tods` | Builds the beam-matched ring TODs and measures the change from the old ones | `model_comparison/beam_matched_tods/` | current |
| 19 | `19_power_law_tod_prediction` | Power-law and curved fits per ring sample, and the predicted ring TOD at any frequency (Tests 1–5); RHINO predictions at 60–80 MHz (Test 6) | `tod_prediction/` | current |
| 20 | `20_zero_level_offsets` | Zero-level offsets as shared fit parameters, with priors from the published systematics (TRIS II Table 3 for TRIS) | `model_comparison/zero_level_offsets/` | current; reads 12's cache for the TRIS rings |
| 21 | `21_robustness_check` | Is the Haslam gain and offset robust? A fixed grid of survey subsets, spectral models, sample subsets and zero-level treatment | `model_comparison/robustness_check/` | current; reads 12's cache for the TRIS rings and 20's Table 3(a) |

## What each notebook can do

### 01 — The data the pipeline consumes
- **§01 TRIS:** the four archive products, where `tris.zero_level_sigma_k = 0.5` comes from, the temperature convention,
  the beam, the GSM2008 prior template, and what 120 samples can support.
- **§02–03:** ARCADE 2 at 3.15 and 3.41 GHz, including how its unobserved pixels are marked; Haslam 408 MHz.
- **§04 TRIS maps:** plotting checks, how unseen sky is marked, whether the Haslam and TRIS masks match, the annulus,
  and nside 8 against nside 16.
- **§05:** every map in `res/` side by side, with frequencies, grids and beams.
- **§06:** the TRIS TOD each map predicts, against TRIS, scaled by the mean or by a power law, and β against RA.
- **§07:** the nside of every file loaded.

### 02 — The Haslam side and the prior problem (superseded)
- The region operators; why the GSM2008 prior is circular (log-space r = 0.993); the independent Berkhuijsen prior.
  Read only.

### 03 — Diagnostics
- **Stage 1 (cells 1–18):** has the solve converged, does the model describe the ring, the zero level and its
  degeneracy, prior → posterior shrinkage, and whether the product is fit for stage 3.
- **Prior-only runs (19–38):** what the data moved, where on the sky, high-noise runs, per-pixel σ and RMS maps, and
  residuals with error bars.
- **Elliptical beam for the map route (39–60):** the three kernels, E/H profiles, the kernel each survey gets, a point
  source through the smoothing code, Haslam and ARCADE 2 under both kernels, and every map-route map.
- **Harmonic-space checks H1–H6 (61–76):** C_ℓ of every map, the TRIS beam's C_ℓ, limTOD pointing against the exact
  beam, the effect on the TODs, Haslam's harmonic rotation, and harmonic against real-space smoothing.
- **Convolve → deconvolve (77–91):** which scales come back, every map, and the effect along the ring with the
  pipeline's noise.
- **Two-step against one-step beam matching (92–99).**
- **Beam variants (100–114):** how much T against RA depends on the beam's roll and on how the E/H cuts are
  interpolated; the roll check, Fig A, Fig B, and Tables 3–4. Writes to `outputs/beam_variants/`, with its TOD cache
  in `cache/`.

### 04 — Coordinate frames
- Checks what each input's header declares against the sky.
- Covers pixel ordering, the rotation method per map (with a round trip), known sources after rotation, and whether
  ARCADE 2's mask survives the rotation.
- The code is in `src/moonrunes/frames.py`.

### 05 — Beam matching (circular kernel)
- `smooth_to_target` on Haslam and ARCADE 2 in the equatorial frame, with sanity checks, before/after maps and the
  usable overlap.

### 06 — Per-pixel SED fit (TRIS + ARCADE 2)
- A power law per TRIS-stripe pixel from 2–4 points, evaluated at 408 MHz against Haslam, with χ² where ARCADE 2 is
  present.
- Includes sensitivity to the two upstream choices. No correction is applied.

### 07 — SED calibration of Haslam
- **TOD route:** a power law per ring sample, β against LST with points left out, and the corrected TOD turned into
  a map with the Wiener filter.
- **Map route:** a power law per pixel.
- Compares the two routes with each other and with the benchmarks.

### 08 — Update plots
- **Figs 1–4:** the maps used, the target and the TRIS stripe, slices along the stripe, and the TRIS maps with their
  TODs.
- **Figs 5–7g:** SED fits per pixel and per sample, best fit − data, best-fit Haslam from both routes, and every map
  scaled to 408 MHz.
- **Figs 8–11:** residuals, χ² and β maps, and the corrected Haslam at full resolution.
- **Figs 12–13b:** T–T plots and β matrices.
- **Figs 14–14b:** the TRIS beam (dB map with E/H crosshairs and cuts), and the beam pointed by limTOD.
- **Figs 15–16:** how strongly TRIS anchors the Haslam gain and offset. Haslam = g P + z for three predictions P
  (TRIS only, without TRIS, all six surveys), the TRIS-only fit against a fixed extrapolation β, and the Fisher
  σ(ln A at 408 MHz) with each survey left out. Tables `fig15_part*`, `fig16_part*`.

### 09 — Beam FWHM tests
- **Unit and formula checks:** the σ ↔ FWHM constant, an elliptical Gaussian through the FWHM function, and checks
  for degree/radian and σ/FWHM slips.
- **The cut beam:** its numerical FWHM, φ = 0 / 90° giving the E/H cuts, and a round trip on a known beam.

### 10 — Paper plots
- Eleven figures written to the paper's `figures/` folder (PDF + PNG; the paper includes the PDFs) at journal widths, ten
  of them in `main.tex` (not `rhino_tods_tris_beam`), named `<section>_<content>` after the paper section (since 2026-11-06):
  - all maps over the TRIS slice;
  - the TRIS maps;
  - the sky-model TODs at 600/820 MHz;
  - the β matrix;
  - β against frequency;
  - model-comparison fits at the TRIS frequencies;
  - the TRIS beam;
  - how strongly TRIS anchors the Haslam gain and offset, the paper version of 08 Fig 15. Its
    cell redoes 08's fits for the g–z covariances and checks them against 08's CSVs, so run 08 first.
  - the predicted TODs at 60–80 MHz against the sky models, through the TRIS beam and through the RHINO beam, the
    paper versions of notebook 19's Test 6 figures. They read 19's CSVs, so run 19 first. Only the RHINO-beam one is
    in the paper.
  - the 408 MHz models against the six-survey power-law prediction along the ring, through the TRIS beam (Figure 11,
    added 2026-11-06). It reads notebook 12's fit cache and notebook 18's TODs, and needs only cell 1.
- Each figure has a "Featured in" note.

### 11 — Model comparison diagnostics
- Follows `comparison.statistic` and `data.datasets` in `model_comparison/config.yaml`.
- **Part 1:** per-sample score differences, absolute scores, distributions and win counts.
- **Part 2:** posteriors, best fits and the Occam split.
- **Part 3:** sensitivity to the C prior.
- **Part 4:** convergence (UltraNest against the exact integral, seeds, live points).
- **Part 5:** fits at 600/820 MHz and χ² of every fit.
- **Part 6:** robustness (Cygnus mask, Occam share, AIC/BIC, held-out prediction, stacked residuals, error-bar
  scaling).

### 12 — Sky models against the fitted 408 MHz sky
- Fits the power law and curved power law per ring sample with UltraNest, keeping the posteriors.
- Builds 408 MHz TODs for the sky models (GSM08/16, LFSM, Carter, Irfan, GINAN).
- Ranks every pair of sky model and fit by χ².
- Its cache feeds 13, 15, 16, 17 and 20.

### 13 — Sky-model residuals at 408 MHz
- Plots 1–3: legacy and corrected models against the power-law prediction, and χ² per model.
- Reads 12's cache only.

### 14 — T–T β along the ring
- **Plots 4–4f:** York fits for every survey pair, β against pair frequency or mean frequency, and RA-split and
  whole-ring β matrices.
- **Plots 6–6d:** sky models and stretched corrections against TRIS at 600/820 MHz, models calibrated to TRIS, the
  suggested GSM correction, and the corrected GSMs against Haslam at 408 MHz.
- **Plot 5:** data β against model-implied β.

### 15 — Model-choice decision sheet
- **The five criteria:** ΔlnZ, χ²/AIC/BIC, leave-one-out prediction, stacked residuals and error-bar sensitivity.
- **Figures:** the decision sheet, ΔlnZ along RA with the Cygnus windows, and a summary table.
- It does not choose a model, and stops if 12's cache is stale.

### 16 — Haslam against the power law at 408 MHz
- York scale–offset fit (T_model = a·T_PL + b) for Haslam and six sky models, with and without Cygnus.
- The T–T plot, the scale–offset plane, and the ratio along the ring.

### 17 — Haslam gain–offset calibration
- One gain g and offset z for Haslam against the power-law prediction, compared with the published corrections.
- Shows the calibration plane and the ring view.

### 18 — Beam-matched ring TODs
- Builds the beam-matched TODs for every survey and sky model (the `tris_cuts` kernel with the native beam taken out,
  nside 64 sky), the same as `sed_prototype.Inputs.ring_tod`.
- Fig 1 shows the beam-matched minus old TODs.

### 19 — Power-law TOD prediction
- Power-law and curved fits per ring sample (`FIT_MODELS`, `DATASETS` in the settings cell), compared with the Irfan
  & Puglisi β.
- **Tests 1–5:** parameters along the ring, models against TRIS (residuals in K), χ² per frequency, held-out Haslam
  408 MHz, and predicted TODs at other frequencies.
- A prediction function for any frequency.
- **Test 6 (RHINO):** the power-law prediction at 60, 65, 70, 75 and 80 MHz, deconvolved from the TRIS beam (Wiener
  solve, stage 1's prior width, GSM2016 prior) and observed with the RHINO
  `DryHornGround` beam (`../rhino-beams/npz_beams/`) at the Jodrell Bank zenith. Two figures, before and after the
  RHINO beam, one panel per frequency, against GSM2008, GSM2016, LFSM, GSM2008 + Carter and GSM2008 + Irfan (carried
  from 408 MHz with GSM2008's spectrum) and GSM2016 + GINAN through the same beams.

### 20 — Zero-level offsets
- One additive offset per survey shared by all samples, with Gaussian priors from the published zero-level
  systematics. The per-sample fits are nested inside the offset fit.
- **Tables 1–3:** the offsets, χ² per configuration, and the effect on the science.
- **Figures A–C:** the profile of F, example SEDs, and stacked residuals per survey.

### 21 — Robustness check of the Haslam gain and offset
- Haslam = g P + z (York, as 08 Figs 15–16), with P the spectral fit's A at 408 MHz, against the baseline (six-survey
  power law, Cygnus A mask ±25°, 103 samples).
- **The grid:** each survey left out, both LWA1 bands out, the best-measured three, TRIS only; power law and (with ≥ 4
  surveys) curved power law; Cygnus A mask ±15° / ±35° and the RA halves split at 200°; three rows copied from 20's
  Table 3(a).
- Counts the spectral fits that do not converge. One table (`robustness_table.csv`), differences by subtraction, and
  the two-panel g–z figure.

## Not notebooks: standalone test folders

| Folder | What | Script |
|---|---|---|
| `outputs/fwhm_scaling_test/` | FWHM-scaling test on Haslam (two-step vs one-step, 4 width rules) | `scripts/fwhm_scaling_test.py` |
| `beam_averaging_blobs/HASLAM/`, `HASLAM_DESTRIPED_ONLY/` | β-blob beam averaging on Haslam | `beam_averaging_blobs/beam_averaging_blobs.py`, `replot_from_csv.py` |
| `beam_tests_lwa35_stockert1420/fwhm_scaling/` | the FWHM-scaling test on LWA1 35 MHz and Stockert 1420 MHz | `beam_tests_lwa35_stockert1420/fwhm_scaling_test_maps.py` |
| `beam_tests_lwa35_stockert1420/beam_averaging_blobs/` | the β-blob test on LWA1 35 MHz and Stockert 1420 MHz | `beam_tests_lwa35_stockert1420/beam_averaging_blobs_maps.py` |
| `model_comparison/results/` | the model comparison itself (evidence / max likelihood per ring sample) | `model_comparison/run_comparison.py` (run manually) |
