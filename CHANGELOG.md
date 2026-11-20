# Changelog

The decisions and changes that determine the current results, newest first. Each entry says what was decided, why,
where it is set, and which products it affects. Code history that no longer affects any result is not kept here.

## Decisions in force

### The Irfan correction split into an offset part and a spectrum part (2026-11-18)

Diagnostic model "<base>, Irfan spectrum only" = the base model's 408 MHz map (CMB-free RJ) × `irfan_scaling`, i.e.
the Irfan-corrected model without its −(8.9 K − T_CMB,RJ(408)) = −6.184 K offset. `pipeline.irfan_split` gives
total = corrected − base, offset part = corrected − spectrum only, spectrum part = spectrum only − base, and asserts
offset + spectrum = total; the offset part is also checked against −6.184 K × the beam-averaged scaling.
Products: step 2 `irfan_split_tris.csv` and `modeltris::<base>, Irfan spectrum only::<band>`; step 7
`rhino_checks.csv` ("Irfan split, …" columns) and the extra model in `rhino_tods_*.csv`. No existing value changes.

- Step 7, median along the ring, per cent of our prediction, 60 → 80 MHz: TRIS beam total −25.7 → −25.3
  (offset −20.4 → −20.9, spectrum −3.3 → −3.6); RHINO beam total −20.1 → −20.1 (offset −20.5 → −20.3,
  spectrum +2.7 → +0.1).
- Step 2, range along the ring at 600.5 / 817.8 MHz: GSM2008 total 0.752 / 0.668 K (offset 0.392 / 0.364, spectrum
  1.117 / 1.008, anticorrelated, r ≈ −0.97); GSM2016 total 1.923 / 1.499 K (offset 0.392 / 0.364, spectrum
  2.269 / 1.818, r ≈ −0.82).

### Step 7 reports the closure and the limits of its uncertainty (2026-11-18)

`rhino_checks.csv` gives the median closure |A x − d| (32.3, 22.1, 15.3, 10.7, 7.9 K at 60–80 MHz) and the
RHINO-beam σ/T three ways: as quoted, diag(R Σ Rᵀ) with independent errors (0.36–0.45 %); with the ring errors fully
correlated (0.82–0.94 %); with the prior errors fully correlated (3.6–4.0 %). The quoted value treats correlated
errors as independent and is a lower bound.

### T–T fits and Figure 1 on the pipeline's data sets (2026-11-18)

`tests_analysis` Part E and `10_paper_plots` (Figures 1, 5, 6) take their surveys from `pipeline.fit_data_sets`
(step 5's data sets) plus Haslam: 14 ring TODs and 90 pairs. LWA1 45 vs Maipu 45 is skipped (equal frequencies give
no β). Figure 1 shows the twelve maps in two columns.

### The Irfan correction is carried with Irfan & Puglisi's spectrum (2026-11-17)

At 408 MHz the correction is T_tot − 8.9 K (their eq. 6; −6.184 K on the CMB-subtracted maps). Off 408 MHz it is
carried with their synchrotron model, T(p, ν) = T(p, ν₀)(ν/ν₀)^(β_s(p) + c_s(p) ln(ν/ν₀)), ν₀ = 45 MHz (their eqs. 1,
2, 11), with their β_s and c_s maps (`external/`, licence not cleared). Set in `configs/run_pipeline.yaml`
`pipeline.sky_models.irfan_*`; code `pipeline.irfan_scaling`. Affects step 2 (TRIS frequencies) and step 7.
The Carter and GINAN corrections are carried with the base model's own per-pixel spectrum.

### The end-to-end pipeline (2026-11-16)

`run_pipeline.py --step all`: 0 TRIS maps, 1 beam-matched maps, 2 ring TODs, 3 model comparison (TRIS excluded),
4 selection, 5 our spectral model (the selected model fitted with TRIS), 6 the 408 MHz prediction, 7 the RHINO
prediction (`docs/PIPELINE.md`). Data sets: every map in `res/` that `sed.exclude_maps` leaves in with data on the
ring (11 surveys, including the LWA1 38–74 MHz maps), plus TRIS 600 and 820 in step 5 (13 points).

### Model comparison on the evidence, TRIS excluded (2026-10-28)

`configs/model_comp_config.yaml`: `comparison.statistic: evidence` (UltraNest ln Z summed over the ring samples),
`comparison.exclude_tris: true`; ties within `pipeline.select.tie_lnz` = 1 go to the model with fewer parameters.
The alternative `statistic: max_likelihood` (no Occam penalty) is implemented. The summed ln Z counts each
independent beam about six times (correlated samples), so it ranks models but overstates the preference.

### Ring TODs are beam-matched with the measured TRIS beam (2026-10-27)

`sed.tod_beam: beam_matched`: every map is observed at the 120 ring pointings with the `tris_cuts` kernel (the
archive's E/H cuts, 19.155° × 23.366°, E plane 7° east) with the map's own beam taken out; samples kept where the
observed beam weight is ≥ 0.5. `beam_matching.kernel: circular_gaussian` (23.366°) is used only by stage 2 and the
map route.

### Survey uncertainties from the publications (2026-10-21)

σ = scale·|T| ⊕ zero level ⊕ noise per survey (`configs/error_bars.csv`, `sed_prototype.PUBLISHED_ERRORS`):
LWA1 20 % (middle of 15–25 %), Maipu 1.7 % and 39 K (1σ of Monsalve et al. 2021's corrections, which are not
applied), Stockert 5 % and 0.5 K, CHIPASS 5 % and 0.5 K (Stockert's), Rhodes 5 % and 80 mK (upper limits), Haslam
10 % and 3 K. TRIS: statistical (floored at the ring's smallest positive value) ⊕ the larger side of the zero-level
uncertainty. Scale and zero-level terms are treated as independent per point.

### Maps left out of every fit (2026-10-15, Juno 2026-10-27)

`sed.exclude_maps`: ARCADE 2 (3.15, 3.41 GHz: ~7 % of the sky, 16–17 of 120 ring samples, and its 6–7.5 mK errors
let those few points pull the fits) and Juno 600 and 1248 MHz (unit, frame and beam not stated). CHIPASS and
Rhodes/HartRAO drop out because less than half of the beam falls on observed sky at every sample.

### Stage 1 (step 0) settings (2026-10-02 to 2026-10-14)

- Everything on nside 16, including `tris.nside_hires`. The cost, measured: reduced χ² 16.1 / 5.9 instead of
  3.09 / 3.22 with `nside_hires: 64`, and 23 / 7 instead of 388 / 201 pixels narrowed by more than 5 %.
- Pixels kept where the beam response exceeds 0.01 of its peak (1272 pixels, 41.4 % of the sky, 99.6 % of the beam
  power). Across thresholds 0.003–0.03 the maps agree to 0.40 σ; the fitted zero level moves by ~0.1 K.
- The solve is limTOD's (`TRISMapMakingInputs.solve`). The exact cross-check (`tris.solver.cross_check`) is off;
  when last on, it agreed to ~1e-10 σ.
- Zero level: one symmetric 0.5 K prior on both rings (the archive quotes +0.43 / −0.30 K at 820 MHz).
- Prior: GSM2008. The fitted zero levels (+1.74, +0.52 K) measure GSM2008's deficit along the ring.
