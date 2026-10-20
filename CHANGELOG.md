# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Stages refer to `tris_haslam_pipeline.md`, and the blockers to its "three blockers to
resolve before running anything".

## [Unreleased]

### Added

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
