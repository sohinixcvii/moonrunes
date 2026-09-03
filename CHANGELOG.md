# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Stages refer to `tris_haslam_pipeline.md`, and the blockers to its "three blockers to
resolve before running anything".

## [Unreleased]

### Added

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
  the Woodbury acceleration that is wrong by eight orders of magnitude at two channels is
  never on that path.
- 🔍 **A convergence check that is not "did it finish"** — `woodbury_reference` computes
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
  refuses to start if the pinned `beam.nside_new = 8` has drifted from it, so the grid
  stage 3 hands to `bayesian_skymap` cannot move under the pipeline unnoticed.
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
  expects), and a preview of the degrade to `nside_new = 8` (432 of 768 pixels touched,
  364 fully covered).
- ✅ **Tests** (`tests/test_stage_io.py`, 14 of them). The algebra is pinned against dense
  references — the Krylov solve for all three methods, the Woodbury posterior including
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
