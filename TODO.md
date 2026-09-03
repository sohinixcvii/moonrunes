# TODO

Everything below came out of building stages 1 and 2 — either a decision that was
deliberately deferred, a substitute that stands where a real quantity should, or a
discrepancy found in `bayesian_skymap` while wiring against it. Items are ordered by what
they block, not by effort.

Cross-references: `tris_haslam_pipeline.md` for the stage definitions, the stage manifests
in `outputs/*/`, and `CHANGELOG.md` for what is already done.

---

## Blocking — results are not interpretable until these are resolved

- [ ] **B1 — Obtain the Berkhuijsen (1972) 820 MHz survey map.**
  *Why:* the only genuinely Haslam-independent prior basis. Everything else in the
  pipeline is built and waiting on it.
  *Where:* set `paths.berkhuijsen_map`, then `haslam.prior.source: berkhuijsen_1972`.
  The reader, the 820 → 408 MHz extrapolation and the manifest flag already exist.
  *Done when:* `outputs/stage2/stage2_manifest.json` reports
  `prior_is_haslam_derived: false` with `source: berkhuijsen_1972`.

- [ ] **B2 — Replace the GSM2008 prior in stage 1 as well.**
  *Why:* stage 1 uses GSM2008 as *both* the prior template and the reference truth
  (`prior_equals_truth: true`), so its maps cannot be checked against anything
  independent. Stage 2 already has three prior sources; stage 1 has one.
  *Where:* `tris.prior.template` currently only implements `gsm2008` — `stage1_tris_maps.py`
  raises on anything else.
  *Done when:* stage 1 accepts the same source switch stage 2 has, and its manifest can
  report `prior_equals_truth: false`.

- [ ] **B3 — Decide what to do about `gain`, which cannot exist for real data.**
  *Why:* `gibbs_spectral_nside.py` reads `data["gain"]` — the **true** gain — and uses it
  for `covG` (line 89) and as `mu_g`, the prior mean of the gain. There is no true gain for
  Haslam. Stage 2 emits a `covG_diag` substitute evaluated at `g = 0`, but nothing yet
  supplies `mu_g`.
  *Done when:* stage 4 passes stage 2's `covG_diag` and an explicit `mu_g` (presumably
  zero) instead of reaching into the npz for a truth.

---

## Stage 3 owes (assembly only — it should invent nothing)

- [ ] **S3-1 — Rotate stage 1's TRIS maps into galactic.**
  Stage 1 works equatorial, stage 2 galactic (`haslam.frame`). Stage 2's
  `_tris_at_408()` already does exactly this resampling (rotate each galactic pixel centre
  into equatorial, take the stage 1 pixel it lands in — nearest neighbour, no smearing at
  the band edge); factor it out rather than writing it twice.

- [ ] **S3-2 — Decide how the 68 partially-covered coarse pixels are filled.**
  Degrading the TRIS band from nside 64 to `beam.nside_new = 8` touches 432 of 768 pixels:
  364 fully covered, **68 partial**, mixing solved sky with pixels the ring never reached.
  `notebooks/03` previews this but explicitly leaves the decision to stage 3.

- [ ] **S3-3 — Build `pixel_mask` at `nside_new` from the rotated footprint.**
  `bayesian_skymap` uses it to zero out unobserved low-resolution pixels. It must match
  whatever S3-2 decides.

- [ ] **S3-4 — Write `sky_haslam` from stage 2's `s0_init_k`, not a truth.**
  The scripts read that key unconditionally as `true_s`. `assemble.sky_haslam_is_substitute`
  is already in the config; the manifest must record the substitution so nothing downstream
  mistakes it for a truth map.

- [ ] **S3-5 — Assemble `freq`, `ref_freq`, `freq_beam`, `lmax_beam`.**
  `freq` must have 408 MHz **last** — `precompute_conv` assumes `ref_freq == freq[-1]` —
  and `ref_freq` must be an element of `freq`, because the scripts look it up with
  `np.where`. With two TRIS channels this gives `sky_lowres` shape `(2, npix_new)`.
  Note the tension: `tris.ref_freq_mhz` is 600, but the Haslam channel is last; check which
  the scripts actually require before setting it.

---

## Stage 4 prerequisites (patches to `bayesian_skymap`)

- [ ] **S4-1 — `iter=0` on the `hp.smoothing` calls. There are THREE, not two.**
  The pipeline document says "both `hp.smoothing` calls". The fork actually has three:
  `bayesian_func.py:472` (`transpose_beam`), `:483` (`apply_beam`), and `:693`
  (`precompute_conv`). Patching only the first two leaves the asymmetry that makes the
  solver diverge in the third path.

- [ ] **S4-2 — Force `gmres`, never `cgs`.**
  `cgs` is documented to diverge on this operator. `gibbs.solver: gmres` is already the
  recorded decision; the script still has `cgs` branches at lines 223 and 276.

- [ ] **S4-3 — Keep the Woodbury fast path off (Blocker 3).**
  `gibbs.use_woodbury: false`. It drops cross-frequency coupling and is wrong by ~8 orders
  of magnitude at 2 channels, which is exactly the TRIS regime.

- [ ] **S4-4 — Feed stage 2's `covA` in instead of the hardcoded truth prior.**
  `gibbs_spectral_nside.py:158` hardcodes `covA = (0.1 * true_s)**2`. Stage 2 produces the
  real `covA`; it currently has no way in.

- [ ] **S4-5 — Decide whether `sprior_mean` stays 0.**
  Line 160 sets it to zero, so only `covA` is consumed and stage 2's `prior_mean_k` is
  ignored. That may be right, but it should be a decision rather than an inherited default —
  a zero-mean prior on a strictly positive sky is a strong statement.

- [ ] **S4-6 — Confirm `gibbs.nsample: 2500` is intended.**
  Changed from 25 in the working tree and not yet committed. At 2500 emcee samples per
  Gibbs iteration × 200 iterations this is a large run; worth a timing estimate first.

---

## Open decisions carried from stage 1

- [ ] **S1-1 — Reconsider `tris.nside: 64`.**
  At nside 64 the band holds 25704 pixels against 120 samples, and the map is
  ~entirely prior: median posterior/prior σ ratio **0.99998**, with only **2 of 25704**
  pixels tightened by more than 5%. Since stage 3 degrades to nside 8 anyway, solving
  directly at nside 8 (400 band pixels) would be the more honest grid. One config line;
  deliberately left alone because it is a recorded decision.

- [ ] **S1-2 — Model the asymmetric 820 MHz zero level properly.**
  The archive quotes **+0.430 / −0.300 K**; stage 1 covers both rings with one symmetric
  0.5 K prior (`tris.zero_level_sigma_k`) rather than silently symmetrising. A two-sided or
  skewed treatment would use the published numbers as published.

- [ ] **S1-3 — The E/H beam asymmetry is lost in the recalibration forward model.**
  Stage 1 carries the real 19.155° / 23.366° beam. `bayesian_skymap`'s forward model is
  azimuthally symmetric (one `fwhm` per `hp.smoothing`), so `beam.beam_deg` is the E/H mean
  (21.2605°). Blocker 1 is resolved in the sense that the grid is consistent, but the
  asymmetry itself cannot be represented downstream. Quantify the error this introduces.

- [ ] **S1-4 — `bicgstab` stalls on the stage 1 operator; document or drop it.**
  It is in `_SOLVERS` and passes the well-conditioned unit test, but on the real nside-8
  problem it returns `info=1000` at every tolerance tried and sits 5.5e-2 posterior σ from
  the exact answer. The cross-check catches it, but the option invites a bad run.

- [ ] **S1-5 — `tris.mapmaker: beamwidth_binning` is not wired up.**
  Stage 1 raises on it. `limTOD.tris.tris_beamwidth_binning_map` exists if the simplified
  scheme is ever wanted as a comparison.

- [ ] **S1-6 — Revisit the 2.5 GHz channel.**
  Dropped as a point set (6 points, not a ring) and flagged poor quality. Stage 1 refuses it
  with an explicit error. Worth confirming that judgement against the literature review
  rather than inheriting it.

---

## Open decisions carried from stage 2

- [ ] **S2-1 — Reconsider `region_geometry: theta_bands`.**
  Equal-spaced colatitude bands leave a **7.8×** pixel-count imbalance across the 12 gain
  regions (3280 in the polar bands vs 25600 in the central ones), so the polar gain
  parameters are far less constrained. `equal_area` is implemented and reduces the
  imbalance to <1.2×; `theta_bands` is the default only because it is `bayesian_skymap`'s
  own convention.

- [ ] **S2-2 — Sanity-check `n_gain_regions: 12` against 2 TRIS frequencies.**
  12 gain regions were chosen to match `bayesian_skymap`'s generator, but that was for a
  ~20-channel synthetic problem. With 2 usable TRIS channels, 12 free gain parameters may
  be more than the data supports — worth a rank or degeneracy check before stage 4.

- [ ] **S2-3 — `op_alm = 1` (spherical-harmonic gain) is unimplemented.**
  Stage 2 raises. That parametrisation builds no operator matrix at all, so supporting it
  means deciding what stage 2 produces instead, not just relaxing a check.

- [ ] **S2-4 — Verify the Haslam zero-level and scale numbers used downstream.**
  `gibbs.sig2_mk: 910` comes from the EDGES-based ~0.91 K zero-level estimate on the Sky Map
  Priors page, and stage 5 checks against ~3% temperature scale. Both are quoted from
  memory of that page and should be pinned to the actual reference.

---

## Housekeeping

- [ ] **H1 — Create the environment the repo describes.**
  `environment.yml` defines `tris-haslam-cal`; all work so far has run in the pre-existing
  `limtod` env with `limTOD` imported from the sibling checkout via `paths.limtod`. The
  documented setup path is therefore untested.

- [ ] **H2 — `jupyterlab`/`nbclient` are not in the env used to run the notebooks.**
  The committed notebooks had to be executed with another environment's `nbclient` driving
  the `limtod` kernel. Anyone re-running them from `environment.yml` should be fine, but
  that is exactly what H1 has not verified.

- [ ] **H3 — `build_s0_init()` takes `prior_mean_k` and never uses it.**
  `stage2_haslam_prep.py:369`. Harmless, but it implies a dependency that is not there.

- [ ] **H4 — Add a stage 2 diagnostics section to `notebooks/03`, or split it.**
  `03_diagnostics.ipynb` is stage 1 only. Stage 2's diagnostics currently live in
  `notebooks/02`, which is really the prior-construction notebook.

- [ ] **H5 — Decide the CI story for the end-to-end tests.**
  6 of the 26 tests skip themselves without the TRIS archive, the cached Haslam map, or
  `pygdsm`. The algebra tests run anywhere; the pipeline tests only run where the data is.
