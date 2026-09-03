# TRIS → Haslam Calibration: Literal Pipeline
### Combining `limTOD.tris` (map-making) with `bayesian_skymap` (recalibration)

**Primary aim:** make maps from TRIS TOD data, use them to calibrate the
Haslam maps. This document is the concrete algorithm, stage by stage,
using the actual functions/keys from both repos — not a restatement of
the READMEs.

---

## The core insight tying the two repos together

`bayesian_skymap`'s data model needs exactly two things: a high-res,
badly-calibrated map (`data2`/`sky_gain` — Haslam) and a low-res but
**absolutely calibrated** map (`data1`/`sky_lowres`). TRIS *is* an
absolute radiometer. So:

```
limTOD.tris  (Stage 1)  produces  →  calibrated TRIS sky map(s)
                                       ↓
                              becomes bayesian_skymap's `data1`/`sky_lowres`
                                       ↓
bayesian_skymap (Stage 2-4)  jointly infers  →  true sky (s), Haslam gain
                                                  corrections (g), spectral
                                                  index (β)
```

---

## ⚠️ Three blockers to resolve BEFORE running anything

These aren't implementation detail — each one changes what the pipeline
actually does if left unresolved.

### Blocker 1 — beam width mismatch -- NOW RESOLVED
`bayesian_skymap` hard-codes `beam_deg ∈ {10, 30}` only (any other value
now raises immediately, per the recent fix). TRIS's real beam is
**19.155° (E) / 23.366° (H)** — neither value. **This needs a decision:**
extend `bayesian_skymap` to accept a `beam_deg≈20` → new `nside_new`
mapping (real code change, needs care since `nside_new` propagates through
the Woodbury solver's shapes), or approximate TRIS's beam as the nearer
of the two fixed options (30°, accepting a real mismatch between the
beam actually used to make the TRIS map and the beam `bayesian_skymap`
assumes when forward-modelling it).

### Blocker 2 — no "truth" exists for real data
`bayesian_skymap`'s scripts expect `nstart=0` to **"initialise `s0` at the
true sky"** — fine for the paper's synthetic validation, meaningless for
real TRIS+Haslam data, since there is no true sky to initialise from.
Likewise the fixed prior `covA = (0.1 * true_s)**2` is truth-informed by
construction. **Both need a real-data substitute before this can run on
anything but synthetic data:**
- `s0` initialisation: use Haslam itself (crude but not circular for
  *initialisation*, only for the *prior* — see next point), or the
  TRIS map extrapolated across pixels.
- The prior `covA`: **this is exactly what the Sky Map Priors research
  already solved** — use Berkhuijsen (1972, 820 MHz) as the
  Haslam-independent basis for a real prior covariance, not
  `0.1 * true_s`. This is the single most important connection between
  your two active research threads right now — the prior-sourcing
  problem `bayesian_skymap` has is the exact problem that research
  answered.

### Blocker 3 — Woodbury fast path is documented-broken at low channel count
The Woodbury acceleration (`solve_ally_parallel` etc.) is verified correct
with ~20 low-res channels but **wrong by 8 orders of magnitude with 2
channels** (cross-frequency coupling dropped). TRIS realistically gives
**2-3 usable frequencies** (600, 820 MHz; 2.5 GHz flagged as poor quality
in your own TRIS literature review). **Decision: do not use the Woodbury
fast path for this application — use the direct/brute Krylov solve
instead.** Compute cost isn't a real concern at 2-3 channels anyway, so
there's no reason to use the accelerated path that's documented wrong in
exactly this regime.

---

## Stage 1 — Produce calibrated TRIS maps (`limTOD.tris`)

1. Read raw TRIS archive data via `limTOD.tris`'s strict readers, for
   600 MHz and 820 MHz (drop 2.5 GHz per the known quality issue, unless
   there's a specific reason to include it).
2. Build the beam from the archive's own E/H cuts (already verified in
   your walkthrough: HPBW 19.155°/23.366°, `selfrot=-7°`).
3. Apply the verified zenith/roll geometry (`RA=LST`, `dec=lat=42°26′`).
4. Run map-making. Two real options, worth deciding rather than
   defaulting: `limTOD.tris`'s own prior-regularized map-maker (the one
   from your walkthrough notebook — has the known zero-level/monopole
   degeneracy), or the simplified beam-width scheme already implemented.
   **Given Blocker 2's prior problem, the prior-regularized map-maker is
   probably the right choice here** — it's the one built to consume an
   external prior, which is exactly what the Berkhuijsen-based prior from
   Blocker 2 is for.
5. **Output:** calibrated TRIS maps at each frequency, at whatever
   `nside` Stage 1 produces (needs reconciling against Blocker 1's
   `nside_new` decision — these must match before Stage 3).

## Stage 2 — Prepare Haslam (`bayesian_skymap`'s `data2`)

1. Load the Jodrell Bank reprocessed Haslam map, regrid to `nside=128`
   if not already (this is `bayesian_skymap`'s fixed full-resolution
   grid).
2. Build the `operator` matrix (gain-parameter → pixel mapping). Decide
   `op_alm=0` (region-based) vs `op_alm=1` (spherical harmonic,
   `lmax_gain=5`) — genuinely different parametrisations of the same
   unknown, worth picking deliberately rather than defaulting.
3. Build the Berkhuijsen-based prior covariance (resolves Blocker 2).

## Stage 3 — Assemble `bayesian_skymap` inputs

Build the `.npz` with the required keys, substituting real TRIS/Haslam
data for the synthetic-test keys the repo's own data-generation script
would normally produce (that script isn't in the repo — you're building
this file by hand):
- `sky_lowres` ← Stage 1's TRIS maps, shape `(Nfreq-1, npix_new)`
- `sky_gain` ← Stage 2's Haslam map
- `operator` ← Stage 2's gain-parameter mapping
- `freq`, `ref_freq`, `freq_beam` ← TRIS's actual frequencies (ref_freq
  likely 600 or 820 MHz, whichever has the better-quantified absolute
  calibration per your TRIS literature review)
- No `sky_haslam` "truth" key available — see Blocker 2 for the
  initialisation substitute

## Stage 4 — Run the Gibbs sampler

1. Use `gibbs_nside.py` if fixing β, or `gibbs_spectral_nside.py` if
   jointly inferring the spectral index too (recommended, given
   synchrotron β is genuinely uncertain, not a nuisance to fix).
2. Set `sig1` (TRIS map noise) from your literature review's quoted
   TRIS precision (66 mK systematic + 18 mK statistical at 600 MHz).
   Set `sig2` (Haslam) from a real calibration-uncertainty estimate
   (your Sky Map Priors page already has one: Haslam's temperature scale
   uncertain at ~3%, zero-level at ~0.91 K per the EDGES-based analysis
   quoted there).
3. Apply the two easy, already-documented fixes before running: pass
   `iter=0` to both `hp.smoothing` calls (fixes the asymmetry causing
   `cgs` to diverge), and force `gmres` (or `cg` with `iter=0`) rather
   than the default `cgs`, which is documented to diverge on exactly this
   operator.
4. Do **not** use the Woodbury fast path (Blocker 3).

## Stage 5 — Validate

1. Check the log-posterior trace and per-pixel χ² for convergence — not
   just "did it finish," given `bayesian_skymap`'s own documented history
   of silently-divergent solves being accepted as valid samples.
2. Compare recovered Haslam gain corrections `g` against the literature
   benchmark already in hand (~3% temperature scale, ~0.91 K zero-level,
   from the Sky Map Priors page) — this is your actual sanity check that
   the pipeline is doing something real, not just running.
3. Sanity-check recovered β against the expected Galactic synchrotron
   range (~−2.5 to −3.0) — `bayesian_skymap` itself initialises 6 regions
   spanning exactly this range, so a wildly different recovered value is
   a red flag, not a discovery.

---

## What to actually decide before Stage 1 starts

In order of how much they block everything downstream:
1. Beam-width handling (Blocker 1) — this affects which frequencies and
   what `nside_new` are even valid.
2. Prior source (Blocker 2) — Berkhuijsen-based, but the actual
   covariance construction still needs building.
3. `op_alm` parametrisation choice (Stage 2.2) — affects the whole gain
   model, not a detail.
4. Whether to fix or jointly infer β (Stage 4.1).
