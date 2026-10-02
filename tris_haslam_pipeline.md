# TRIS → Haslam Calibration: Literal Pipeline
### Pixel-by-pixel SED fitting using TRIS + ARCADE2 to calibrate Haslam
### (supersedes the earlier `bayesian_skymap` Gibbs-sampler design)

**Primary aim:** make maps from TRIS TOD data, use them (plus ARCADE2) to
calibrate the Haslam map, via a per-pixel SED fit rather than a joint
Gibbs sampler.

**Why this replaced the original plan:** the `bayesian_skymap` Gibbs
approach hit three real blockers (beam-width hardcoding, no real-data
prior, a documented Woodbury solver bug at low channel count). A
pixel-by-pixel least-squares SED fit sidesteps all three at once — no
truth-informed prior needed, no Woodbury solver, and the beam mismatch
becomes a smoothing/regridding problem rather than a code-patching one.
The old blockers are kept below, marked resolved-by-redesign, so the
reasoning isn't lost.

---

## The core idea

At each low-resolution sky pixel, assemble a small spectral energy
distribution (SED) — brightness temperature vs. frequency — from every
dataset available at that pixel:

| Frequency | Source | Trusted (absolutely calibrated)? |
|---|---|---|
| 408 MHz | Haslam | **No** — this is what's being calibrated |
| 600 MHz | TRIS | Yes |
| 820 MHz | TRIS | Yes |
| 3150 MHz | ARCADE2 | Yes (real spatial map, confirmed) |
| 3410 MHz | ARCADE2 | Yes (real spatial map, confirmed) |

Fit a power law $T(\nu) = A(\nu/\nu_0)^\beta$ to the **4 trusted points
only**, evaluate the fit at 408 MHz to get a predicted "true" Haslam
value at that pixel, and compare against Haslam's actual raw value to
derive a per-pixel gain/offset correction.

**Why ARCADE2's bands matter beyond "more data":** with only TRIS's 2
bands, a 2-parameter power law fit to 2 points is exactly determined —
zero residual, no way to check the power-law assumption itself. With
ARCADE2's 2 bands added, the fit is over-determined (4 points, 2
parameters), giving a genuine per-pixel chi-square goodness-of-fit check.

**Physical caveat to keep in view:** 408 MHz–3.41 GHz spans nearly a
decade; free-free emission's relative contribution rises with frequency,
which can introduce real curvature a pure power law won't capture. If
per-pixel chi-square is poor across many pixels, that's evidence the
single power-law model is too simple, not a fitting failure to force
through.

---

## Coordinate system mismatch — resolved (reprojection to Equatorial)

Haslam ships in **Galactic** coordinates; TRIS work to date has been in
**Equatorial** (RA/Dec, LST-based). ARCADE2 is **Galactic** too — its
headers say so as `SKYCOORD='Galactic'`, a keyword healpy does not read.
The TRIS stage 1 products declare no frame; they are equatorial by
construction, and the rings peak at the Cygnus plane crossing, as they
should. `src/moonrunes/frames.py` rotates Haslam (harmonic space, full
sky) and ARCADE2 (pixel space, mask-aware) to Equatorial; the checks —
landmarks at their literature positions, the ARCADE2 exact-zero mask
surviving the rotation — are in `notebooks/04_coordinate_frames.ipynb`.
Not yet wired into stage 2's beam matching.

---

## Old blockers — resolved by the architecture change, not deleted

### Former Blocker 1 — beam width mismatch -> now a smoothing/regridding problem
No longer about patching `bayesian_skymap`'s hardcoded `beam_deg`. See
**Step 0** below for the actual beam-matching procedure.

### Former Blocker 2 — no "truth" for a real-data prior -> moot
The SED fit needs no prior at all — it's anchored directly by TRIS/ARCADE2
data at each pixel. The Berkhuijsen-prior research is no longer a blocking
dependency, though it remains a good independent validation cross-check.

### Former Blocker 3 — Woodbury solver bug at low channel count -> moot
No Gibbs sampler, no Woodbury solver, in this design at all.

---

## Step 0 — Beam-matching (new, real requirement)

Three datasets, three native resolutions: Haslam (~1 deg), ARCADE2 (a few
degrees — confirm exact native FWHM before computing the numbers below),
TRIS (~19-23 deg, asymmetric). **TRIS is the coarsest — smooth the other
two down to match it, never the reverse.**

1. **Pick a common target FWHM.** TRIS's beam is asymmetric (19.155 deg E /
   23.366 deg H); `healpy`'s smoothing only takes one symmetric FWHM. Use
   the **larger axis (23.366 deg)** as the common target — conservative,
   so nothing ends up under-smoothed in either direction. (A proper
   elliptical-beam convolution is more correct but substantially more
   work; flag this as a deliberate simplification, not an oversight.)

2. **Compute the *additional* smoothing needed, not the raw target FWHM.**
   Gaussian beams add in quadrature: if a map already has native
   resolution FWHM_native and the target is FWHM_target, the extra kernel
   needed is
   ```
   FWHM_extra = sqrt(FWHM_target^2 - FWHM_native^2)
   ```
   Smoothing directly to the target FWHM without this correction
   over-smooths.

3. **Handle ARCADE2's mask *before* smoothing, not after.** This is the
   most important practical gotcha given the masking issue just found and
   fixed: `hp.smoothing` doesn't know about `UNSEEN` — the spherical
   harmonic transform is global, so naively smoothing a masked map bleeds
   incorrect values across the mask boundary into the ring itself.
   Standard fix:
   - Set `UNSEEN` pixels to `0` (not left as `UNSEEN`) before smoothing.
   - Separately smooth a **binary weight map** (1 where observed, 0 where
     masked) with the *same* kernel.
   - Divide the smoothed data map by the smoothed weight map to correct
     for the boundary leakage.
   - Re-mask any pixel where the smoothed weight falls below some
     threshold (e.g. 0.5) — these are unreliable edge pixels, not
     recoverable by the division above.

4. **Smooth before regridding, never after.** Smoothing on an
   already-downgraded (coarse) map risks aliasing. Order: smooth each map
   at a sufficiently fine native `nside` first, *then* `ud_grade` down to
   the common coarse grid.

5. **Regrid all three to TRIS's own map-making pixelization** (nside=16
   since 2026-10-02, ~3.66 deg pixels; originally nside=8, ~7.33 deg pixels, per your walkthrough's "~3 pixels per beam FWHM"
   convention) — so Step 1's per-pixel SEDs are being assembled on a grid
   that actually matches the physical resolution of the trusted
   calibrators, not an arbitrarily finer one.

```python
import healpy as hp
import numpy as np

def smooth_to_target(m, native_fwhm_deg, target_fwhm_deg, weight_mask=None):
    """
    Smooth a map from its native resolution to a coarser target FWHM,
    correctly handling a binary observed/unobserved mask if given.
    """
    extra_fwhm_rad = np.radians(
        np.sqrt(target_fwhm_deg**2 - native_fwhm_deg**2)
    )

    if weight_mask is None:
        return hp.smoothing(m, fwhm=extra_fwhm_rad)

    m_zeroed = np.where(weight_mask, m, 0.0)
    m_smooth = hp.smoothing(m_zeroed, fwhm=extra_fwhm_rad)
    w_smooth = hp.smoothing(weight_mask.astype(float), fwhm=extra_fwhm_rad)

    with np.errstate(invalid='ignore', divide='ignore'):
        corrected = m_smooth / w_smooth
    corrected[w_smooth < 0.5] = hp.UNSEEN
    return corrected
```

## Stage 1 — Produce calibrated TRIS maps (`limTOD.tris`)

Unchanged from the original plan:
1. Read raw TRIS archive data for 600 MHz and 820 MHz (2.5 GHz dropped —
   known quality issue).
2. Build the beam from the archive's own E/H cuts (HPBW 19.155/23.366 deg,
   `selfrot=-7 deg`).
3. Apply the verified zenith/roll geometry (`RA=LST`, `dec=lat=42d26m`).
4. Run map-making (the simplified beam-width scheme, or the
   prior-regularized map-maker if the Berkhuijsen-based prior is wanted
   as an independent cross-check even though it's no longer required).
5. Output: calibrated TRIS maps at each frequency.

## Stage 2 — Beam-match and regrid (Step 0's procedure, applied)

Apply Step 0 to Haslam and both ARCADE2 bands, regridding everything
(including TRIS's own maps) onto the common nside=16 grid, in a common
coordinate system (resolves the Open TODO above — do this *before*
Stage 3, not after).

## Stage 3 — Assemble per-pixel SEDs and fit

1. At each common-grid pixel, build the 5-point SED (408, 600, 820, 3150,
   3410 MHz).
2. Fit the power law to the 4 trusted points (600-3410 MHz).
3. Evaluate the fit at 408 MHz -> predicted true Haslam value.
4. Record the per-pixel chi-square of the fit — this is the diagnostic
   for where the power-law assumption breaks down.

## Stage 4 — Derive and apply the calibration correction

1. Decide the error model explicitly: multiplicative gain g, additive
   offset c, or both jointly (the literature — your own TRIS/Sky Map
   Priors research — suggests Haslam has *both* a scale uncertainty
   (~3%) and a real zero-level offset, so a joint fit is likely the more
   physically honest choice).
2. Apply the resulting per-pixel correction back onto full-resolution
   Haslam — uniformly within each low-res patch as a first pass; smoothing
   the correction across patches is a natural refinement to avoid blocky
   artifacts at patch boundaries.

## Stage 5 — Validate

1. Check the per-pixel chi-square distribution from Stage 3 — flags where
   the power-law SED assumption is failing, not just whether the fit ran.
2. Compare recovered g/c against the literature benchmarks already in
   hand (~3% scale, ~0.91 K zero-level) — same sanity check as before,
   still valid.
3. Sanity-check the fitted spectral index beta against the expected
   Galactic synchrotron range (~-2.5 to -3.0).

---

## What to actually decide, in order

1. **Coordinate system reprojection** (Open TODO) — blocks everything
   downstream if wrong.
2. **Common target FWHM and mask-handling in Step 0** — the beam-matching
   choices propagate into every later pixel value.
3. **Error model** (Stage 4.1) — multiplicative, additive, or both.
4. Confirm ARCADE2's actual native beam FWHM before computing Step 0's
   quadrature-difference smoothing kernel precisely.