# How the TRIS maps were made

**For collaborators.** This records exactly which `limTOD` code, which settings and which
conventions produced `outputs/stage1/tris_maps_tris_haslam_v0.npz` — the 600 and 820 MHz
maps that become two of the four trusted points in the SED fit.

Everything below was read back out of the code and the run manifest, not from memory. The
manifest (`outputs/stage1/stage1_manifest.json`) is the machine-readable version; this
document is the annotated one, and it flags the places where a number is a *decision* we
made rather than something the archive or limTOD told us.

**The short version:** the maps are a prior-regularized MAP reconstruction of two
120-sample drift rings, solved on a 25704-pixel declination band at nside 64, with GSM2008
standing in as the prior. They are **prior-dominated by construction** — 120 samples
carrying roughly 15 independent numbers against 25705 free parameters — and the fitted zero
levels (+1.736 K and +0.515 K) are a measurement of GSM2008's offset, not of the TRIS
archive's zero point. Read the [caveats](#what-to-check-before-you-trust-a-number) before
quoting anything from them.

---

## 1. Provenance

| | |
|---|---|
| **Product** | `outputs/stage1/tris_maps_tris_haslam_v0.npz` |
| **Run name** | `tris_haslam_v0` |
| **Written** | 2026-09-03 19:48:09 UTC |
| **Config** | `configs/run_config.yaml` (the exact block used is copied into the manifest) |
| **Input data** | the four LAMBDA TRIS text products — see `DATA_SOURCES.md` |
| **Archive path at run time** | `../limTOD/downloads/TRIS` |

### Code

| component | version |
|---|---|
| `limTOD` | sibling source checkout, **commit `dbe720b`** ("Merge: make the TRIS map figures say what the numbers already said", 14 Aug 2026) |
| `bayesian_skymap` | `external/bayesian_skymap` (the patched fork), used only to validate `nside_new` |
| Python / numpy / healpy | 3.14.6 / 2.5.2 / 1.20.0 |

`limTOD` was **not** imported as an installed package. The manifest records the path
`/Users/user/Documents/Codes/project3/limTOD/limTOD/tris/__init__.py`, which is the
`paths.limtod` fallback in `configs/run_config.yaml`.

> ⚠️ **The manifest records the path, not the commit.** The commit above was recovered by
> inspecting the checkout afterwards, so it is accurate for this product but is not
> captured automatically. See [known gaps](#known-gaps).

> ⚠️ **That checkout has uncommitted changes**, but none of them are in the map-making
> path. The dirty files are a new untracked `limTOD/tris/binning.py` (the beam-width
> binning map-maker), its export block in `tris/__init__.py`, its tests, and some
> docs/notebooks. Every module stage 1 calls — `archive`, `sky`, `beam`, `geometry`,
> `noise`, `mapmaking` — is clean at `dbe720b`. The maps came from committed code.

---

## 2. What limTOD was asked to do

Five entry points, in call order.

| # | call | module | what it produced |
|---|---|---|---|
| 1 | `read_tris_ring` | `archive` | the 600 / 820 MHz rings, strict-parsed |
| 2 | `read_tris_beam_cuts` | `archive` | the 55-row E/H principal-plane table |
| 3 | `to_tris_temperature_convention` | `sky` | the GSM2008 template put into TRIS's convention |
| 4 | `tris_prior_from_template` | `mapmaking` | per-pixel prior mean and width |
| 5 | `build_tris_mapmaking_inputs` | `mapmaking` | operator, data vector, noise model, pixel selection |

Call 5 is where most of the physics is. Internally it builds the beam with
`tris_cut_beam_map(cuts, nside=64, normalization="peak")`, takes pointing from
`tris_zenith_geometry`, applies `tris_horizon_mask`, selects the band with
`tris_ring_pixels`, and calls `limTOD.simulator.generate_sky2sys_projection` to build the
sky→sample operator. The zero level is appended as a column of ones, making it an explicit
nuisance parameter rather than a rank-1 noise term — algebraically equivalent, but it hands
back the fitted offset.

The solve itself is **not** limTOD's. Stage 1 does it in-repo, deliberately using the same
Krylov machinery as the recalibration step — see [§5](#5-the-solve).

---

## 3. Every setting, and where it came from

`decision` = we chose it, and the reasoning is in `configs/run_config.yaml` next to the key.
`archive` = stated or measured in the TRIS files. `limTOD default` = we did not pass it.

### Data selection

| setting | value | source |
|---|---|---|
| `tris.frequencies_mhz` | `[600, 820]` | **decision** — 2.5 GHz dropped: it is a 6-point set, not a ring, and is flagged poor quality in the literature review |
| `tris.ref_freq_mhz` | `600` | **decision** — better-quantified absolute calibration, so it anchors the spectral model |
| effective frequencies | 600.5 / 817.8 MHz | archive — the *effective* values are used everywhere, never the file labels |
| samples per ring | 120 | archive |

### Geometry

| setting | value | source |
|---|---|---|
| boresight | `RA = LST`, `dec = latitude` | limTOD default (`tris_zenith_geometry`), parked-zenith |
| site latitude | 42.4333° (42°26′) | limTOD constant `TRIS_SITE_LATITUDE_DEG` |
| E-plane roll | `selfrot = −7°` | limTOD constant, from the archive's `NOTE2` ("E-plane tilted 7 degrees Eastwards") |
| declination band | ±45° about the boresight → **−2.57° to +87.43°** | **decision** (`tris.dec_half_width_deg`) |

> Note the boresight sits at **42.4333°**, the site latitude — not the archive's "+42"
> *label*. limTOD keeps those as separate constants and uses the latitude for pointing.

### Beam

| setting | value | source |
|---|---|---|
| beam model | built from the archive's E/H cuts, `normalization="peak"` | limTOD, via `cuts=` |
| measured HPBW | 19.155° (E) / 23.366° (H) | **derived by limTOD from the cut table** — the ring headers say "18 degrees", which is 6% narrow |
| horizon mask | **on** | **decision** (`beam.apply_horizon_mask`) — the cut beam has ~1.2e-4 of its power below the horizon, ≈0.035 K against 300 K ground, comparable to the published 0.066 K systematic |
| `nside_hires` | 64 | **decision** — equal to the working nside; the cut-built beam is already well sampled there (0.92° pixels for a 19° FWHM) |
| beam normalisation in the projection | `normalize_beam=True` | limTOD default via `build_tris_mapmaking_inputs` |

### Pixelisation

| setting | value | source |
|---|---|---|
| `tris.nside` | **64** (49152 pixels; 0.92°/pixel) | **decision** — the grid the map is *solved* on |
| pixels in the band | **25704** (52.3% of sky) | derived from the band cut |
| free parameters | **25705** (25704 sky + 1 zero level) | — |
| ordering / frame | RING, **equatorial** | verified in `notebooks/01_explore_data.ipynb` |

### Noise and the zero level

| setting | value | source |
|---|---|---|
| statistical noise | per-row, from the ring files | archive |
| `tris.uncertainty_floor_k` | `null` → **0.004 K** at both frequencies | **decision** — each ring contains exactly one row with a 0.000 K statistical error, so a floor is mandatory; `null` means "use that ring's smallest strictly positive error" |
| `tris.zero_level_sigma_k` | **0.5 K**, symmetric, both rings | **decision** — see below |
| archive zero levels | 0.066 K at 600 MHz; **+0.430 / −0.300 K** at 820 MHz | archive |

> The 820 MHz zero level is **asymmetric** and limTOD deliberately refuses to symmetrise
> it. 0.5 K is a single Gaussian width wide enough to cover both rings without claiming a
> precision the 820 MHz ring does not have. This is our choice, not the archive's.

### Prior (this is the big caveat)

| setting | value | source |
|---|---|---|
| template | **GSM2008** via `pygdsm`, at the ring's effective frequency | **decision**, and a substitute — see below |
| template pipeline | `ud_grade` → rotate G→C → add CMB monopole back | matches the limTOD walkthrough's order |
| `relative_sigma` | 0.10 | **decision** — "trust the template to 10% of local sky" |
| `absolute_sigma_k` | 0.3 K | **decision** — flat term added in quadrature |
| `floor_sigma_k` | 1e-3 K | **decision** — keeps the prior off a delta function |
| zero-level prior | mean 0, σ = 0.5 K | **decision** |

> ⚠️ **GSM2008 is used as both the prior template and the reference "truth"**
> (`prior_equals_truth: true` in the manifest). There is no true sky for real data and
> Berkhuijsen (1972) is not in hand. Agreement between a stage 1 map and `truth_k` in the
> npz is therefore **circular and is not evidence of anything.**

---

## 4. Temperature convention

The archive states only "Sky Brightness Temperature (K)" — neither the temperature scale
nor whether the CMB monopole is included. limTOD settles it from the data: the two rings
share an RA grid, so the Galactic spectral index between them is a per-sample discriminator,
and only one treatment gives synchrotron (−2.5 to −3.1).

**TRIS temperatures are Rayleigh–Jeans and INCLUDE the CMB monopole.** Everything the
pipeline forward-models is put on that footing by `to_tris_temperature_convention`, which
adds the RJ-equivalent monopole back at the ring's effective frequency. Get this wrong and
every kelvin downstream is wrong.

One assumption nobody can discharge from the files: **the coordinate epoch is unstated.**
B1950→J2000 would shift RA by ~0.6° at this declination — negligible against a 23° beam,
but it is an assumption we own, not a fact.

---

## 5. The solve

Done in-repo (`moonrunes.stage1_tris_maps`), not by limTOD, and deliberately using the same
machinery as the recalibration step so that both fail the same way if they fail at all:
a `scipy.sparse.linalg.LinearOperator` around the MAP normal-equation matvec, plus a
Jacobi-preconditioned Krylov solve.

| setting | value | why |
|---|---|---|
| method | `gmres` | **decision** |
| `rtol` | **1e-12** | not the usual 1e-9: the prior width spans an order of magnitude across the band, so a 1e-9 residual still leaves the map ~1e-3 σ off |
| `maxiter` | 1000 | |
| restart | `null` → **160** | `n_samples + 40`; the data term has rank 120, so a restart longer than that terminates inside one cycle instead of stalling at the default 20 |
| Woodbury fast path | **never** | it is documented wrong by ~8 orders of magnitude at low channel count, and TRIS gives 2 |

**Because "it finished" is not "it converged"**, every solve is cross-checked against the
exact Woodbury identity at a single frequency — a 120×120 inner matrix, cheap, and *not*
the multi-frequency approximation above. The same identity supplies the posterior σ, which
a Krylov solve alone cannot give. Stage 1 refuses to write a product that fails the check.

| | 600 MHz | 820 MHz |
|---|---|---|
| gmres info / iterations | 0 / 114 | 0 / 72 |
| relative residual | 7.9e-13 | 3.1e-14 |
| agreement with exact solve | 1.6e-07 σ | 1.8e-09 σ |
| tolerance | 1e-04 σ | 1e-04 σ |

---

## 6. What came out

| | 600 MHz (νeff 600.5) | 820 MHz (νeff 817.8) |
|---|---|---|
| fitted zero level | **+1.7361 ± 0.0108 K** | **+0.5148 ± 0.0069 K** |
| reduced χ² (dof 120) | 33.62 | 13.25 |
| residual RMS | 0.083 K | 0.075 K |
| prior → posterior σ (median) | 1.1553 → 1.1550 K | 0.7044 → 0.7042 K |
| pixels whose σ shrank >5% | **2** of 25704 | **2** of 25704 |
| minimum beam coverage | 0.9993 | 0.9993 |
| monopole degeneracy | ≈1.0 | ≈1.0 |

Those fitted offsets reproduce the GSM2008 deficits measured independently in
`notebooks/01_explore_data.ipynb` (1.92 K and 0.56 K cold) — a real check that the
geometry, beam and convention translation came through intact.

---

## What to check before you trust a number

1. **The maps are prior-dominated.** The beam suppresses harmonics above *m* ≈ 8, so 120
   samples carry roughly **15 independent numbers** — against 25705 free parameters. The
   σ-shrinkage row above is the honest measure: the posterior σ is within 0.002% of the
   prior σ for all but **2 pixels**. Structure you see between ring samples is the prior
   interpolating, not TRIS measuring.
2. **The zero levels measure GSM2008, not the archive.** The monopole degeneracy is ≈1, so
   with a per-pixel prior this tight the nuisance offset is the only place a template error
   can go.
3. **The reduced χ² of 33.6 is expected, not a failure.** It is the signature of a template
   that is kelvins cold, forward-modelled against absolutely calibrated data.
4. **Nothing here is Haslam-independent yet** — GSM2008 is built from surveys that include
   Haslam 408 MHz.
5. **2.5 GHz is not in these maps** at all.

## Known gaps

- **The manifest records limTOD's path, not its commit.** A product cannot currently be
  traced to the code that made it, which is the one thing the manifest convention exists to
  guarantee. Recording `git rev-parse HEAD` plus a dirty flag would close this — and given
  the checkout *is* dirty, it is not a hypothetical.
- **`beam.selfrot_deg: -7.0` in the config is never read.** Stage 1 does not pass
  `geometry`, so limTOD's own `TRIS_E_PLANE_EAST_OF_MERIDIAN_DEG` is what acts. The two
  agree today and nothing checks that they do.
- **`beam.e_plane_hpbw_deg` / `h_plane_hpbw_deg` are only checked against each other**, to
  confirm `beam.beam_deg` is their mean. Neither is ever compared against
  `cuts.half_power_full_width_deg()`, which is what limTOD actually measures from the file.

## Reproducing this

```bash
git submodule update --init                      # external/bayesian_skymap
python run_pipeline.py --stage 1
```

Requires the TRIS archive at `paths.tris_archive_dir`, `pygdsm` for the GSM2008 template,
and `limTOD` (installed, or the sibling checkout at `paths.limtod`). Stage 1 refuses to
overwrite an existing product unless `run.overwrite` is set. It takes about 12 s, of which
the two Krylov solves are 0.45 s.

Read the products back through the one contract everything downstream uses:

```python
from moonrunes.stage1_tris_maps import load_stage1
products = load_stage1()
products["sky_full_k"]               # (Nfreq, 49152), NaN outside the band
products["band_mask"]                # the reliable mask -- use this, not the NaNs
products.diagnostics(600)["solver"]  # info, iterations, agreement with the exact solve
```
