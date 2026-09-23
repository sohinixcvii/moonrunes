# How the TRIS maps were made

**For collaborators.** This records exactly which `limTOD` code, which settings and which
conventions produced `outputs/stage1/tris_maps_tris_haslam_v0.npz` — the 600 and 820 MHz
maps that become two of the four trusted points in the SED fit.

Everything below was read back out of the code and the run manifest, not from memory. The
manifest (`outputs/stage1/stage1_manifest.json`) is the machine-readable version; this
document is the annotated one, and it flags the places where a number is a *decision* we
made rather than something the archive or limTOD told us.

**The short version:** the maps are a prior-regularized MAP reconstruction of two
120-sample drift rings, solved by **limTOD's own solver** on a deliberately coarse
**nside 16** grid, over the **1272 pixels the beam actually saw**, with GSM2008 standing in
as the prior. They remain **prior-dominated** — 120 samples carrying roughly 15 independent
numbers — and the fitted zero levels (+1.744 K and +0.510 K) are a measurement of GSM2008's
offset, not of the TRIS archive's zero point. Read the
[caveats](#what-to-check-before-you-trust-a-number) before quoting anything from them.

---

## 1. Provenance

| | |
|---|---|
| **Product** | `outputs/stage1/tris_maps_tris_haslam_v0.npz` |
| **Run name** | `tris_haslam_v0` |
| **Written** | 2026-09-21 16:37 UTC |
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
| `tris.nside` | **16** (3072 pixels; 3.66°/pixel) | **decision** — this is *supposed* to be a coarse map: the beam is 19–23°, so nside 16 is already ~6× finer than the resolution behind it, and stage 3 degrades to `beam.nside_new = 8` regardless |
| `tris.beam_response_threshold` | **0.01** | **decision** — a pixel is retained if the beam response exceeds 1% of the peak for **at least one observation** |
| pixels retained | **1272** (41.4% of sky) | measured, not geometric — see below |
| beam power retained | **99.604%** | what the 1% cut costs |
| free parameters | **1273** (1272 sky + 1 zero level) | — |
| ordering / frame | RING, **equatorial** | verified in `notebooks/01_explore_data.ipynb` |

### Pixel selection — on the beam, not on a declination band

A pixel is retained when the beam response at it exceeds `tris.beam_response_threshold` of
the peak response **for at least one observation**. The comparison is against that
observation's own peak, not a global maximum.

This replaced the earlier `dec_half_width_deg: 45.0` band cut, which selected on geometry
alone: it kept pixels inside the band the beam never illuminated, and cut pixels in the
beam's wings that it did. The criterion now matches the physical statement — keep what the
horn actually saw.

**Convergence over the threshold**, as the collaborator asked, at 0.003 / 0.01 / 0.03:

| threshold | pixels | sky fraction | beam power kept | zero level 600 | zero level 820 | reduced χ² 600 / 820 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.003 | 1469 | 47.8% | 99.865% | +1.7127 K | +0.4921 K | 3.10 / 3.21 |
| **0.01** | **1272** | **41.4%** | **99.604%** | **+1.7442 K** | **+0.5101 K** | **3.09 / 3.22** |
| 0.03 | 1099 | 35.8% | 98.919% | +1.8245 K | +0.5569 K | 3.12 / 3.21 |

**The map is converged; the zero level is not.** On the 1099 pixels common to all three
selections, the temperatures agree to a maximum of **0.40 σ** (0.003 vs 0.03 against the
0.01 posterior width), median |ΔT| of 0.013–0.027 K — immaterial, as expected.

The fitted zero level is another matter: it moves **+0.112 K at 600 MHz** and **+0.065 K at
820 MHz** across the decade of threshold, monotonically increasing as the cut tightens.
That is **2.6 σ and 2.4 σ** of its own fitted uncertainty, though only 22% and 13% of the
0.5 K prior width.

The direction makes sense — a looser cut admits more wing pixels, the sky model absorbs
more of the signal, and less is left for the offset to carry. Treat the fitted zero level
as carrying a **systematic of order 0.1 K from the pixel-selection choice**, on top of its
quoted statistical error. It does not affect the maps.

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

**Current code — limTOD's solver.** `TRISMapMakingInputs.solve` wraps
`limTOD.wiener_filter_map`, which forms the dense normal matrix
`AᵀN⁻¹A + S⁻¹ + εI` and factorises it directly. limTOD builds the prior itself from
`prior_map` / `prior_sigma_k` and takes the zero-level column's width from the
`zero_level_sigma_k` given to `build_tris_mapmaking_inputs`.

| setting | value |
|---|---|
| `tris.solver.regularization` | `1e-12`, added to the diagonal before inversion — a bias on the answer as well as a numerical guard, so it stays as small as conditioning allows |
| cross-check | unchanged: every solve is compared against the exact Woodbury identity and the product is refused if it misses by more than `cross_check_max_sigma` |

**Measured cost.** The dense path scales as the cube of the parameter count. At the
configured nside 16 it is effectively free:

| `tris.nside` | parameters | solve time | peak memory |
|---|---:|---:|---:|
| **16 (configured)** | **1,273** | **0.06 s** | trivial |
| 32 | 6,461 | 4.82 s | ~1.6 GB |
| 64 | 25,705 | — | ~25 GB, did not complete on 16 GB |

nside 64 is recorded here only to show where the ceiling is: `wiener_filter_map` allocates
`S_inv`, `AtNA`, `regularization * np.eye(n)`, `covariance_inv` and the factorisation, each
n² float64. At the coarse resolution this map is supposed to have, none of that binds.

**The previous in-repo solver** (`krylov_map_solve`) is still present and still pinned by
the test suite against a dense reference, but it no longer makes the maps.

**Because "it finished" is not "it converged"**, every solve is cross-checked against the
exact Woodbury identity at a single frequency — a 120×120 inner matrix, cheap, and *not*
the multi-frequency approximation above. The same identity supplies the posterior σ, which
a Krylov solve alone cannot give. Stage 1 refuses to write a product that fails the check.

| | 600 MHz | 820 MHz |
|---|---:|---:|
| agreement with exact Woodbury | 3.96e-10 σ | 1.07e-10 σ |
| tolerance | 1e-04 σ | 1e-04 σ |
| solve time | 0.06 s | 0.06 s |

---

## 6. What came out

| | 600 MHz (νeff 600.5) | 820 MHz (νeff 817.8) |
|---|---:|---:|
| fitted zero level | +1.7442 ± 0.0425 K | +0.5101 ± 0.0269 K |
| reduced χ² (dof 120) | 3.09 | 3.22 |
| residual RMS | 0.025 K | 0.036 K |
| prior → posterior σ (median) | 1.1846 → 1.1559 K | 0.7162 → 0.7053 K |
| pixels whose σ shrank >5% | **388 of 1272** | **201 of 1272** |
| minimum beam coverage | 0.9958 | 0.9958 |
| monopole degeneracy | ≈1.0 | ≈1.0 |

Those fitted offsets reproduce the GSM2008 deficits measured independently in
`notebooks/01_explore_data.ipynb` (1.92 K and 0.56 K cold) — a real check that the
geometry, beam and convention translation came through intact.

---

## What to check before you trust a number

1. **The maps are still prior-dominated, but much less so.** The beam suppresses harmonics
   above *m* ≈ 8, so 120 samples carry roughly **15 independent numbers** — against 1,273
   free parameters. At the old nside 64 the posterior σ beat the prior by >5% for just
   **2** of 25,704 pixels; at nside 16 it does so for **388 of 1272** at 600 MHz and **201 of
   1272** at 820 MHz. Coarsening the grid is what bought that.
2. **The zero levels measure GSM2008, not the archive.** The monopole degeneracy is ≈1, so
   with a per-pixel prior this tight the nuisance offset is the only place a template error
   can go.
3. **The reduced χ² of ~3.1 is expected, not a failure.** It is the signature of a template
   that is kelvins cold, forward-modelled against absolutely calibrated data.
4. **The fitted zero level carries a ~0.1 K systematic from the pixel-selection
   threshold**, on top of its quoted statistical error. See the convergence table.
5. **Nothing here is Haslam-independent yet** — GSM2008 is built from surveys that include
   Haslam 408 MHz.
6. **2.5 GHz is not in these maps** at all.

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
