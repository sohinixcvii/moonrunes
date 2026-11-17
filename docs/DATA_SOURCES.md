# Data sources

Every external dataset this pipeline reads, where it comes from, and what consumes it.
**Nothing in this repo touches the network** — each stage reads a local file (or an
astropy cache populated by `pygdsm`) and raises with a clear message when it is absent.
So this file is the download list.

Every product below except the last section is hosted by **LAMBDA** (NASA/GSFC's Legacy Archive for
Microwave Background Data Analysis).

---

## TRIS — the absolutely calibrated low-resolution data

<https://lambda.gsfc.nasa.gov/product/tris/tris_prod_table.html>

TRIS is an absolute radiometer. Stage 1 makes maps from the published drift rings.

| | |
|---|---|
| **Consumed by** | Stage 1 (`moonrunes.stage1_tris_maps`), via `limTOD.tris`'s strict readers |
| **Config** | `paths.tris_archive_dir` (default `../limTOD/downloads/TRIS`) |
| **Files** | `TRIS_absolute_600.txt`, `TRIS_absolute_820.txt`, `TRIS_absolute_2500MHz.txt`, `TRIS_Beam_Profile.txt` |
| **Format** | Plain text; the rings carry temperature with separate statistical and systematic errors, and the beam file the E/H principal-plane cuts |

Notes that matter downstream:

* **2.5 GHz is downloaded but not used.** `tris.frequencies_mhz` is `[600, 820]` — the
  2.5 GHz ring is flagged poor quality in the TRIS literature review. Keep the file
  anyway; the reader is strict and the notebooks look at it.
* The beam file is where the **19.155° (E) / 23.366° (H)** HPBW comes from, and hence
  `beam.beam_deg = 21.2605`.
* Each published ring contains exactly one row with a **0.000 K statistical error**, which
  is why `tris.uncertainty_floor_k` exists at all.
* The archive's zero-level uncertainties are **asymmetric at 820 MHz**
  (+0.430 / −0.300 K, against 0.066 K at 600 MHz); `limTOD` refuses to symmetrise it,
  which is why `tris.zero_level_sigma_k` is a single deliberately loose 0.5 K.

---

## Haslam 408 MHz — the high-resolution map being recalibrated

<https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_get.html>

The high-resolution map whose scale and zero level the pipeline tests against TRIS.

| | |
|---|---|
| **Consumed by** | Stage 2 (`moonrunes.stage2_beam_matching`), and `diagnostics.ipynb` Parts A–B (were notebooks 01 and 04; 05–06 removed 2026-11-16) |
| **Config** | `paths.haslam_map_fits` = `res/HASLAM.fits`. The superseded `stage2_haslam_prep` read `paths.haslam_map`, where `null` meant "use the copy `pygdsm` already caches" |
| **File** | `haslam408_dsds_Remazeilles2014.fits` — the **destriped, desourced reprocessed** map (Remazeilles et al. 2015), not the original 1982 release |
| **Format** | HEALPix FITS, `COORDSYS = GALACTIC`, RING, nside 512, `BEAMSIZE` = 56′, K (RJ, CMB included) |

Notes:

* `null` in `paths.haslam_map` is **not** an open decision (unlike other nulls in the
  config): it means stage 2 pulls the file from the astropy cache that `pygdsm` populates,
  which is the same reprocessed map. Set the path to override with a local download.
* The map is total-intensity, so stage 2 **removes the CMB monopole** (2.7157 K, the
  Rayleigh–Jeans value at 408 MHz) before anything else.
* This is also the source of the circularity warning throughout the repo: GSM2008 is
  built from eleven surveys, **one of which is this map**, so the default prior is not
  Haslam-independent. See the prior table in `README.md`.

---

## ARCADE 2 — the two extra SED points

<https://lambda.gsfc.nasa.gov/product/arcade/arcade_maps_get.html>

ARCADE 2 is a balloon-borne, absolutely calibrated radiometer. Under the pixel-by-pixel
SED design in `tris_haslam_pipeline.md` its two lowest bands are trusted points 3 and 4:
TRIS alone gives two, and a two-parameter power law through two points is exactly
determined, so ARCADE 2 is what makes the per-pixel fit over-determined and its χ²
meaningful.

| | |
|---|---|
| **Consumed by** | Stage 2 (`moonrunes.stage2_beam_matching`, `paths.arcade2_map_3150` / `_3410`), and the SED fit prototype in `notebooks/06` (removed 2026-11-16) as two of the four trusted points |
| **Looked at in** | `notebooks/diagnostics.ipynb` Part A (was `01_explore_data`), sections 02 and 05, and Part B (was `04`); `05`–`06` archived |
| **Files** | `res/ARCADE2_315.fits`, `res/ARCADE2_341.fits` — 3.15 and 3.41 GHz |
| **Format** | HEALPix FITS, **`NESTED` at nside 16** (a 3.7° pixel), **K thermodynamic with the CMB monopole included** (only the dipole removed — LAMBDA `arcade_maps_info`), `BEAMSZ = 11.6°` in HDU 0 |

Three properties measured in the notebook, each of which changes how the maps must be
handled:

* **unobserved sky is stored as exact zero** — no `hp.UNSEEN`, no `NaN`, and ~93% of
  pixels (2852/3072 at 3.15 GHz, 2840/3072 at 3.41 GHz) are exactly 0. The mask has to be
  rebuilt from the zeros, and rebuilt *before* smoothing, because `hp.smoothing` works
  through a global spherical harmonic transform and would bleed those zeros across the
  boundary. Step 0.3 of the pipeline document is the procedure.
* **`NESTED` ordering** — `hp.read_map` reorders to `RING` on read, silently.
* **Galactic, stated as `SKYCOORD='Galactic'`**, in both HDUs — not as `COORDSYS`,
  which is the keyword healpy reads, so healpy reports no frame. Confirmed against the
  sky as well as the header: the brightest pixels sit on b ≈ 0 in the file's own frame,
  and after rotation to equatorial the maps correlate at 0.95 with rotated Haslam
  (`notebooks/diagnostics.ipynb` Part B, was `04_coordinate_frames`).

One caveat that is design, not format: 3.15/3.41 GHz against Haslam's 408 MHz is close to
a decade of lever arm, over which free-free emission's rising contribution can put real
curvature in the SED that a single power law will not capture.

---

## Maipu + MU 45 MHz, and Stockert + Villa-Elisa 1420 MHz — reference maps

<https://lambda.gsfc.nasa.gov/product/foreground/fg_diffuse.html>

Two more all-sky surveys in `res/`. **No stage reads them.** They are shown alongside the
others in `notebooks/diagnostics.ipynb` Part A (was `01_explore_data`), sections 05–06. Neither header gives a
frequency, unit or beam, so those come from the LAMBDA product pages.

| | Maipu + MU 45 MHz | Stockert + Villa-Elisa 1420 MHz |
|---|---|---|
| **File** | `res/MAIPU_45mhZ.fits` | `res/STOCKERT_VILLA_1400MHZ.fits` (the survey is **1420** MHz, not the 1400 in the name) |
| **Format** | HEALPix, `COORDSYS = G`, NESTED, nside 64 | HEALPix, `COORDSYS = G`, NESTED, nside 256 |
| **Units** | K | **mK** (full-beam brightness temperature) |
| **Beam** | 4.6° × 2.4° (Maipu, south) / 3.6° (MU, north) | 35.4′ |
| **Unobserved** | −32768 (FITS null), 1877 px: caps around both celestial poles, data only Dec −79° to +71° | −32768, 47 px |
| **Reference** | Guzmán et al. 2011, A&A 525, A138 (LAMBDA `fg_maipu_info`) | Reich 1982; Reich & Reich 1986; Testori et al. 2001 (LAMBDA `fg_stockert_villa_info`) |

## CHIPASS 1400 MHz — not in the repository

<https://lambda.gsfc.nasa.gov/product/foreground/fg_chipass_info.html>

At 104.9 MB the file is over GitHub's 100 MB limit, so `res/CHIPASS_1400MHz.fits` is
git-ignored. Download it from the LAMBDA page above and save it under that name. It is a
FITS image in the HPX projection, in mK, with a 14.4′ beam and the CMB included
(Calabretta et al. 2014). `moonrunes.sky_maps` converts it back to HEALPix.

## Uncertainties (noise model) for every map

*Updated 2026-10-28.* No map file in `res/` carries an error, variance or standard-deviation column (only
temperature, plus `N_OBS` for ARCADE 2), so every map's uncertainty comes from its own publication; TRIS's come from
the archive files. This section is the single reference for what each map is given and where.

### The error model

For a survey map, the 1σ uncertainty of one point (a ring sample or a stripe pixel) with CMB-free Rayleigh–Jeans
temperature T is

    σ = sqrt( (scale · |T|)² + zero_level² + noise² )

(`sed_prototype.sigma_for`). The three terms, all 1σ:

* **scale** — the temperature-scale (gain) error, a fraction of T;
* **zero level** — the error of the map's absolute offset, in K;
* **noise** — the rms receiver noise, in K.

A term the publication does not quote is 0. The values are in `moonrunes.sed_prototype.PUBLISHED_ERRORS` (one entry
per file, with its source and any judgement call) and, for the model comparison, in `configs/error_bars.csv`
(kind `published`, the same numbers). ARCADE 2 is the exception: a constant calibration error
(`sed_prototype.ARCADE_CAL_K`). A map with no entry falls back to 10% of T (`ASSUMED_FRACTION`); no
fitted map uses it (until 2026-10-21 every survey did).

For TRIS (`sed_prototype.ring_sigma`), per ring sample:

    σ = sqrt( statistical² + zero_level² )

* **statistical** — the archive's per-sample column. Each ring has exactly one row with 0.000 K, so values of 0 are
  replaced by the ring's smallest positive value (0.004 K at both frequencies);
* **zero level** — the archive's per-ring value, 0.066 K at 600 MHz and +0.43/−0.30 K at 820 MHz, of which the larger
  side (0.43 K) is used.

### Per map

| Map (file) | ν [MHz] | Scale | Zero level | Noise | Typical σ on the TRIS ring | Source, and judgement calls | Used in |
|---|---|---|---|---|---|---|---|
| LWA1 35 MHz (`LWA_35MHz.fits`) | 35 | 20% | — | — | 20% (≈ 2 840 K of 14 100 K) | LAMBDA `fg_lwa1_radio_maps_info` (Dowell et al. 2017): "a conservative estimate of the calibration accuracy … 15–25 per cent", for structure and sky average. **20% is the middle of that range.** No zero level or noise is quoted. The paper's per-pixel uncertainty maps are not in `res/`. | SED fits, model comparison, T–T fits |
| Maipu + MU 45 MHz (`MAIPU_45mhZ.fits`) | 45 | 1.7% | 39 K | — | 1.8% (≈ 134 K of 7 500 K) | Guzmán et al. 2011 quote no error. **The 1σ halves of the corrections Monsalve et al. 2021 find against EDGES** (×1.076 ± 0.034 and −160 ± 78 K, 2σ; LAMBDA `fg_maipu_info`). **The corrections themselves are not applied** (TODO N-16). | SED fits, model comparison, T–T fits |
| LWA1 80 MHz (`LWA_80MHz.fits`) | 80 | 20% | — | — | 20% (≈ 428 K of 2 140 K) | As LWA1 35 MHz. | SED fits, model comparison, T–T fits |
| Haslam 408 MHz, destriped (`HASLAM_DESTRIPED_ONLY.fits`) and original (`HASLAM.fits`) | 408 | 10% | 3 K | — | 14% (≈ 4.2 K of 29 K) | Haslam et al. 1982. Never a fit point (it is the map being tested). | T–T fits (`notebooks/tests_analysis/` Part E, was notebook 14; `10_paper_plots.ipynb`), Haslam's error band in plots, the noise test of `notebooks/tests_analysis/` Part K2 (was notebook 03) |
| TRIS 600 MHz (ring data) | 600.5 | — | 0.066 K | per sample (0.004–0.020 K, median 0.010 K) | 0.5% (0.066–0.069 K, median 0.067 K, of 12.4 K) | The archive file `TRIS_absolute_600.txt` (Zannoni et al. 2008). The zero-level term dominates. | everything with TRIS ring data |
| TRIS 820 MHz (ring data) | 817.8 | — | +0.43/−0.30 K (0.43 K used) | per sample (0.004–0.041 K, median 0.011 K) | 8.5% (0.430–0.432 K, median 0.430 K, of 5.1 K) | The archive file `TRIS_absolute_820.txt`. **The larger side of the asymmetric zero level is used.** | everything with TRIS ring data |
| CHIPASS 1400 MHz (`CHIPASS_1400MHz.fits`) | 1400 | 5% | 0.5 K | 30 mK | — (does not reach the ring) | Calabretta et al. 2014 §4: "a typical rms error of 30 mK"; zero level tied to the 1420 MHz map. **It quotes no scale or zero-level error, so Stockert's are used.** | map route only (the southern edge of the stripe; never ≥ half a beam on the ring) |
| Stockert + Villa-Elisa 1420 MHz (`STOCKERT_VILLA_1400MHZ.fits`) | 1420 | 5% | 0.5 K | 17 mK | 47% (≈ 0.50 K of 1.07 K) | LAMBDA `fg_stockert_villa_info`; Testori et al. 2001: gain scale ≈ 5%, zero level 0.5 K, sensitivity 50 mK = 3 × rms. **Noise = 50 mK / 3.** | SED fits, model comparison, T–T fits |
| Rhodes/HartRAO 2326 MHz (`RHODES_HARTRAO_2326MHz.fits`) | 2326 | 5% | 80 mK | 30 mK | — (does not reach the ring) | Jonas, Baart & Nicolson 1998: scale < 5%, zero level < 80 mK, rms < 30 mK. **The upper limits are used as 1σ.** | map route only (the southern edge of the stripe) |
| ARCADE 2 (`ARCADE2_315.fits`, `ARCADE2_341.fits`) | 3150, 3410 | — | — | 7.5, 6.0 mK (constant) | — | Fixsen et al. 2011 Table 3, used as a constant σ. | nothing: excluded by `sed.exclude_maps` |
| Juno 600, 1248 MHz (`Juno_600MHz.fits`, `Juno_1248MHz.fits`) | 600, 1248 | — | — | — | — | No entry (unit, frame and beam unknown, so not usable); would fall back to 10%. | nothing: unusable, and in `sed.exclude_maps` |
| WMAP K, Ka (`WMAP_K.fits`, `WMAP_Ka.fits`) | — | — | — | — | — | No entry: unusable (frequency contradicted, monopole removed). | nothing |

"Typical σ on the TRIS ring" is the median over the 120 beam-matched ring samples (notebook 10's input table,
2026-10-28). The percentages for the scale-only surveys are exactly their scale term.

TRIS **maps** (the map route only) use a different error: the stage 1 Wiener map's posterior σ per pixel, plus the
fitted zero level's σ in quadrature when the zero level is added (`configs/model_comp_config.yaml`
`zero_level.add_zhat_to_tris`).

### Where the errors enter

* **SED fits** (`sed_prototype.tod_route` / `map_route`; notebooks 07, 08, both removed 2026-11-16): `sigma_for` for the surveys,
  `ring_sigma` (or the stage 1 σ) for TRIS.
* **Model comparison** (`scripts/model_comparison/run_comparison.py`; `notebooks/tests_analysis/` Parts C, D, were notebooks 11, 15): `error_bars.csv` for the surveys (the
  same numbers), `ring_sigma` for TRIS, in a Gaussian likelihood.
* **T–T fits** (`notebooks/tests_analysis/` Part E, was notebook 14; `10_paper_plots.ipynb`): `sigma_for` on both axes, including Haslam (10% ⊕ 3 K).
* **Sky models against TRIS** (`notebooks/tests_analysis/` Part E Plots 6a–6d and Part B, were notebooks 14 and 12; notebook 13 removed 2026-11-16): only the TRIS `ring_sigma`; the models carry
  no error.

### What the error model leaves out

* **Correlations.** The scale and zero-level terms are the same offset or factor for every point of a survey (and the
  TRIS zero level for the whole ring), but each point is treated as independent (TODO N-15). Ring samples 3° apart
  under a ~20° beam are also correlated, so fitted uncertainties are underestimated by roughly ×2.6.
* **Noise after beam matching.** The noise term is the survey's per-pixel rms at its native resolution, applied
  unchanged to the beam-averaged ring TOD. Averaging over the ~20° TRIS beam would reduce it a great deal; it is
  negligible next to the scale and zero-level terms anyway (17–30 mK).
* **Main-beam vs full-beam scale.** Every map is treated as full-beam; for LWA1 and Maipu, which state no scale, that is
  an assumption (TODO N-19).
* **Stockert's weight.** At 1420 MHz the 0.5 K zero level is about half the CMB-free sky on the ring, so Stockert
  carries little weight in any fit.

## Irfan & Puglisi 2025 — synchrotron spectral index and curvature maps

| | |
|---|---|
| **Files** | `external/irfan_puglisi_2025/final_beta_1deg.npy` (β_s at 45 MHz), `final_cs_1deg.npy` (curvature c_s); git-ignored |
| **Paper** | Irfan & Puglisi 2025 (`Irfan_and_puglisi.pdf` in the repo root): eqs. 1–2, model T(p, ν) = T(p, ν₀) (ν/ν₀)^(β_s(p) + c_s(p) ln(ν/ν₀)), ν₀ = 45 MHz; 1° resolution, Nside 256 (their §5) |
| **Format** | 786 432 float64 values = HEALPix nside 256. Ordering and frame are not stated in the paper; **RING** (neighbouring pixels 7× smoother than read as NEST) and **Galactic** (small-scale structure correlates with Haslam's at +0.23 for β_s, −0.18 for c_s read as Galactic; +0.02 / −0.02 read as equatorial), determined 2026-11-17 |
| **Licence** | not in the repo; to be cleared with the authors (`HANDOVER.md` §5) |
| **Used by** | the pipeline's Irfan correction off 408 MHz (`pipeline.irfan_scaling`, steps 2 and 7; `run_pipeline.yaml` `pipeline.sky_models.irfan_*`); `notebooks/tests_analysis/` Part J (β comparison along the ring) |

## Not on LAMBDA

* **Berkhuijsen (1972), 820 MHz** — the one genuinely Haslam-independent prior basis, and
  the open item blocking every interpretable result (`TODO.md` B1). Not available from
  the archives above; `paths.berkhuijsen_map` and `haslam.prior.source:
  berkhuijsen_1972` are already implemented and waiting for the file.
* **GSM2008** — not a download at all: generated by `pygdsm`'s `GlobalSkyModel`, which
  fetches its own coefficient data into the astropy cache on first use. Stand-in prior
  and stand-in "truth" for stages 1 and 2.
