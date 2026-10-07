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

The reason the two codebases combine at all: `bayesian_skymap` needs an absolutely
calibrated map, and TRIS is an absolute radiometer. Stage 1 makes maps from the published
drift rings.

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
  (+0.430 / −0.300 K, against 0.066 K at 600 MHz); `limTOD` refuses to symmetrise that
  for you, which is why `tris.zero_level_sigma_k` is a single deliberately loose 0.5 K.

---

## Haslam 408 MHz — the high-resolution map being recalibrated

<https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_get.html>

The badly-calibrated high-resolution map (`data2` / `sky_gain`) whose gain and zero level
the Gibbs sampler is being asked to recover.

| | |
|---|---|
| **Consumed by** | Stage 2 (`moonrunes.stage2_beam_matching`), and notebooks 01, 04–06 |
| **Config** | `paths.haslam_map_fits` = `res/HASLAM.fits`. The superseded `stage2_haslam_prep` read `paths.haslam_map`, where `null` meant "use the copy `pygdsm` already caches" |
| **File** | `haslam408_dsds_Remazeilles2014.fits` — the **destriped, desourced reprocessed** map (Remazeilles et al. 2015), not the original 1982 release |
| **Format** | HEALPix FITS, `COORDSYS = GALACTIC`, RING, nside 512, `BEAMSIZE` = 56′, K (RJ, CMB included) |

Notes:

* `null` in `paths.haslam_map` is **not** an open decision (unlike other nulls in the
  config): it means stage 2 pulls the file from the astropy cache that `pygdsm` populates,
  which is the same reprocessed map. Set the path to override with a local download.
* The map is total-intensity, so stage 2 **removes the CMB monopole** (2.7157 K, the
  Rayleigh–Jeans value at 408 MHz) before anything else — `bayesian_skymap` models
  `sky_gain` as a pure power law with no monopole term.
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
| **Consumed by** | Stage 2 (`moonrunes.stage2_beam_matching`, `paths.arcade2_map_3150` / `_3410`), and the SED fit prototype in `notebooks/06` as two of the four trusted points |
| **Looked at in** | `notebooks/01_explore_data.ipynb`, sections 02 and 05; `notebooks/04`–`06` |
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
  (`notebooks/04_coordinate_frames.ipynb`).

One caveat that is design, not format: 3.15/3.41 GHz against Haslam's 408 MHz is close to
a decade of lever arm, over which free-free emission's rising contribution can put real
curvature in the SED that a single power law will not capture.

---

## Maipu + MU 45 MHz, and Stockert + Villa-Elisa 1420 MHz — reference maps

<https://lambda.gsfc.nasa.gov/product/foreground/fg_diffuse.html>

Two more all-sky surveys in `res/`. **No stage reads them.** They are shown alongside the
others in `notebooks/01_explore_data.ipynb`, sections 05–06. Neither header gives a
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

## Not on LAMBDA

* **Berkhuijsen (1972), 820 MHz** — the one genuinely Haslam-independent prior basis, and
  the open item blocking every interpretable result (`TODO.md` B1). Not available from
  the archives above; `paths.berkhuijsen_map` and `haslam.prior.source:
  berkhuijsen_1972` are already implemented and waiting for the file.
* **GSM2008** — not a download at all: generated by `pygdsm`'s `GlobalSkyModel`, which
  fetches its own coefficient data into the astropy cache on first use. Stand-in prior
  and stand-in "truth" for stages 1 and 2 (see Blocker 2).
