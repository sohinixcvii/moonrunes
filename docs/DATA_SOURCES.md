# Data sources

Every dataset the code reads: where it comes from, its format, what reads it, and the uncertainty it is given. Unless
noted, the maps are hosted by **LAMBDA** (NASA/GSFC's Legacy Archive for Microwave Background Data Analysis). Every map
in `res/` is read by `moonrunes.sky_maps.discover`, which records its frequency, unit, frame, beam and CMB convention
(from the file, or from `sky_maps.EXTERNAL_METADATA` where the file does not say). The code reads only local files;
the one exception is `pygdsm`, which downloads its sky-model data into the astropy cache on first use.

---

## TRIS — the absolutely calibrated rings

<https://lambda.gsfc.nasa.gov/product/tris/tris_prod_table.html>; Zannoni et al. 2008 (ApJ 688, 12), Tables 10–11.

| | |
|---|---|
| **Files** | `TRIS_absolute_600.txt`, `TRIS_absolute_820.txt`, `TRIS_absolute_2500MHz.txt`, `TRIS_Beam_Profile.txt` (copies in `res/`; the code reads `paths.tris_archive_dir`, `../limTOD/downloads/TRIS`) |
| **Format** | plain text: 120 samples per ring (RA every 3°, Dec +42°), temperature (K, CMB included), per-sample statistical error, and the zero-level uncertainty in the header; the beam file has the E- and H-plane cuts |
| **Read by** | `limTOD.tris` readers: step 0 (the maps), step 2 (the ring data), step 7 (the geometry) |

- **Effective frequencies** 600.5 and 817.8 MHz (the files say 0.6 and 0.82 GHz): from limTOD's table
  (`../limTOD/limTOD/tris/archive.py`), attributed there to Zannoni et al. 2008 Sect. 2 and Tab. 1.
- **The beam:** the archive cuts give half-power widths of 19.155° (E plane) and 23.366° (H plane); the E plane is
  rolled 7° east of the meridian. The `tris_cuts` kernel interpolates between the two cuts.
- **Zero level:** 0.066 K at 600 MHz; asymmetric, +0.430 / −0.300 K, at 820 MHz.
- **One sample per ring has a statistical error of 0.000 K**; the code replaces it with the ring's smallest positive
  value (0.004 K).
- **The 2.5 GHz file** (six points, not a ring) is not used.

## Haslam 408 MHz — the map being tested

<https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_get.html>; Haslam et al. 1982; Remazeilles et al. 2015.

| | |
|---|---|
| **Files** | `res/HASLAM_DESTRIPED_ONLY.fits` (Remazeilles et al.'s reprocessing, destriped but not desourced: the pipeline's target, `pipeline.target_file`); `res/HASLAM.fits` (used by `--stage 2` and the beam tests) |
| **Format** | HEALPix, Galactic, RING, nside 512, 56′ beam, K (RJ, CMB included) |
| **Read by** | step 2 (Haslam's ring TOD; never a fit point), step 6 (compared with the prediction), the T–T fits, the beam tests |

GSM2008 is built from surveys that include this map, so GSM2008 is not independent of Haslam.

## The fitted surveys

| Survey | File(s) | Format | Beam | Reference |
|---|---|---|---|---|
| LWA1 35, 38, 40, 45, 50, 60, 70, 74, 80 MHz | `LWA_35MHz.fits`, `LWA1_38.fits` … `LWA1_74.fits`, `LWA_80MHz.fits` | HEALPix, equatorial | 4.65° at 35 MHz to 2.04° at 80 MHz (from the files) | Dowell et al. 2017; LAMBDA `fg_lwa1_radio_maps_info` |
| Maipu + MU 45 MHz | `MAIPU_45mhZ.fits` | HEALPix, Galactic, NESTED, nside 64, K; unobserved = −32768 (caps around both poles; data Dec −79° to +71°) | 4.6° × 2.4° (Maipu) / 3.6° (MU); 4.6° used | Guzmán et al. 2011; LAMBDA `fg_maipu_info` |
| Stockert + Villa-Elisa 1420 MHz | `STOCKERT_VILLA_1400MHZ.fits` (the survey is at 1420 MHz) | HEALPix, Galactic, NESTED, nside 256, mK (full-beam brightness temperature); 47 unobserved pixels | 35.4′ | Reich 1982; Reich & Reich 1986; Testori et al. 2001; LAMBDA `fg_stockert_villa_info` |

Neither the Maipu nor the Stockert header gives a frequency, unit or beam; those come from the LAMBDA pages
(`EXTERNAL_METADATA`). All three surveys include the CMB, which the code removes in Rayleigh–Jeans.

## Maps in `res/` that are not fitted

| Map | Why not | Notes |
|---|---|---|
| CHIPASS 1400 MHz (`CHIPASS_1400MHz.fits`) | southern: less than half of the TRIS beam falls on observed sky at every ring sample | 104.9 MB, over GitHub's limit, so git-ignored: download it from <https://lambda.gsfc.nasa.gov/product/foreground/fg_chipass_info.html> and save it under that name. HPX projection, mK, 14.4′ beam, CMB included (Calabretta et al. 2014); `sky_maps` converts it to HEALPix |
| Rhodes/HartRAO 2326 MHz (`RHODES_HARTRAO_2326MHz.fits`) | southern, as CHIPASS | CMB removed (measured: 0.17 K at \|b\| > 40°); Jonas et al. 1998 |
| ARCADE 2 3.15, 3.41 GHz (`ARCADE2_315.fits`, `ARCADE2_341.fits`) | in `sed.exclude_maps`: ~7 % of the sky as a thin ring, data at 16–17 of the 120 ring samples, and its 6–7.5 mK errors would let those few points pull the fits | HEALPix NESTED, nside 16, K thermodynamic with the CMB monopole (dipole removed), `BEAMSZ = 11.6°`; unobserved sky stored as exact 0 (~93 % of the pixels), and the frame given as `SKYCOORD='Galactic'` (not `COORDSYS`, so healpy reports none); Fixsen et al. 2011. `--stage 2` beam-matches it |
| Juno 600, 1248 MHz (`Juno_600MHz.fits`, `Juno_1248MHz.fits`; git-ignored) | in `sed.exclude_maps`: unit, frame, beam and provenance unknown | the header has only NSIDE (512) and ORDERING. **Do not redistribute** (`QUICKSTART.md` §5) |
| WMAP K, Ka (`WMAP_K.fits`, `WMAP_Ka.fits`) | unusable: the headers say V band, contradicting the file names, and the monopole is removed (no absolute zero level) | WMAP 7-year band maps, mK thermodynamic |

## Sky models

- **GSM2008, GSM2016, LFSM** — generated by `pygdsm` (`GlobalSkyModel`, `GlobalSkyModel16`, `LowFrequencySkyModel`),
  CMB excluded, RJ K; native beams 1.0°, 0.8°, 5.1° (`pipeline.sky_models.native_fwhm_deg`). Step 0 also uses GSM2008
  as the prior template of the TRIS map-maker.
- **The corrections** (`pipeline.sky_models`): Carter (Carter et al. 2025: 1.029 T + 0.91 K at 408 MHz); Irfan
  (T_tot − 8.9 K at 408 MHz, Irfan & Puglisi 2025 eq. 6); GINAN (McKay et al. 2025, as cited in
  `configs/run_pipeline.yaml`; `pygdsm.mckay25.generate_mckay25_correction`, built at 350 MHz on GSM2016). How each is
  carried to other frequencies: `docs/PIPELINE.md` and `CHANGELOG.md`.

## Irfan & Puglisi 2025 — synchrotron spectral index and curvature maps

| | |
|---|---|
| **Files** | `external/irfan_puglisi_2025/final_beta_1deg.npy` (β_s at 45 MHz), `final_cs_1deg.npy` (curvature c_s); git-ignored |
| **Model** | T(p, ν) = T(p, ν₀) (ν/ν₀)^(β_s(p) + c_s(p) ln(ν/ν₀)), ν₀ = 45 MHz (their eqs. 1–2); 1° resolution, nside 256 (their §5) |
| **Format** | 786 432 float64 values = HEALPix nside 256. Ordering and frame are not stated in the paper; **RING** (neighbouring pixels are 7× smoother than read as NEST) and **Galactic** (small-scale structure correlates with Haslam's at +0.23 for β_s and −0.18 for c_s read as Galactic, against +0.02 / −0.02 read as equatorial), determined from the data |
| **Licence** | not cleared with the authors; **do not redistribute** |
| **Read by** | `pipeline.irfan_scaling` (steps 2 and 7: the Irfan correction off 408 MHz); Part J (β along the ring) |

## RHINO beam and horizon mask

`../rhino-beams/npz_beams/DryHornGround.npz` and `../rhino-beams/horizon_masks/JodrellBank.fits`, outside the repo
(`pipeline.rhino.beam_file`, `horizon_file`), read by step 7 only. Not to be redistributed.

## Not available

- **Berkhuijsen 1972, 820 MHz** — the one high-resolution map at a TRIS frequency independent of the 408 MHz survey;
  not available in digital form.

---

## Uncertainties for every map

No map file carries an error column, so every survey's uncertainty comes from its publication; TRIS's from the archive.

### The error model

For a survey point (a ring sample) with CMB-free Rayleigh–Jeans temperature T:

    σ = sqrt( (scale · |T|)² + zero_level² + noise² )

(`sed_prototype.sigma_for`), with scale the temperature-scale (gain) error as a fraction of T, zero_level the error of
the map's absolute offset (K), and noise the rms receiver noise (K); a term the publication does not quote is 0. The
values are in `sed_prototype.PUBLISHED_ERRORS` (one entry per file, with its source and any judgement call) and in
`configs/error_bars.csv` (the same numbers, read by the model comparison). A map with no entry would fall back to 10 %
of T (`ASSUMED_FRACTION`); no fitted map uses it.

For TRIS (`sed_prototype.ring_sigma`), per ring sample:

    σ = sqrt( statistical² + zero_level² )

with the statistical error floored as above and the larger side of the zero-level uncertainty (0.43 K at 820 MHz).

### Per map

| Map | Scale | Zero level | Noise | Median σ on the ring | Source, and judgement calls |
|---|---|---|---|---|---|
| LWA1 35–80 MHz (all nine) | 20 % | — | — | 20 % | Dowell et al. 2017: "a conservative estimate of the calibration accuracy … 15–25 per cent", for structure and sky average. **20 % is the middle of that range.** The per-pixel uncertainty maps of that paper are not used |
| Maipu + MU 45 MHz | 1.7 % | 39 K | — | 1.8 % | Guzmán et al. 2011 quote no error. **The 1σ halves of the corrections Monsalve et al. 2021 find against EDGES** (×1.076 ± 0.034 and −160 ± 78 K, 2σ). **The corrections themselves are not applied** |
| Haslam 408 MHz | 10 % | 3 K | — | 14 % | Haslam et al. 1982. Never a fit point |
| TRIS 600 MHz | — | 0.066 K | per sample (0.004–0.020 K) | 0.067 K (0.5 %) | the archive; the zero level dominates |
| TRIS 820 MHz | — | 0.43 K (larger side of +0.43 / −0.30) | per sample (0.004–0.041 K) | 0.430 K (8.5 %) | the archive |
| Stockert + Villa-Elisa 1420 MHz | 5 % | 0.5 K | 17 mK | 47 % (0.50 K of 1.07 K) | Testori et al. 2001: gain ≈ 5 %, zero level 0.5 K, sensitivity 50 mK = 3 × rms, so noise = 50 mK / 3 |
| CHIPASS 1400 MHz | 5 % | 0.5 K | 30 mK | — | Calabretta et al. 2014 §4 ("a typical rms error of 30 mK"; zero level tied to the 1420 MHz map). **No scale or zero-level error of its own, so Stockert's are used** |
| Rhodes/HartRAO 2326 MHz | 5 % | 80 mK | 30 mK | — | Jonas et al. 1998: scale < 5 %, zero level < 80 mK, rms < 30 mK. **The upper limits are used as 1σ** |
| ARCADE 2 | — | — | 7.5, 6.0 mK | — | Fixsen et al. 2011 Table 3 (excluded from the fits) |

The step 0 TRIS **maps** (map route only) use the Wiener posterior σ per pixel, plus the fitted zero level's σ when
the zero level is added (`configs/model_comp_config.yaml` `zero_level.add_zhat_to_tris`).

### What the error model leaves out

- **Correlations.** The scale and zero-level terms are one factor or offset for all of a survey's points (and the TRIS
  zero level for the whole ring), but each point is treated as independent. Ring samples 3° apart under a ~20° beam
  are also correlated: there are about 18 independent beams on the ring, so formal uncertainties are underestimated by
  roughly √(120/18) ≈ 2.6, and summed evidences are overstated (`QUICKSTART.md` open item 1).
- **Noise after beam matching.** The noise term is the per-pixel rms at native resolution, applied unchanged to the
  beam-averaged TOD; it is negligible next to the scale and zero-level terms (17–30 mK).
- **Main-beam vs full-beam scale.** Every map is treated as full-beam; for LWA1 and Maipu, which do not say, that is an
  assumption.
- **Stockert's weight.** At 1420 MHz the 0.5 K zero level is about half the CMB-free sky on the ring, so Stockert
  carries little weight in any fit.
