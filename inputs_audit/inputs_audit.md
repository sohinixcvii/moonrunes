# Inputs audit: the SED fit and the zero level

> **Update 2026-10-21:** the 10% placeholder described below has been replaced by each
> survey's published scale, zero-level and noise errors (`sed_prototype.PUBLISHED_ERRORS`;
> table in `DATA_SOURCES.md`, "Survey uncertainties used in the SED fit"). The rest of this
> audit describes the code as it stood on 2026-10-20, and its line numbers in
> `sed_prototype.py` are now offset by the new table.

Audited 2026-10-20 on branch `fwhm-scaling-test`, from the working tree. No code was run
except read-only loads to read values, which are marked **(read from products)**. Line
numbers are for the files as they stand today.

**Scope.** "The SED fit" means the current stage 3 prototype,
`src/moonrunes/sed_prototype.py`. It has a **TOD route** (`tod_route`, fitted per TRIS ring
sample) and a **map route** (`map_route`, fitted per TRIS-stripe pixel). Notebooks 08 and 10
call it, and notebook 07 holds a step-by-step copy of the same calculation. Notebook 06 (the
four-point TRIS + ARCADE 2 fit) is mentioned only where it differs.

Claims marked **[inferred]** are my reading of the code, not something the code or docs
state outright.

---

## 1. Surveys and maps in the SED fit

**How the list is built.** Every `*.fits` in `res/` is read by `sky_maps.discover`
(`sky_maps.py:286-295`). The fit uses the maps with `usable` set, `role == "fit"` and not in
`sed.exclude_maps` (`sed_prototype.py:113-115`; `configs/run_config.yaml:449`). The TRIS
points come from the archive rings (TOD route, `sed_prototype.py:124-126, 255-258`) or the
stage 1 maps (map route, `sed_prototype.py:140-147`).

**The target, Haslam, is not a fit point.** `Inputs` takes
`target_file="HASLAM_DESTRIPED_ONLY.fits"` (`sed_prototype.py:93`). Its role is `target`
(`sky_maps.py:47-49`), so it is compared with the fit's prediction, never fitted.

| # | Survey / map | Frequency used | Reference (as recorded in the repo) | In the fit? |
|---|---|---|---|---|
| 1 | TRIS 600 MHz | 600.5 MHz effective (`ring.effective_frequency_mhz`, `sed_prototype.py:256`; map route `effective_freq_mhz`, `:143`) | LAMBDA TRIS products (`DATA_SOURCES.md:13-26`); no paper is cited in the repo | yes |
| 2 | TRIS 820 MHz | 817.8 MHz effective (same lines) | same | yes |
| 3 | LWA1 35 MHz (`LWA_35MHz.fits`) | 35.0, header `RESTFREQ` (`sky_maps.py:173`) | Dowell et al. 2017; LAMBDA `fg_lwa1_radio_maps_info` (`sky_maps.py:54-56`) | yes |
| 4 | Maipu + MU 45 MHz (`MAIPU_45mhZ.fits`) | 45.0, from the **filename** (`sky_maps.py:187-189`) | Guzmán et al. 2011, A&A 525, A138 (`sky_maps.py:50-53`; `DATA_SOURCES.md:124`) | yes |
| 5 | LWA1 80 MHz (`LWA_80MHz.fits`) | 80.0, header `RESTFREQ` | Dowell et al. 2017 (`sky_maps.py:57-59`) | yes |
| 6 | CHIPASS 1400 MHz (`CHIPASS_1400MHz.fits`) | 1400, from the **filename** | Calabretta et al. 2014 (`sky_maps.py:64-68`; `DATA_SOURCES.md:133`) | listed, but **no data on the TRIS ring** (`PIPELINE.md:911-912`) |
| 7 | Stockert + Villa-Elisa 1420 MHz (`STOCKERT_VILLA_1400MHZ.fits`) | 1420, external override; the filename says 1400 (`sky_maps.py:60-63`) | Reich 1982; Reich & Reich 1986; Testori et al. 2001 (`DATA_SOURCES.md:124`; `sky_maps.py:63`) | yes |
| 8 | Rhodes/HartRAO 2326 MHz (`RHODES_HARTRAO_2326MHz.fits`) | 2326, header comment (`sky_maps.py:182-186`) | Jonas et al. 1998; Platania et al. 2003 (`sky_maps.py:69-72`) | listed, but **no data on the TRIS ring** (`PIPELINE.md:911-912`) |
| — | ARCADE 2 3.15 / 3.41 GHz | 3150 / 3410, header `FREQ` | errors from Fixsen et al. 2011 Table 3 (`sed_prototype.py:38-39`); LAMBDA `arcade_maps_info` (`sky_maps.py:73-78`) | **excluded** since 2026-10-15 (`run_config.yaml:441-449`) |
| (target) | Haslam 408 MHz, destriped only | 408, header `FREQ` | Remazeilles et al. 2015, "destriped but not desourced" (`sky_maps.py:47-49`) | no: the target |

- **The WMAP maps are marked unusable** (`sky_maps.py:79-95`), so they never enter the fit.
- **[inferred] CHIPASS and Rhodes reach the TOD route only through the coverage cut.** That
  cut is `MIN_COVERAGE = 0.5` (`sed_prototype.py:40, 138, 261`). On the ring they are
  dropped by it. In the map route they would enter wherever the beam-matched map has a
  value on the stripe (`sed_prototype.py:208-209`); I have not checked how many stripe
  pixels that is.

---

## 2. Each input's uncertainty, beam scale and FWHM

### 2a. Where each uncertainty comes from

| Input | Uncertainty used in the fit | Code | Placeholder? |
|---|---|---|---|
| TRIS 600 / 820, **TOD route** (ring data) | Published per-sample statistical error, combined in quadrature with the published zero-level uncertainty. For the asymmetric 820 MHz value (+0.430 / −0.300 K) the **larger side, 0.430 K**, is used. The one row published as 0.000 K is replaced by the smallest positive error in that ring. | `ring_sigma`, `sed_prototype.py:80-86`; used at `:258` | Not a placeholder, but two **choices** are made: the floor (`:83`) and taking the larger side (`:85`). **[inferred]** No calibration-scale (gain) term is included: the ring object carries only `statistical_uncertainty_k` and `zero_level_uncertainty_k` **(read from products)**, and `ring_sigma` uses only those. |
| TRIS 600 / 820, **map route** (stage 1 maps) | The stage 1 posterior σ of each pixel, combined in quadrature with the posterior σ of ẑ | `sed_prototype.py:145-146` | Not a placeholder. **⚠ Flag:** the quadrature sum leaves out the pixel–ẑ covariance. `PIPELINE.md:766` gives σ² = Σ_pp + Σ_zz + 2Σ_pz, and notebook 06 computes that full form (cell 5, lines 29 and 33-35). |
| LWA1 35, Maipu 45, LWA1 80, CHIPASS 1400, Stockert 1420, Rhodes 2326 | **10% of the CMB-removed temperature**, `0.10 · |T_excess|` | `sigma_for`, `sed_prototype.py:72-77`, with `ASSUMED_FRACTION = 0.10` at **`sed_prototype.py:37`** ("Error assumed for maps with no published per-pixel error", `:36`). Used at `:211` (map route) and `:263` (TOD route) | **⚠ PLACEHOLDER**, the same 10% for every survey. A copy of the constant is in `notebooks/07_sed_calibration.ipynb` cell 2, line 32, and its `sigma_for` in cell 4, lines 14-18. |
| ARCADE 2 (currently excluded) | A fixed 7.5 mK (3.15 GHz) or 6.0 mK (3.41 GHz) | `ARCADE_CAL_K`, `sed_prototype.py:38-39`; `sigma_for` at `:75-76` | Not a placeholder: cited to Fixsen et al. 2011 Table 3. **[inferred]** It is only the calibration term. No per-pixel noise is added. |

Further flags:

- **⚠ Stockert: the docs and the code disagree.** `PIPELINE.md:910-911` says the pipeline's
  noise is "2.8% for Stockert". `sigma_for` has no Stockert case, so the code uses 10%
  (`sed_prototype.py:77`). Neither `src/` nor `configs/` contains a 2.8% value **[from
  searching both]**.
- **[inferred] Survey error bars scale only.** The placeholder is proportional to T, so it
  has no additive zero-level term for each survey, although each survey has its own
  published zero-level uncertainty.
- **[inferred] The error is a fraction of T *after* CMB removal.** `t_excess` comes from
  `sky_maps.rj_excess` (`sky_maps.py:310-317`), so the absolute error is slightly smaller
  than 10% of the map value.

### 2b. Main-beam vs full-beam scale

**It is not recorded anywhere in the code, for any map.**

- **`EXTERNAL_METADATA`** (`sky_maps.py:43-96`) has no such field, and `describe` returns no
  beam-scale key (`sky_maps.py:276-284`).
- **A search of `src/`, `configs/` and the top-level docs** for main beam, full beam, beam
  efficiency or η_mb finds **one** mention: `DATA_SOURCES.md:121`, which lists Stockert's
  units as "mK (full-beam brightness temperature)". It is documentation only, and no code
  reads it.

| Input | Main- or full-beam scale recorded? |
|---|---|
| TRIS 600 / 820 | no |
| LWA1 35 / 80 | no |
| Maipu + MU 45 | no |
| CHIPASS 1400 | no |
| Stockert 1420 | **only in docs:** full-beam (`DATA_SOURCES.md:121`); not in code |
| Rhodes 2326 | no |
| ARCADE 2 | no |
| Haslam (target) | no |

### 2c. Beam FWHM

**The rule:** `fwhm_deg` is the native beam, and for an elliptical beam it is the **larger
axis**, "so nothing is left under-smoothed" (`sky_maps.py:41-42`; parser
`sky_maps.py:193-211`, larger axis at `:204`).

**Where it is used:**
- **Map route:** `match_to_tris_beam(rotated, observed, m["fwhm_deg"], ...)`
  (`sed_prototype.py:171-173`).
- **TOD route: not at all.** `_native_tod` (`sed_prototype.py:186-199`) observes the native
  map through the TRIS beam after degrading it to nside 64 and then nside 16 (`:191-192,
  196`). **[inferred]** The survey's own beam is not modelled there, which is reasonable only
  while it is much narrower than the 19–23° TRIS beam.

| Input | FWHM used | From |
|---|---|---|
| TRIS 600 / 820 | none: TRIS is the target beam, 19.155° (E) × 23.366° (H) | `run_config.yaml:230-231`; the circular kernel uses 23.366° (`run_config.yaml:224`) |
| LWA1 35 | 4.65° | header `BMAJ`/`BMIN` (`sky_maps.py:201-204`) |
| Maipu + MU 45 | 4.6°, the larger axis of Maipu's 4.6° × 2.4°; MU is 3.6° | external, `sky_maps.py:51-52` |
| LWA1 80 | 2.04° | header `BMAJ`/`BMIN` |
| CHIPASS 1400 | 14.4′ = 0.24° | external, `sky_maps.py:65` |
| Stockert 1420 | 35.4′ = 0.59° | external, `sky_maps.py:61-62` |
| Rhodes 2326 | 20′ = 0.333° | header comment "Resolution" (`sky_maps.py:205-210`) |
| ARCADE 2 | 11.6° | header `BEAMSZ` (`sky_maps.py:198-200`); also `run_config.yaml:272-274` |
| Haslam (target) | 56′ = 0.933° | header `BEAMSIZE` (`sky_maps.py:195-197`); also `run_config.yaml:270-271` |

**[inferred]** `beam_matching.haslam_native_fwhm_deg` and `arcade2_native_fwhm_deg`
(`run_config.yaml:271, 274`) are not read by the SED code, which takes `m["fwhm_deg"]` from
the headers. The values agree today (0.9333° and 11.6°).

---

## 3. The zero-level offset ẑ

### Where it is defined

- **Stage 1 fits it.** It is an explicit nuisance parameter in the TRIS map-maker
  (`stage1_tris_maps.py:22` in the docstring, `:493-494` the zero-level column), with a
  Gaussian prior of σ = `tris.zero_level_sigma_k` = **0.5 K** (`run_config.yaml:125`, read
  at `stage1_tris_maps.py:456`, passed on at `:470`).
- **Fitted values** are `zero_level_k` and `zero_level_sigma_k` (`stage1_tris_maps.py:532-533`),
  written to the manifest (`:609-612`) and the product (`:848-851`).
- **Current product** (`outputs/stage1/tris_maps_tris_haslam_v0.npz`) **(read from products)**:
  ẑ = **+1.737 ± 0.043 K** at 600 MHz and **+0.516 ± 0.027 K** at 820 MHz, prior σ 0.5 K.
- **What it means:** "the fitted zero level absorbs a real model offset -- read it as a
  model-mismatch amplitude, not as the archive's zero point" (`run_config.yaml:148-151`).
  `PIPELINE.md:223-226` says the same: it is GSM2008's deficit, and the sky monopole and z
  are degenerate.
- **Not ẑ:** the 0.91 K "zero level" in `validate.zero_level_expected_k`
  (`run_config.yaml:471`) is Haslam's literature zero-level uncertainty, a separate quantity.

### The options (TODO N-3)

`TODO.md:175-178`, **N-3: decide whether ẑ is sky or instrument.**
- **As sky: add ẑ to the TRIS temperatures.** T = x̂_p + ẑ − T_CMB,RJ, with
  σ² = Σ_pp + Σ_zz + 2Σ_pz (`PIPELINE.md:764-766`).
- **As instrument: leave it out.** Use x̂_p and Σ_pp (`PIPELINE.md:767-768`).
- **What is at stake:** "worth about 6–7 K at 408 MHz". Together with N-2 it moves the
  median predicted − Haslam between −19% and +21% (`TODO.md:176-178`).
- **There is no config switch for it.** The only zero-level keys in `run_config.yaml` are
  `tris.zero_level_sigma_k` (`:125`) and `validate.zero_level_expected_k` (`:471`) **[from
  searching the config]**.

### What the code does by default

| Path | ẑ handling | Code |
|---|---|---|
| **Map route** (`sed_prototype`) | **Adds ẑ** to every TRIS pixel, hard-coded with no switch. σ is the quadrature sum, without the covariance (see 2a). | `sed_prototype.py:144-146`; points built at `:204-206` |
| **TOD route** (`sed_prototype`) | **Does not use ẑ.** It fits the archive's ring temperatures directly, and the archive's *published* zero-level uncertainty enters σ instead. | `:255-258`; `ring_sigma` `:84-86` |
| TOD route, mapping the correction | Builds limTOD's operator with a zero-level column, then **drops that column** (`[:, :-1]`). **[inferred]** The mapped correction therefore has no zero-level term. | `:274-280` |
| Notebook 07 | Repeats the map route's choice: adds ẑ, quadrature σ | cell 16, lines 4-5 |
| Notebook 06 | Builds **both** versions, `with_z` and `sky_only`, with the full-covariance σ | cell 5, lines 29, 33-35. Which one its fits use downstream was **not checked**. |
| Notebook 10 | The TRIS map figure prints ẑ "(not added)": display only | cell 3, lines 43-44 |

**Summary:** by default the map route treats ẑ as sky, and the TOD route avoids the
question by using the ring data. N-3 is open, and no single setting controls it.
