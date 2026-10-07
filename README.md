# moonrunes

**TRIS map-making (`limTOD`) feeding a Gibbs recalibration of the Haslam 408 MHz map
(`bayesian_skymap`).**

The idea that ties the two codebases together is one sentence: `bayesian_skymap` needs a
high-resolution but badly-calibrated map (Haslam) *and* a low-resolution but **absolutely
calibrated** one — and TRIS is an absolute radiometer. So `limTOD.tris` makes maps from
the TRIS drift rings, and those maps become `bayesian_skymap`'s `data1`/`sky_lowres`.

`PIPELINE.md` describes the pipeline in detail as it currently stands.
`tris_haslam_pipeline.md` is the algorithm, stage by stage, and is the document this repo
implements. `configs/run_config.yaml` is the single source of truth for every decision it
flags as "decide rather than default".

---

## Status

> **The design changed.** `tris_haslam_pipeline.md` now specifies a **pixel-by-pixel SED
> fit** using TRIS + ARCADE 2 to calibrate Haslam, superseding the `bayesian_skymap` Gibbs
> sampler this README was written around. Stage 1 carries over unchanged; stage 2 has been
> replaced; the sections below that describe the Gibbs design, its three blockers and its
> prior sourcing are **out of date** and are kept until the rewrite catches up.

| Stage | What it does | State |
|---|---|---|
| **1** | Calibrated TRIS maps from the public archive (`limTOD.tris`) | **implemented** |
| **2** | Match Haslam and ARCADE 2 to the TRIS beam, frame and grid (Step 0) | **implemented**: writes `outputs/stage2/` |
| 3 | Assemble per-pixel SEDs and fit | prototype in `notebooks/06`, not in the CLI |
| 4 | Derive and apply the calibration correction | not started |
| 5 | Validate | not started |

---

## Setup

```bash
conda env create -f environment.yml     # or use an existing env with the same stack
conda activate tris-haslam-cal
pip install -e .
git submodule update --init             # external/bayesian_skymap
```

`DATA_SOURCES.md` is the download list: every external dataset, its LAMBDA page, and
which stage reads it.

Stage 1 additionally needs:

* **the TRIS archive** — the four public LAMBDA text files
  (`TRIS_absolute_600.txt`, `TRIS_absolute_820.txt`, `TRIS_absolute_2500MHz.txt`,
  `TRIS_Beam_Profile.txt`). Point `paths.tris_archive_dir` at your local copy; nothing
  in this repo touches the network.
* **`pygdsm`**, for the GSM2008 template (Blocker 2's stand-in prior). Stage 2 also
  takes the reprocessed Haslam map (destriped/desourced, Remazeilles et al. 2015) from
  the same astropy cache, so `paths.haslam_map` can stay `null`.
* **`limTOD`** — a normal pip dependency, but if it is not installed the pipeline falls
  back to the sibling source checkout at `paths.limtod` (`../limTOD` by default), which
  is what the walkthrough notebook uses.

---

## Running

```bash
python run_pipeline.py --stage 1
python run_pipeline.py --stage 2
```

Stage 1 takes two overrides, so a variant run needs no config edit. Whatever they are
set to is what that run's manifest records:

```bash
python run_pipeline.py --stage 1 --nside 32 --run-name fine --output-dir outputs/stage1_nside32
```

`--nside` overrides `tris.nside` (power of two; `nside_hires` is lifted with it if
needed) and `--run-name` overrides `run.name`, which names the product
`tris_maps_<run name>.npz`. Note the **manifest filename is fixed** at
`stage1_manifest.json`, so two runs in the same directory overwrite each other's
manifest — give each variant its own `--output-dir`.

### Stage 1

The whole stage runs in about 12 s (the two Krylov solves are 0.45 s of
that; the rest is building the operator and the GSM2008 template). Products land in `outputs/stage1/` (gitignored):

* `tris_maps_<run name>.npz` — the maps, their posterior σ, the prior, the reference
  template, the fitted zero levels, residuals and beam coverage, stacked over frequency.
* `stage1_manifest.json` — the exact config block used, the solver diagnostics, and the
  three blockers' resolutions, so a product can always be traced back to the decisions
  that produced it.

Read them back through the one contract everything downstream uses:

```python
from moonrunes.stage1_tris_maps import load_stage1
products = load_stage1()
products["sky_k"]                    # (Nfreq, npix_band)
products.diagnostics(600)["solver"]  # info, iterations, agreement with the exact solve
```

---

## What stage 1 actually does

The forward model is the one **verified numerically** in limTOD's own walkthrough
(`examples/TRIS/tris_limtod_walkthrough.ipynb`), and nothing here re-derives it: the beam
comes from the archive's E/H principal-plane cuts (HPBW 19.155° / 23.366°), the pointing
is the parked-zenith geometry (`RA = LST`, `dec = lat = 42°26′`) with the E plane rolled
7° east via `selfrot = −7`, below-horizon response is masked, and temperatures are
Rayleigh–Jeans **including** the CMB monopole. The zero level is carried as an explicit
nuisance parameter with a Gaussian prior — exactly equivalent to the rank-1 noise term,
and it hands back the fitted offset.

**The solve is the direct Krylov solve, per Blocker 3.** Not a dense factorisation, and
not the Woodbury acceleration that is documented wrong by eight orders of magnitude at
low channel count: a `scipy.sparse.linalg.LinearOperator` around the MAP normal-equation
matvec plus a Jacobi-preconditioned GMRES, the same pattern as `bayesian_skymap`'s
amplitude step, so stage 1 and stage 4 fail the same way if they fail at all.

Because *"it finished"* is not *"it converged"* — this project's own history has silently
divergent solves accepted as valid samples — every solve is cross-checked against the
exact Woodbury identity at a single frequency (the inner matrix is only 120×120, so it is
cheap and is **not** the multi-frequency approximation Blocker 3 rules out). The same
identity supplies the posterior σ, which a Krylov solve alone cannot give. Stage 1
refuses to write a product that fails either check.

### Results at a glance

```
600 MHz (νeff 600.5): limtod solve, 5.4e-10 σ from exact, 0.07 s
        beam response > 0.01 of peak: 1272 pixels (41.4% of sky), 99.61% of beam power
        reduced χ² 16.14   zero level +1.7369 ± 0.0428 K
820 MHz (νeff 817.8): limtod solve, 1.5e-10 σ from exact, 0.06 s
        reduced χ² 5.88    zero level +0.5158 ± 0.0271 K
```

These come from the current config, where **everything is nside 16**, the grid the TRIS
stage 1 maps are made on, so every other product is put on the TRIS grid. That includes
`tris.nside_hires`, the grid limTOD rotates the beam on. That last choice has a measured
cost. The same run with `nside_hires: 64` gives reduced χ² **3.09 / 3.22**, and the data
tighten **388 / 201** pixels by more than 5%, against **23 / 7** at 16. The 64 numbers are
the product `TRIS_MAP_PROVENANCE.md` documents. Stage 1 also warns that `beam.nside_new`
(16) differs from what `bayesian_func.nside_for_beam` derives for the TRIS beam; that
comparison only mattered for the superseded Gibbs path, so it is recorded, not enforced.

The solve is **limTOD's own** (`TRISMapMakingInputs.solve`), on a deliberately coarse
nside 16 grid, over the pixels the beam actually saw. `TRIS_MAP_PROVENANCE.md` has the
full settings and the threshold convergence study.

Those fitted offsets reproduce the GSM2008 deficits measured independently in
`notebooks/01` (1.92 K and 0.56 K cold) — a real check that the geometry, beam and
convention translation came through intact.

**Read them as a measurement of GSM2008, not of the archive's zero point.** The monopole
degeneracy is ≈1, so with a per-pixel prior this tight the nuisance offset is the only
place a template error can go. This is the honest limit of one absolutely-calibrated ring
plus a template that is kelvins off, not a defect in the map-maker.

### Stage 2 — beam matching

Step 0 of `tris_haslam_pipeline.md`. Haslam and both ARCADE 2 bands are matched **to
TRIS**: its beam, its frame and its grid. An SED is only meaningful if every point in it
came through the same beam. The TRIS maps themselves are not touched. Per map:

1. read it, and check the frame the file declares (`moonrunes.frames.declared_frame`);
2. rotate it into equatorial, the frame of the TRIS maps: harmonic space for Haslam (full
   sky), mask-aware pixel space for ARCADE 2;
3. smooth it to the TRIS beam (`beam_matching.target_fwhm_deg` = 23.366°, the larger TRIS
   axis), with the quadrature kernel √(target² − native²), mask-aware for ARCADE 2;
4. keep a pixel only if it was **observed** and its smoothed weight reached
   `beam_matching.weight_floor` (0.5);
5. regrid onto `beam.nside_new` (= `tris.nside` = 16). A coarse pixel gets a value only if
   all of its sub-pixels were observed (`beam_matching.coarse_pixel_rule`).

```bash
python run_pipeline.py --stage 2
```

```
haslam_408     G nside 512  native  0.9333° -> extra kernel 23.3474°
               observed 3145728 -> 3145728 pixels after smoothing (0 extrapolated removed); 3072 on the nside 16 grid
arcade2_3150   G nside 16   native 11.6000° -> extra kernel 20.2832°
               observed 220 -> 170 pixels after smoothing (4 extrapolated removed); 170 on the nside 16 grid
arcade2_3410   G nside 16   native 11.6000° -> extra kernel 20.2832°
               observed 234 -> 192 pixels after smoothing (5 extrapolated removed); 192 on the nside 16 grid
```

**The product** is `outputs/stage2/beam_matched_<run name>.npz` plus `stage2_manifest.json`
(the config used, each file's declared frame, the rotation method, kernels and pixel
counts). Read it back with:

```python
from moonrunes.stage2_beam_matching import load_stage2
products = load_stage2()
products.map("arcade2_3150")        # nside 16, equatorial, hp.UNSEEN where unobserved
```

On the TRIS grid, rotated ARCADE 2 correlates with rotated Haslam at 0.98, and Haslam's
mean is preserved exactly (34.405 K). The weight floor costs ARCADE 2 the narrow
junctions of its observed ring (`TODO.md` N-7).

---

## What the old stage 2 did *(superseded)*

> This section describes `stage2_haslam_prep`, the Gibbs-design stage 2, which has been
> deleted. Its config blocks are kept in `configs/run_config.yaml` under a SUPERSEDED
> banner, and its circularity measurement is still worth reading in `notebooks/02`.

**The frame is galactic** and everything from here on lives in it. Haslam is native there,
so the high-resolution map is never rotated, and the gain and β regions become galactic
latitude bands — which is what `bayesian_skymap`'s own `region_operator` assumes and what
is physical for synchrotron β. Stage 3 rotates stage 1's equatorial TRIS band in.

The CMB monopole is removed from Haslam before anything else. `bayesian_skymap` models
`sky_gain` as `gain × s × (freq/ref)^β`, a pure power law with no monopole term, so leaving
2.7 K in it makes the forward model wrong everywhere.

**12 gain regions, 6 β regions.** This corrects what looks like a conflation in the
original config, which set 6 gain regions "to match the six regions bayesian_skymap
initialises" — those six are the *spectral index* regions (`BETA_REGIONS`); that
repository's own generator uses `ngain = 12`. The gain is meant to have more freedom than
β. Our `region_operator` is pinned byte-for-byte against theirs in the test suite.

### ⚠️ The prior is not independent, and that is the whole problem

Blocker 2 wants Berkhuijsen (1972) at 820 MHz because it is the one basis **independent of
Haslam**. It is not available offline, and this repo never touches the network. Three
sources are implemented and switched by `haslam.prior.source`:

| source | independent? | what it is |
|---|---|---|
| `gsm2008` *(default)* | **no** | GSM2008 at 408 MHz. It is built by fitting principal components to eleven input surveys, **one of which is Haslam 408 MHz** — so at 408 MHz it is a smoothed reconstruction of Haslam. `notebooks/02` measures the damage: log-space correlation **0.993**, median ratio **1.005**, 80% of pixels within 5%. |
| `tris_stage1` | **yes, where TRIS looked** | Stage 1's TRIS maps extrapolated 600 → 408 MHz. TRIS is an absolute radiometer, so this is genuinely independent — but only across the 52% of the sky its declination band covers. Outside it the width falls back to 10× the local sky, i.e. deliberately uninformative. |
| `berkhuijsen_1972` | **yes** | The real answer. Raises with a clear message until `paths.berkhuijsen_map` points at a file. |

Stage 2 sets `prior_is_haslam_derived` in its manifest and prints a warning when the prior
is circular. **With the shipped default, stage 5's comparison of a recovered gain against
the ~3% / ~0.91 K literature benchmark tests the plumbing, not the sky.** That is a
deliberate choice — it lets stages 3–5 be built and exercised — but a number that comes out
of it is not a result.

Two other truth-informed quantities have recorded substitutes: `covG`, which `gibbs_*.py`
builds from the *true* gain as `(0.15 × (1 + g_true))²` (stage 2 emits the same width at
`g = 0`), and `sky_haslam`, which is read as the true sky (stage 2 emits `s0_init_k`).

---

## Notebooks

All are committed **with their outputs**, so they can be read without re-running.
**03, 05 and 06 predate the switch to nside 16** and have not been re-run on it yet; 01's sections 05–06 follow `tris.nside`.

| Notebook | Contents |
|---|---|
| `01_explore_data.ipynb` | The inputs, one section per dataset. **01** the four TRIS products, the beam, the temperature convention, the GSM2008 template — and where `zero_level_sigma_k`, `uncertainty_floor_k` and `nside_new` come from. **02** a first look at ARCADE 2 (3.15/3.41 GHz). **03** the Haslam 408 MHz download. **04** the stage 1 TRIS maps. **05** every map in `res/`, discovered automatically: full-sky panels and Dec +30–55° strip slices with the TRIS strip, plus a frequency / frame / nside / beam table. **06** the TRIS TOD each map predicts (limTOD `generate_TOD_sky`), compared by shape and by spectral index against the TRIS 600 MHz ring. |
| `02_build_berkhijsen_prior.ipynb` | **Superseded, does not run** — imports the deleted `stage2_haslam_prep`. Kept for section 3's circularity measurement. The stage 2 products: the Haslam map, the region operators, the three prior sources compared — and section 3, which *measures* the GSM2008/Haslam circularity instead of asserting it. |
| `03_diagnostics.ipynb` | The stage 1 results: convergence, χ² and residuals, the zero level against the template deficit, prior→posterior shrinkage and z-scores, the band maps, β between the two solved maps, a preview of the degrade to `nside_new`, the prior-only and high-noise comparison runs, residuals with error bars (with the derivation of their expected scatter) and the RMS maps. |
| `04_coordinate_frames.ipynb` | What each input's header declares, checked against the sky; the rotation to equatorial (harmonic for Haslam, mask-aware pixel rotation for ARCADE 2); Sgr A*, Cas A and Tau A at their literature positions; the ARCADE 2 zero mask after rotation. |
| `05_beam_matching.ipynb` | Haslam and ARCADE 2 smoothed to 23.366° in the equatorial frame; UNSEEN counts and the ring's break at its narrow junctions; the final ARCADE 2 mask and the four-way overlap with TRIS. |
| `06_sed_fit.ipynb` | Stage 3 prototype with the four trusted points only (TRIS 600/820 + ARCADE 2), at nside 16 from the stage 1 and 2 products, fitted at **every** TRIS-stripe pixel (4 points where ARCADE 2 exists, 2 elsewhere): full maps of best-fit Haslam, β, χ², σ(β) and points per fit, and the sensitivity to ẑ and TRIS beam matching. With ARCADE 2: +3.4 K (+9.6%); TRIS only: +31%, driven by ẑ. |
| `07_sed_calibration.ipynb` | Stage 3 prototype on nside 16, with every usable map in `res/`: the **TOD route** (limTOD TODs → per-sample SED fit → correction TOD → map via the Wiener filter) and the **map route** (per-pixel fit on the TRIS grid), old vs corrected Haslam, and a hybrid map. The TOD route finds Haslam ~18% low along the ring; the map route is circular through stage 1's GSM2008 prior. |
| `08_paper_plots.ipynb` | Paper figures, target = destriped-only Haslam: the maps used and the target with the TRIS stripe, stripe slices, the TRIS maps and TODs, pixel-by-pixel and sample-by-sample SED fits with residuals, fit-residual maps, old vs best-fit Haslam (both routes), χ² and β maps, full-resolution corrected Haslam (both routes), T–T corner plots (raw and beam-convolved). Saved to `outputs/paper_plots/`. |

---

## The three blockers

| | Status |
|---|---|
| **1 — beam width mismatch** | Resolved. The fork accepts arbitrary `beam_deg`; the TRIS E/H mean (21.2605°) is passed through as itself. Stage 1 still compares the pinned `nside_new` with `bayesian_func.nside_for_beam`, but since the SED redesign a mismatch is recorded in the manifest and warned about, not fatal. `nside_new` is now 16. |
| **2 — no "truth" for real data** | **Open, in both stages.** Stage 1 uses GSM2008 as both prior and reference truth (`prior_equals_truth: true`). Stage 2's prior defaults to GSM2008 too, which is not Haslam-independent (`prior_is_haslam_derived: true`, measured at r = 0.993 in `notebooks/02`). `tris_stage1` is the honest interim — independent across the 52% of sky TRIS observed. The real answer stays Berkhuijsen (1972), and swapping it in is one config line. |
| **3 — Woodbury fast path broken at low channel count** | Resolved by not using it. `gibbs.use_woodbury` is `false` and stage 1 solves with the direct Krylov path. |

---

## Layout

```
PIPELINE.md                 the whole pipeline in detail, as it currently stands
TODO.md                     open items from stages 1 and 2, ordered by what they block
DATA_SOURCES.md             every external dataset, its archive page, and what reads it
TRIS_MAP_PROVENANCE.md      the limTOD code, settings and conventions behind the stage 1 maps
configs/run_config.yaml     every flagged decision, with the reason next to it
run_pipeline.py             the single entry point (--stage N)
src/moonrunes/
    config.py               loader; `require()` raises on null = "still an open decision"
    stage1_tris_maps.py     stage 1 + its product loader
    stage2_beam_matching.py stage 2 -- Step 0 beam matching
    frames.py               declared frames; Galactic -> equatorial rotation (full-sky and mask-aware)
    sky_maps.py             reads every map in res/ (any layout), with frequency, frame, beam and CMB convention
    sed_prototype.py        stage 3 prototype: the TOD route and the map route of the SED calibration
    stage3..stage5          empty
notebooks/                  analysis and plotting, committed with outputs
external/bayesian_skymap    submodule (patched fork)
outputs/                    stage products, gitignored
res/                        input maps: TRIS text products, Haslam, ARCADE 2, Maipu 45 MHz, Stockert 1420 MHz
tests/                      pytest: test_stage_io.py, test_frames.py
```

Two rules the stages rely on: a `null` in the config is an *open decision* and the stage
that needs it raises rather than defaulting; and every stage records the config block it
consumed in its manifest.

---

## Tests

```bash
pytest
```

49 tests. The algebra ones — the Krylov solve against a dense reference, the Jacobi
diagonal against `bayesian_func.estimate_diag_precond`, the Woodbury posterior against a
dense inverse, stage 2's beam quadrature and masked-smoothing behaviour, and the frame rotations — run
anywhere. The end-to-end tests skip themselves when the TRIS archive, the FITS downloads,
or `pygdsm` is absent.

---

## License

MIT — see `LICENSE`.
