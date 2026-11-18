# moonrunes

**TRIS map-making (`limTOD`) and a TRIS-anchored spectral model of the Haslam 408 MHz map.**

TRIS is an absolute radiometer: its drift rings at 600 and 820 MHz anchor a spectral fit to
the other surveys along the TRIS ring, which tests the Haslam 408 MHz map and predicts the
sky at other frequencies.

`PIPELINE.md` describes the pipeline in detail as it currently stands.
`tris_haslam_pipeline.md` is the algorithm, stage by stage, and is the document this repo
implements. `configs/run_pipeline.yaml` is the single source of truth for every decision it
flags as "decide rather than default".

---

## Status

> **Design.** `tris_haslam_pipeline.md` specifies a **pixel-by-pixel SED
> fit** using TRIS + ARCADE 2 to calibrate Haslam. Stage 1 is the TRIS map-maker; stage 2 is
> the beam matching; `docs/PIPELINE.md` describes the end-to-end pipeline built on them.

| Stage | What it does | State |
|---|---|---|
| **1** | Calibrated TRIS maps from the public archive (`limTOD.tris`) | **implemented** |
| **2** | Match Haslam and ARCADE 2 to the TRIS beam, frame and grid (Step 0) | **implemented**: writes `outputs/stage2/` |
| 3 | Assemble per-pixel SEDs and fit | prototype in `notebooks/06` (removed 2026-11-16; in git at `a5c827c`), not in the CLI |
| 4 | Derive and apply the calibration correction | not started |
| 5 | Validate | not started |

---

## Setup

```bash
conda env create -f environment.yml     # or use an existing env with the same stack
conda activate tris-haslam-cal
pip install -e .
```

`DATA_SOURCES.md` is the download list: every external dataset, its LAMBDA page, and
which stage reads it.

Stage 1 additionally needs:

* **the TRIS archive** — the four public LAMBDA text files
  (`TRIS_absolute_600.txt`, `TRIS_absolute_820.txt`, `TRIS_absolute_2500MHz.txt`,
  `TRIS_Beam_Profile.txt`). Point `paths.tris_archive_dir` at your local copy; nothing
  in this repo touches the network.
* **`pygdsm`**, for the GSM2008 template (stage 1's stand-in prior). Stage 2 also
  takes the reprocessed Haslam map (destriped/desourced, Remazeilles et al. 2015) from
  the same astropy cache, so `paths.haslam_map` can stay `null`.
* **`limTOD`** — a normal pip dependency, but if it is not installed the pipeline falls
  back to the sibling source checkout at `paths.limtod` (`../limTOD` by default), which
  is what the walkthrough notebook uses.

---

## Running

The pipeline runs end to end from one command (`PIPELINE.md` Part I describes every step, with a flowchart):

```bash
python run_pipeline.py --step all              # steps 0-7, then paper_outputs/pipeline/PIPELINE_REPORT.md
python run_pipeline.py --step fit              # one step: tris_maps, beam_match, ring_tods, compare, select,
                                               #   fit, predict_408, rhino (each reads the steps before it)
python run_pipeline.py --from-step compare     # this step and every step after it
python run_pipeline.py --step report           # rewrite the report
```

Steps: 0 TRIS maps → 1 beam-match every allowed map to the TRIS beam → 2 ring TODs (120 samples) of every data set
→ 3 model comparison without TRIS (UltraNest) → 4 pick the best model → 5 fit it with TRIS (our spectral model)
→ 6 predict Haslam at 408 MHz through the TRIS beam → 7 predict RHINO (60–80 MHz). Products go to
`paper_outputs/pipeline/`; settings are in `configs/run_pipeline.yaml` (data, beam, ring, `pipeline:`) and `configs/model_comp_config.yaml`
(models, priors, sampler, error bars, the comparison statistic); no setting is in both.
Then run `notebooks/10_paper_plots.ipynb` for the figures.

The earlier interface still works: `python run_pipeline.py --stage 1` (= step 0) and `--stage 2` (the older beam
matching of Haslam and ARCADE 2 alone, which nothing reads). Stage 1 takes two overrides, so a variant run needs no
config edit. Whatever they are set to is what that run's manifest records:

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
* `stage1_manifest.json` — the exact config block used, the solver diagnostics, the
  beam settings and the notes on the GSM2008 prior (the `blocker_*` keys), so a product can always be traced back to the decisions
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

**No convergence cross-check runs at the moment.** Stage 1 can compare each solve against
an exact reference solve (the matrix inversion lemma at a single frequency; the inner matrix
is only 120×120) and refuse to write a product that disagrees. That robustness is not
needed yet, so it is **switched off** (`tris.solver.cross_check: false`, 2026-10-14): the
reference is not computed, nothing can stop the run, and the manifest records the
comparison as `null`. The maps and their posterior σ are limTOD's answer, unverified.
When it was last on, limTOD agreed with the reference to ~1×10⁻¹⁰ σ. Set the key to
`true` to switch it back on (`TODO.md` S1-8).

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
the product `TRIS_MAP_PROVENANCE.md` documents.

The solve is **limTOD's own** (`TRISMapMakingInputs.solve`), on a deliberately coarse
nside 16 grid, over the pixels the beam actually saw. `TRIS_MAP_PROVENANCE.md` has the
full settings and the threshold convergence study.

Those fitted offsets reproduce the GSM2008 deficits measured independently in
`notebooks/diagnostics.ipynb` Part A (1.92 K and 0.56 K cold) — a real check that the geometry, beam and
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
3. smooth it to the TRIS beam with the kernel set by `beam_matching.kernel`, mask-aware
   for ARCADE 2. **In use: `circular_gaussian`**, a 23.366° circular Gaussian
   (`beam_matching.target_fwhm_deg`, the larger TRIS axis) via `hp.smoothing`, with the
   quadrature kernel √(target² − native²): 23.35° for Haslam, 20.28° for ARCADE 2. It
   leaves every map 22% too wide along the TRIS E axis.
   * **Implemented but not in use:** `tris_cuts`, the measured, elliptical TRIS beam
     (19.155° north–south × 23.366° east–west, E plane 7° east of north, the beam the TOD
     route uses), and `elliptical_gaussian`. Both are applied as a direct sum on the sphere
     (`moonrunes/elliptical_beam.py`, `PIPELINE.md` Appendix A.9). They were the default
     on 2026-10-13/14 and were switched back on 2026-10-14.
   * Steps 3–5 are `match_to_tris_beam()`, which the map route also calls;
4. keep a pixel only if it was **observed** and its smoothed weight reached
   `beam_matching.weight_floor` (0.5);
5. regrid onto `beam.nside_new` (= `tris.nside` = 16). A coarse pixel gets a value only if
   all of its sub-pixels were observed (`beam_matching.coarse_pixel_rule`).

```bash
python run_pipeline.py --stage 2
```

Output with the circular kernel (the last run; the elliptical default has not been run yet,
and prints the kernel shape, its E × H widths and why it was chosen instead):

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

## Notebooks

Three notebooks, committed **with their outputs**, so they can be read without re-running. `NOTEBOOK_INDEX.md` lists
what each Part does, what it reads and writes, and where every old notebook went; `PLOT_INDEX.md` lists every plot;
`HANDOVER.md` says which cell makes which paper figure and table.

| Notebook | Contents | Run |
|---|---|---|
| `diagnostics.ipynb` | Checks of the inputs and of stage 1: **A** every input file (TRIS products, ARCADE 2, Haslam, the TRIS maps, every map in `res/`, the TOD each map predicts); **B** coordinate frames, checked against the sky; **C** stage-1 convergence, residuals, the zero level, data against prior, and the prior-only / high-noise runs. | any time after stage 1 |
| `tests_analysis/` (one notebook per Part; `README.md` there) | Every analysis and test behind the paper, in dependency order: **A** beam-matched ring TODs; **B** UltraNest fits per ring sample and the 408 MHz sky models; **C** model-comparison diagnostics; **D** the model-choice decision sheet; **E** T–T indices and the sky models against TRIS; **F** Haslam against the power-law prediction (scale–offset); **G** Haslam with one gain and offset; **H** zero-level offsets; **I** robustness to survey subsets and more; **J** the 13-point fit and the prediction at any frequency, including RHINO; **K** beam tests (FWHM, beam matching, harmonic checks, beam variants, and the script tests below). | after stage 1 and `scripts/model_comparison/run_comparison.py` |
| `10_paper_plots.ipynb` | The paper's figures at MNRAS size, written to `Papers_and_reports/moonrunes_paper/figures/` (PDF + PNG) as `fig1_…`–`fig10_…` (the names `main.tex` includes) and `rhino_tods_tris_beam` (not in the paper). | last |

`diagnostics.ipynb` and the `tests_analysis/` Parts hold the stored outputs of the notebooks they were built from (not re-run). `test_results.ipynb`, the merged notebook the Parts were split from on 2026-11-17, stays as a record and is not maintained.
Those originals are in `../moonrunes_archive_2026-11-16/notebooks/`; six older notebooks are in git at `a5c827c`. **Diagnostics Part C and `tests_analysis` Part K2 predate the switch to nside 16** in places (TODO N-1).

### Where outputs go

**`paper_outputs/`** holds everything a paper figure, table or quoted number is made from; **`test_outputs/`** holds the results of tests and checks the paper does not use; **`outputs/`** keeps the pipeline's stage products; the figures themselves are written to the paper repo. Caches live inside these folders (`cache/` subfolders, `.npz`), so moving a folder moves its cache.

| Folder | What | Written by |
|---|---|---|
| `paper_outputs/pipeline/step1_beam_matched_maps/` | every allowed map beam-matched to the TRIS beam (nside 16) | pipeline step 1 |
| `paper_outputs/pipeline/step2_ring_tods/` | ring TODs and σ of every data set (`ring_tods.npz`); every survey and sky model through the TRIS beam (`beam_matched_tods.npz`, `sky_tods.npz`) | step 2 |
| `paper_outputs/pipeline/step3_model_comparison/` | ln Z per ring sample and model, TRIS excluded | step 3 |
| `paper_outputs/pipeline/step4_model_selection/` | `selection.json`: the best model | step 4 |
| `paper_outputs/pipeline/step5_fit/` | our spectral model: the best model's UltraNest fit at every ring sample, TRIS included | step 5 |
| `paper_outputs/pipeline/step6_predict_408/` | the 408 MHz prediction against Haslam and the 408 MHz models (`tab:408-residuals`, `tab:scale-offset`) | step 6 |
| `paper_outputs/pipeline/step7_rhino/` | the RHINO prediction through the TRIS and RHINO beams (Figure 10) | step 7 |
| `paper_outputs/pipeline/PIPELINE_REPORT.md` | what the last run found, step by step | `--step report` |
| `paper_outputs/model_comparison/diagnostic_plots/` | diagnostics of the step 3 comparison | `notebooks/tests_analysis/` Part C |
| `paper_outputs/model_choice/` | the model-choice decision sheet | Part D |
| `paper_outputs/phil_plots/` | T–T β, sky models against TRIS (`plot6a_models_fitted_to_tris.csv` = `tab:step1`) | Part E |
| `paper_outputs/zero_level_offsets/` | zero-level offsets as fit parameters | Part H |
| `paper_outputs/robustness_check/` | robustness of the Haslam gain and offset | Part I |
| `paper_outputs/model_comparison/results/` | the standalone `run_comparison.py` output (the 2026-10-29 run on 4 points) | `scripts/model_comparison/run_comparison.py` |
| `paper_outputs/nb08_paper_plots/` | the removed notebook 08's figures and tables; `fig15_*.csv` is what Figure 9's cell checks against | notebook 08 (removed) |
| `test_outputs/beam_matched_vs_old/` | step 2's TODs against the old (native-beam) TODs | Part A |
| `test_outputs/sky_models_408/` | curved power-law fits, sky models against both fits (`chi2_sky_models_vs_fits.csv`) | Part B |
| `test_outputs/plot3_haslam_vs_fit/` | the scale-offset fits and figures (a check of step 6) | Part F |
| `test_outputs/plot3_haslam_calibration/` | Haslam with one gain and offset | Part G |
| `test_outputs/tod_prediction/` | 13-point fits and predicted TODs (Tests 1–5); Test 6's last outputs (2026-11-05) | Part J |
| `test_outputs/beam_variants/`, `fwhm_scaling_test/`, `scripts/beam_averaging_blobs/<map>/`, `lwa35_stockert/` | beam tests | Part K2; the beam-test scripts (Parts K3–K5) |
| `outputs/stage1*/` | TRIS maps (the pipeline's step 0 writes `outputs/stage1/`) | step 0 / `--stage 1` |
| `outputs/stage2/` | the older Haslam + ARCADE 2 beam matching (nothing reads it) | `--stage 2` |
| paper repo `figures/` | the paper figures `fig1_…`–`fig10_…` | `10_paper_plots.ipynb` |

---

## Standalone tests

Scripts outside the notebooks, each self-contained. They import the repo's code read-only
and change nothing in `src/`. Run them in the `limtod` conda environment.
Their saved results are shown in `notebooks/tests_analysis/` Parts K3–K5.

| Test | What it checks | Rerun | Outputs |
|---|---|---|---|
| **FWHM scaling** — `scripts/fwhm_scaling_test.py` (2026-10-19) | The two-step vs one-step beam matching of `notebooks/tests_analysis/` Part K2 (`PIPELINE.md` App. E), rerun with four FWHMs given to the pipeline for the complicated survey beam: major axis, minor axis, their arithmetic mean and geometric mean. All three kernels; Haslam 408 MHz, nside 64 → 16. | `python scripts/fwhm_scaling_test.py` (~7.5 min); `--replot-maps` redraws the residual maps from the saved `.npz` in seconds | `test_outputs/fwhm_scaling_test/`: table CSV, comparison plot, residual maps (PNG + `.npz`) |
| **Beam averaging with spectral-index blobs** — `scripts/beam_averaging_blobs/` (2026-10-19/20) | Whether beam averaging changes the *shape* of T vs LST across frequencies when the TRIS beam covers regions of different β. Haslam extrapolated with a blob β map, observed with the TOD route's TRIS beam at 408/600/820 MHz, compared with a uniform-β control. | `python scripts/beam_averaging_blobs/beam_averaging_blobs.py` (seconds); `replot_from_csv.py` redraws the summary plots from the CSVs | `test_outputs/beam_averaging_blobs/<map>/`; details in `scripts/beam_averaging_blobs/README.md` |

**Findings:**

- **FWHM scaling:** telling the pipeline the **major-axis** FWHM gives the smallest
  discrepancy for every kernel. The minor axis gives the largest; the two means fall in
  between.
  - The std of the residual (two-step − one-step) goes from 0.45 to 0.26 K (circular),
    0.50 to 0.30 K (elliptical) and 0.37 to 0.17 K (TRIS cuts), minor axis to major axis.
  - The mean residual is −0.002 to −0.004 K throughout.
- **Beam averaging:** the shape changes by up to 0.27% (600 MHz) and 0.49% (820 MHz), on
  `HASLAM_DESTRIPED_ONLY.fits`; the normalisation changes by less than 0.1%. Without the beam
  the same blobs change the shape by 8% and 15%, so the beam dilutes them strongly. The
  largest effect comes from a blob next to Cygnus A, which the destriped map keeps.

---

## Layout

```
PIPELINE.md                 the whole pipeline in detail, as it currently stands
TODO.md                     open items from stages 1 and 2, ordered by what they block
DATA_SOURCES.md             every external dataset, its archive page, and what reads it
TRIS_MAP_PROVENANCE.md      the limTOD code, settings and conventions behind the stage 1 maps
configs/run_pipeline.yaml     every flagged decision, with the reason next to it
run_pipeline.py             the single entry point (--step NAME | --step all; --stage 1/2 still work)
src/moonrunes/
    config.py               loader; `require()` raises on null = "still an open decision"
    stage1_tris_maps.py     stage 1 + its product loader
    stage2_beam_matching.py stage 2 -- Step 0 beam matching
    frames.py               declared frames; Galactic -> equatorial rotation (full-sky and mask-aware)
    sky_maps.py             reads every map in res/ (any layout), with frequency, frame, beam and CMB convention
    sed_prototype.py        stage 3 prototype: the TOD route and the map route of the SED calibration
    pipeline.py             the end-to-end pipeline, steps 0-7 (run_pipeline.py --step ...)
notebooks/                  diagnostics.ipynb, tests_analysis/ (Parts A–L), 10_paper_plots.ipynb (committed with outputs); test_results.ipynb (record)
scripts/                    standalone tests (fwhm_scaling_test.py); see "Standalone tests"
scripts/beam_averaging_blobs/       the spectral-index blob test: script, CSVs, plots, README
scripts/beam_tests_lwa35_stockert1420/  the beam-matching and blob tests at LWA1 35 and Stockert 1420 MHz: scripts, CSVs, plots, README
paper_outputs/              results that a paper figure, table or quoted number is made from (see "Where outputs go")
test_outputs/               results of tests and checks the paper does not use
outputs/                    pipeline stage products (run_pipeline.py), mostly gitignored
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

73 tests. The algebra ones — the Krylov solve against a dense reference, the exact
reference posterior against a dense inverse, stage 2's beam quadrature and masked-smoothing behaviour, and the frame rotations — run
anywhere. The end-to-end tests skip themselves when the TRIS archive, the FITS downloads,
or `pygdsm` is absent.

---

## License

MIT — see `LICENSE`.
