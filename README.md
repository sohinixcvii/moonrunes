# moonrunes

**A TRIS-anchored spectral model of the low-frequency radio sky, used to test the Haslam 408 MHz map and to predict
the sky seen by RHINO (60–80 MHz).**

TRIS measured the absolute sky brightness at 600 and 820 MHz along a ring at Dec +42.4°. Every other survey is put
through the same beam along the same ring, a spectral model is chosen and fitted at each of the 120 ring samples with
TRIS as the absolute anchor, and the fitted model is compared with the 408 MHz map and with the sky models (GSM2008,
GSM2016, LFSM and their Carter, Irfan and GINAN corrections), and carried to the RHINO band.

| Document | What it covers |
|---|---|
| `QUICKSTART.md` | the handover: how to reproduce the paper, which cell makes which figure, table and number, data restrictions, open items |
| `docs/PIPELINE.md` | every pipeline step and method in detail, with the reasoning and references |
| `docs/DATA_SOURCES.md` | every dataset: source, units, frame, beam, uncertainties |
| `docs/TRIS_MAP_PROVENANCE.md` | the TRIS map-maker (step 0): limTOD settings, conventions, checks |
| `docs/PLOT_INDEX.md` | every plot the notebooks and scripts write |
| `notebooks/tests_analysis/README.md` | the analysis notebooks, Parts A–L |
| `CHANGELOG.md` | the decisions in force, with their dates and effects |

---

## Setup

```bash
conda env create -n limtod -f environment.yml
conda activate limtod
```

- **limTOD** is used from a source checkout next to this repo (`configs/run_pipeline.yaml` `paths.limtod`,
  `../limTOD`; github.com/zzhang0123/limTOD at commit `dbe720b`), which has changes to `limTOD.tris` that the PyPI
  release does not have. Do not `pip install limTOD` (or `pip install -e .`, which pulls it in): an installed copy is
  imported first. The code puts `src/` and the checkout on `sys.path` itself, so nothing needs installing.
- **The TRIS archive** (four LAMBDA text files) is read from `paths.tris_archive_dir` (`../limTOD/downloads/TRIS`).
- **pygdsm** downloads its sky-model data on first use.
- **The RHINO beam and the Jodrell Bank horizon mask** (step 7 only) are read from `../rhino-beams/`, outside the repo
  (`pipeline.rhino.beam_file`, `horizon_file`).
- **Irfan & Puglisi 2025's β_s and c_s maps** (steps 2 and 7) go in `external/irfan_puglisi_2025/` (git-ignored; licence
  not cleared, see `QUICKSTART.md`).
- **CHIPASS** is too large for GitHub; download it as `docs/DATA_SOURCES.md` describes.

## Running

```bash
python run_pipeline.py --step all              # steps 0–7, then paper_outputs/pipeline/PIPELINE_REPORT.md
python run_pipeline.py --step rhino            # one step (it reads the products of the steps before it)
python run_pipeline.py --from-step compare     # this step and every step after it
python run_pipeline.py --step report           # rewrite the report from the steps already run
```

| Step | Name | What it does |
|---|---|---|
| 0 | `tris_maps` | Wiener-filter maps of the TRIS rings with limTOD (`outputs/stage1/`) |
| 1 | `beam_match` | every allowed map beam-matched to the TRIS beam on the nside 16 grid |
| 2 | `ring_tods` | every map, Haslam and the sky models observed along the ring (120 samples) through the TRIS beam |
| 3 | `compare` | model comparison per ring sample, TRIS excluded (UltraNest evidence) |
| 4 | `select` | the best model |
| 5 | `fit` | the best model fitted with TRIS: our spectral model |
| 6 | `predict_408` | our model at 408 MHz against Haslam and the 408 MHz sky models |
| 7 | `rhino` | our model at 60–80 MHz through the TRIS beam and, after Wiener deconvolution, the RHINO beam |

Steps 1–7 write `paper_outputs/pipeline/step*/`. All settings are in `configs/run_pipeline.yaml` (data, beam, ring,
`pipeline:`) and `configs/model_comp_config.yaml` (models, priors, sampler, error bars, comparison statistic); no
setting is in both. A full run takes one to two hours (step 2 alone about 20 minutes; steps 3 and 5 run UltraNest at
every sample). The step products (`.npz`) are git-ignored, so a fresh clone runs the pipeline first.

Then run the notebooks in `notebooks/tests_analysis/` (Parts A–L) and finally `notebooks/10_paper_plots.ipynb`, which
writes the paper's figures. `QUICKSTART.md` §1 gives the order.

`python run_pipeline.py --stage 1` runs step 0 alone and accepts `--nside` and `--run-name` for variant maps (give each
variant its own `--output-dir`: the manifest name is fixed). `--stage 2` beam-matches Haslam and ARCADE 2 on their own
into `outputs/stage2/`; the pipeline does not read that product.

## Layout

```
run_pipeline.py             the entry point
configs/                    run_pipeline.yaml, model_comp_config.yaml, error_bars.csv
src/moonrunes/
    pipeline.py             steps 0–7, their products and the report
    stage1_tris_maps.py     the TRIS map-maker (step 0) and its product loader
    stage2_beam_matching.py beam matching of a map to the TRIS beam, frame and grid (step 1 and --stage 2)
    elliptical_beam.py      the TRIS beam from the archive cuts; beam matching on the sphere
    sed_prototype.py        the TOD-route inputs: ring TODs, survey uncertainties, coverage
    sky_maps.py             reads every map in res/ with its frequency, frame, beam and temperature convention
    frames.py               declared frames; Galactic -> equatorial rotation (full-sky and mask-aware)
    config.py               the config loader (a null setting is an open decision and raises)
    style/                  the plotting style and colours
scripts/
    model_comparison/       models.py (the spectral models), run_comparison.py (standalone comparison, TOD and map routes)
    tests/                  pytest
    fwhm_scaling_test.py, beam_averaging_blobs/, beam_tests_lwa35_stockert1420/   standalone beam tests
notebooks/
    tests_analysis/         Parts A–L: the analyses and checks around the pipeline
    10_paper_plots.ipynb    the paper's figures
    diagnostics.ipynb       checks of the inputs, the coordinate frames and the TRIS maps
res/                        input maps and the TRIS archive files
external/                   third-party maps that may not be redistributed (git-ignored)
paper_outputs/              everything a paper figure, table or quoted number is made from
test_outputs/               results of tests and checks the paper does not use
outputs/                    step 0 / stage products (mostly git-ignored)
docs/                       the documentation listed above, and the early figure plans (PDF)
```

## Tests

```bash
pytest
```

77 tests (`scripts/tests/`). The algebra tests (solver, beam, frame rotations, the Irfan split) run anywhere; the
end-to-end tests skip themselves when the TRIS archive, the FITS files or `pygdsm` are absent.

## Licence

MIT (`LICENSE`). The licence covers the code only, not the third-party data in `res/` and `external/`.
