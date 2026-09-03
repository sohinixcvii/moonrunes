# moonrunes

**TRIS map-making (`limTOD`) feeding a Gibbs recalibration of the Haslam 408 MHz map
(`bayesian_skymap`).**

The idea that ties the two codebases together is one sentence: `bayesian_skymap` needs a
high-resolution but badly-calibrated map (Haslam) *and* a low-resolution but **absolutely
calibrated** one — and TRIS is an absolute radiometer. So `limTOD.tris` makes maps from
the TRIS drift rings, and those maps become `bayesian_skymap`'s `data1`/`sky_lowres`.

`tris_haslam_pipeline.md` is the algorithm, stage by stage, and is the document this repo
implements. `configs/run_config.yaml` is the single source of truth for every decision it
flags as "decide rather than default".

---

## Status

| Stage | What it does | State |
|---|---|---|
| **1** | Calibrated TRIS maps from the public archive (`limTOD.tris`) | **implemented** |
| 2 | Prepare Haslam, the gain operator, the Berkhuijsen prior | not started |
| 3 | Assemble the `bayesian_skymap` `.npz` | not started |
| 4 | Run the Gibbs sampler | not started |
| 5 | Validate | not started |

---

## Setup

```bash
conda env create -f environment.yml     # or use an existing env with the same stack
conda activate tris-haslam-cal
pip install -e .
git submodule update --init             # external/bayesian_skymap
```

Stage 1 additionally needs:

* **the TRIS archive** — the four public LAMBDA text files
  (`TRIS_absolute_600.txt`, `TRIS_absolute_820.txt`, `TRIS_absolute_2500MHz.txt`,
  `TRIS_Beam_Profile.txt`). Point `paths.tris_archive_dir` at your local copy; nothing
  in this repo touches the network.
* **`pygdsm`**, for the GSM2008 template (Blocker 2's stand-in prior).
* **`limTOD`** — a normal pip dependency, but if it is not installed the pipeline falls
  back to the sibling source checkout at `paths.limtod` (`../limTOD` by default), which
  is what the walkthrough notebook uses.

---

## Running

```bash
python run_pipeline.py --stage 1
```

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
600 MHz (νeff 600.5): gmres info=0, 114 iters, 1.6e-07 σ from exact
        reduced χ² 33.6   zero level +1.7361 ± 0.0108 K
820 MHz (νeff 817.8): gmres info=0,  72 iters, 1.8e-09 σ from exact
        reduced χ² 13.3   zero level +0.5148 ± 0.0069 K
```

Those fitted offsets reproduce the GSM2008 deficits measured independently in
`notebooks/01` (1.92 K and 0.56 K cold) — a real check that the geometry, beam and
convention translation came through intact.

**Read them as a measurement of GSM2008, not of the archive's zero point.** The monopole
degeneracy is ≈1, so with a per-pixel prior this tight the nuisance offset is the only
place a template error can go. This is the honest limit of one absolutely-calibrated ring
plus a template that is kelvins off, not a defect in the map-maker.

---

## Notebooks

Both are committed **with their outputs**, so they can be read without re-running.

| Notebook | Contents |
|---|---|
| `01_explore_tris_data.ipynb` | The inputs: the four products, the beam, the temperature convention, the GSM2008 template — and where `zero_level_sigma_k`, `uncertainty_floor_k` and `nside_new` come from. |
| `02_build_berkhijsen_prior.ipynb` | Stage 2 work. Empty. |
| `03_diagnostics.ipynb` | The stage 1 results: convergence, χ² and residuals, the zero level against the template deficit, prior→posterior shrinkage and z-scores, the band maps, β between the two solved maps, and a preview of the degrade to `nside_new`. |

---

## The three blockers

| | Status |
|---|---|
| **1 — beam width mismatch** | Resolved. The fork accepts arbitrary `beam_deg`; the TRIS E/H mean (21.2605°) is passed through as itself. Stage 1 re-derives `nside_new` from `bayesian_func.nside_for_beam` at run time and refuses to start if the pinned value has drifted. |
| **2 — no "truth" for real data** | **Open.** GSM2008 currently stands in as *both* the prior template and the reference truth, which is circular by construction; the manifest records `prior_equals_truth: true` so nothing downstream can read agreement with `truth_k` as a measurement. The real answer is a Haslam-independent prior from Berkhuijsen (1972) — stage 2 work. |
| **3 — Woodbury fast path broken at low channel count** | Resolved by not using it. `gibbs.use_woodbury` is `false` and stage 1 solves with the direct Krylov path. |

---

## Layout

```
configs/run_config.yaml     every flagged decision, with the reason next to it
run_pipeline.py             the single entry point (--stage N)
src/moonrunes/
    config.py               loader; `require()` raises on null = "still an open decision"
    stage1_tris_maps.py     stage 1 + the product loader
    stage2..stage5          empty
notebooks/                  analysis and plotting, committed with outputs
external/bayesian_skymap    submodule (patched fork)
outputs/                    stage products, gitignored
tests/test_stage_io.py      pytest
```

Two rules the stages rely on: a `null` in the config is an *open decision* and the stage
that needs it raises rather than defaulting; and every stage records the config block it
consumed in its manifest.

---

## Tests

```bash
pytest
```

The algebra tests — the Krylov solve against a dense reference, the Jacobi diagonal
against `bayesian_func.estimate_diag_precond`, the Woodbury posterior against a dense
inverse — run anywhere. The end-to-end tests skip themselves when the TRIS archive or
`pygdsm` is absent.

---

## License

MIT — see `LICENSE`.
