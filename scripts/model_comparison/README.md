# Spectral model comparison

The model comparison of the paper is **pipeline step 3** (`python run_pipeline.py --step compare`,
`moonrunes.pipeline.step_compare`): every allowed map, TRIS excluded, UltraNest evidence per ring sample, written to
`paper_outputs/pipeline/step3_model_comparison/`. Step 3 uses `configs/model_comp_config.yaml` and the helpers here:
`models.py` (the models and their priors) and `run_comparison.py` (data loading, fits, CSV output).

`run_comparison.py` can also be run on its own (the same configs and data sets, plus the map route), writing to
`paper_outputs/model_comparison/results/`.

## Files

| File | What it is |
|---|---|
| `../../configs/model_comp_config.yaml` | data route, zero-level switch, models, **all priors in one block**, sampler, `comparison.statistic`, `comparison.exclude_tris`. The data (which maps, the Haslam target, beam, ring) are in `configs/run_pipeline.yaml`; no setting is in both |
| `../../configs/error_bars.csv` | the uncertainty of each input. Kinds: `tris` (computed by the code), `published` (`scale_frac` × \|T\| ⊕ `zero_level_k` ⊕ `noise_k`); `fraction` and `absolute_k` are supported |
| `models.py` | the models (one moment expansion per order in `models.moment_orders`) and the prior transform |
| `run_comparison.py` | loads the data and compares the models per unit, on `comparison.statistic`: UltraNest's ln Z (`evidence`, the configured choice) or the best-fit ln L (`max_likelihood`, no Occam penalty) |

## Running on its own

```bash
python scripts/model_comparison/run_comparison.py [--config configs/model_comp_config.yaml]
```

Building the inputs takes about a minute; the fits about 6 s per ring sample for three models, so the TOD route's
120 samples take about 13 minutes (more with order 3). The map route has far more units.

## What it does

1. **Builds the SED data** with `moonrunes.sed_prototype.Inputs`.
   - **TOD route** (`data.route: tod`): one SED per TRIS ring sample, made of the ring data plus every allowed map
     observed through the TRIS beam (`sed.tod_beam: beam_matched`). The run header prints the TOD method and the
     summary CSV records it (`tod_beam` column).
   - **Data sets:** every map `sed.exclude_maps` leaves in, plus TRIS 600 and 820 unless `comparison.exclude_tris`
     (true, as step 3). A map with no data at a ring sample (CHIPASS, Rhodes) adds no point. With the current configs
     each ring sample has **11 points**: LWA1 35, 38, 40, 45, 50, 60, 70, 74, 80, Maipu 45 and Stockert 1420. The
     summary CSVs record the data sets (`datasets` column).
   - **Map route** (`data.route: map`): one SED per TRIS-stripe pixel, made of the step 0 maps plus every
     beam-matched map.
2. **Uses the CSV's uncertainties for the surveys**; TRIS keeps the code's own (`ring_sigma` for the rings, the step 0
   posterior for the maps).
3. **Keeps a point only if T > 0, σ > 0 and both are finite.** A unit needs at least `data.min_points` (4) points.
4. **Compares the models per unit** with a Gaussian likelihood in T, ln L = −χ²/2 − Σ ln(σ√2π):
   - **`evidence`:** UltraNest's ln Z per unit and model, which includes the Occam penalty (the cost of each extra
     parameter is about ln(prior width / posterior width)). A model's total is the sum over units.
   - **`max_likelihood`:** each model's best fit inside the prior box (minimum χ², L-BFGS-B in the prior's unit cube
     from `comparison.ml_starts` starting points, the first at the box centre), compared on Δln L_max = −Δχ²/2 with a
     ±1 band for "indistinguishable". A model that contains another can never fit worse (the curved power law with
     C = 0 is the power law).
5. **Warns** before and after the results if any point used a placeholder uncertainty (none does now).

**Zero-level switch:** `zero_level.add_zhat_to_tris` adds step 0's fitted zero level to the TRIS map temperatures.
It affects the **map route only**; the TOD route fits the archive ring data.

## Models and priors

With x = ln(ν / 408 MHz) (`models.py`):

- **Power law:** T = A·e^(βx)
- **Curved power law:** T = A·e^(βx + Cx²)
- **Moment expansion (order 2):** T = A·e^(β̄x)·(1 + ω₁x + ω₂x²/2)
- **Moment expansion (order 3):** T = A·e^(β̄x)·(1 + ω₁x + ω₂x²/2 + ω₃x³/6)

`models.moment_orders: [2, 3]` compares both orders, each its own model. `models.moment_beta_fixed: -2.7` fixes β̄ at
−2.7 with ω₁ free (parameters A, ω₁, ω₂[, ω₃]); with `null`, β is free and ω₁ = 0. Order 3 has four parameters: 7
degrees of freedom on 11 points.

Every prior is independent and flat between two limits (`priors:`); none was tuned to the data.

| parameter | used by | prior | range |
|---|---|---|---|
| `A_K`: CMB-removed RJ temperature at 408 MHz | all | flat in log(A) | 1 – 1000 K |
| `beta`: spectral index | power law, curved power law (moment expansion only if `moment_beta_fixed: null`) | flat | −3.5 to −1.5 |
| `C`: curvature | curved power law | flat | −0.5 to +0.5 |
| `omega1` | moment expansion, when `moment_beta_fixed` is set | flat | −1 to +1 |
| `omega2` | moment expansion | flat | −1 to +1 |
| `omega3` | moment expansion (order 3) | flat | −1 to +1 |

- **Matched curvature ranges:** the moment expansion's x² coefficient is ω₂/2, which spans −0.5 to +0.5, the same as
  C, so neither model is favoured by its prior. To second order the two models are the same (exp(Cx²) ≈ 1 + Cx²).
- **Prior width changes ln Z:** a wider C range penalises the curved model more. `notebooks/tests_analysis/` Part C
  tests C ranges of ±0.5, ±0.2 and ±0.1.
- **The moment model can go negative inside its prior** (at LWA1 35 MHz x² ≈ 6, so 1 + ω₂x²/2 < 0 for ω₂ below about
  −0.33). Those parameter sets get a very poor likelihood; part of the ω₂ prior is effectively unused, which lowers
  the moment model's ln Z.
- **No zero-level offset parameter** is fitted (`notebooks/tests_analysis/` Part H tests one per survey).

## Assumptions and caveats

- **Units are treated as independent, so ln Z is summed.** Neighbouring ring samples are 3° apart under a ~20° beam,
  so the sum counts each independent beam about six times: it ranks the models but overstates the preference.
- **Uncertainties are diagonal.** Each survey's scale and zero-level terms are added to every point as if
  independent, although they are common to all of that survey's points.
- **The pivot is 408 MHz**, and A is the CMB-removed Rayleigh–Jeans temperature there.
- **The TRIS points use the effective frequencies** (600.5 and 817.8 MHz).
- **Sampler:** UltraNest with 400 live points and dlogz 0.5 (`sampler:`).
