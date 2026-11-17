# Spectral model comparison (skeleton)

> The model comparison of the paper is **pipeline step 3** (`python run_pipeline.py --step compare`, `moonrunes.pipeline.step_compare`): every allowed map, TRIS excluded, UltraNest evidence, written to `paper_outputs/pipeline/step3_model_comparison/` in the same CSV format. It uses `configs/model_comp_config.yaml` (models, priors, sampler, `configs/error_bars.csv`, `min_points`, `comparison.statistic`, `comparison.exclude_tris`) and the helpers in `run_comparison.py` and `models.py`. `run_comparison.py` itself is a standalone tool (the same two configs and data sets, plus the map route), writing to `paper_outputs/model_comparison/results/`.

A Bayesian comparison of spectral models with UltraNest, on the existing SED fit's
data:

- power law
- curved power law
- moment expansion

Written 2026-10-20 as a skeleton. The structure is in place, but the choices below are
first guesses, not decisions.

## Files

| File | What it is |
|---|---|
| `../configs/model_comp_config.yaml` | the model choices: data route, zero-level switch, models, **all priors in one block**, sampler, `comparison.statistic` and `comparison.exclude_tris`. The data (which maps, the Haslam target, beam, ring) are in `../configs/run_pipeline.yaml`; no setting is in both |
| `../configs/error_bars.csv` | the error bar for each input, with an `is_placeholder` column. Kinds: `tris` (the code's own), `fraction` (`value` × \|T\|), `absolute_k` (`value` K), `published` (`scale_frac` × \|T\| ⊕ `zero_level_k` ⊕ `noise_k`) |
| `models.py` | the models (one moment expansion per order in `models.moment_orders`) and the prior transform |
| `run_comparison.py` | loads the data and compares the models per unit, on the statistic set by `comparison.statistic`: the best-fit ln L (minimum χ², no Occam penalty) or UltraNest's ln Z (`evidence`, the code's default and the configured choice) |
| `../paper_outputs/model_comparison/results/` (written by `run_comparison.py`) | `maxlike_per_unit_<route>.csv` and `maxlike_summary_<route>.csv` (maximum likelihood), or `evidence_per_unit_<route>.csv` and `evidence_summary_<route>.csv` (evidence). Each run overwrites only its own statistic's pair |
| `archive/` | frozen records of past runs with their inputs, code, notebook 11 diagnostics and a write-up; `archive/2026-10-27_beam-matched_4-point_run/` is the run of 2026-10-27 |

## Running

```bash
python model_comparison/run_comparison.py [--config configs/model_comp_config.yaml]
```

Building the inputs takes about a minute. The fits take about 6 s per unit for the first three
models, so the TOD route's 120 ring samples take about 13 minutes (more with order 3). The map route has far
more units.

## What it does

1. **Builds the SED data with the existing loader**, `moonrunes.sed_prototype.Inputs`.
   - **TOD route** (default): one SED per TRIS ring sample, made of the ring data plus every
     fit map observed through the TRIS beam. These are the same points as `tod_route`.
     The map TODs are beam-matched (`sed.tod_beam: beam_matched` in
     `configs/run_pipeline.yaml`; `sed_prototype.Inputs.ring_tod`, the same TODs as `notebooks/test_results.ipynb` Part A).
     `native` switches back to the older limTOD TODs. The run header prints which is used, and
     `evidence_summary_tod.csv` records it in the `tod_beam` column.
   - **Which data sets:** the pipeline's: every map `sed.exclude_maps` in
     `configs/run_pipeline.yaml` leaves in, plus TRIS 600 and 820 unless `comparison.exclude_tris` (true: no TRIS,
     as pipeline step 3). The run header prints the data sets used, and both summary CSVs record them in a
     `datasets` column. A map with no data at a ring sample (CHIPASS and Rhodes do not reach Dec 42.4°) adds no
     point there. With the current configs each ring sample has **11 points**: LWA1 35, LWA1 38–74 ×7, Maipu 45,
     LWA1 80 and Stockert 1420. The stored results of 2026-10-29 used 4 points; the 2026-10-27 4-point run is in
     `archive/2026-10-27_beam-matched_4-point_run/`.
   - **Map route:** one SED per TRIS-stripe pixel, made of the stage 1 maps plus every
     beam-matched map. These are the same points as `map_route`.
2. **Replaces the non-TRIS error bars with the CSV's values.** TRIS keeps the code's own
   errors: `ring_sigma` for ring data, the stage 1 posterior for maps.
3. **Keeps a point only if T > 0, σ > 0 and both are finite,** the rule
   `sed_prototype.fit_power_law` uses. A unit needs at least `min_points` (4) usable points.
4. **Compares the models per unit** with a Gaussian likelihood in T,
   ln L = −χ²/2 − Σ ln(σ√2π). `comparison.statistic` in `configs/model_comp_config.yaml` picks how:
   - **`max_likelihood`: no Occam penalty.** Each model's best fit inside
     the prior box (minimum χ², L-BFGS-B in the prior's unit cube from `comparison.ml_starts` starting
     points, the first at the box centre). Models are compared on ΔlnL_max = −Δχ²/2 against the power
     law, with the same ±1 band as before for "indistinguishable". The priors only set the box. Because
     nothing charges for extra parameters, a model that contains another can never fit worse: the
     curved power law (C = 0 is the power law) always has ΔlnL_max ≥ 0. `maxlike_per_unit_<route>.csv`
     has, per unit and model, χ², ln L_max, degrees of freedom and the best-fit parameters;
     `maxlike_summary_<route>.csv` has per-model medians and counts, and no sums over units.
   - **`evidence`:** UltraNest's ln Z per unit and model, which includes the Occam penalty (the cost of
     each extra parameter is about ln(posterior width / prior width)). A model's total ln Z is the sum
     over units.
5. **Prints a prominent warning,** before and after the results, if any point used a
   placeholder error bar.

**Zero-level switch:** `zero_level.add_zhat_to_tris` adds stage 1's fitted ẑ to the TRIS
temperatures (TODO N-3). It affects the **map route only**. The TOD route fits the archive
ring data, which contain no ẑ.

## Models and priors

With x = ln(ν / 408 MHz) (`models.py`):

- **Power law:** T = A·e^(βx)
- **Curved power law:** T = A·e^(βx + Cx²)
- **Moment expansion (order 2):** T = A·e^(β̄x)·(1 + ω₁x + ω₂x²/2)
- **Moment expansion (order 3):** T = A·e^(β̄x)·(1 + ω₁x + ω₂x²/2 + ω₃x³/6), added 2026-10-27

**Which moment orders are compared** is set by `models.moment_orders` in `configs/model_comp_config.yaml`, now `[2, 3]`: each order is
its own model, `moment expansion (order N)`, with its own ln Z columns in `evidence_per_unit_*.csv` and its own row
in `evidence_summary_*.csv`. Both orders use the same `moment_beta_fixed` (−2.7), so order 3 has four parameters
(A, ω₁, ω₂, ω₃): 7 degrees of freedom on the current 11-point ring samples (none on the 4-point samples of the
2026-10-29 run, where its evidence mostly reflected its prior). The older single `moment_order: N` is still read if `moment_orders` is absent. `test_results.ipynb` Part C
diagnoses both orders; Part D uses only the order-2 model.

**The moment expansion's β is fixed at −2.7.** This is set by `models.moment_beta_fixed` in
`configs/model_comp_config.yaml`. With a number there, β̄ is fixed at that value and ω₁ is free (parameters A, ω₁,
ω₂). With `null`, β is free and ω₁ = 0 (parameters A, β, ω₂), the setup before 2026-10-22.
To test another value, change that one line and rerun. The code that uses it is
`build_models` in `models.py`. The value is printed in the run header and written to
`evidence_summary_*.csv`. `test_results.ipynb` Part C caches these fits under a `beta<value>` tag, so
different values don't overwrite each other.

Every prior is independent and flat between two limits, set in `configs/model_comp_config.yaml` (`priors:`).
None was tuned to the data.

| parameter | used by | prior | range |
|---|---|---|---|
| `A_K`: CMB-removed RJ temperature at 408 MHz | all three | flat in log(A) | 1 – 1000 K |
| `beta`: spectral index | power law, curved PL (moment expansion only if `moment_beta_fixed: null`) | flat | −3.5 to −1.5 |
| `C`: curvature | curved power law | flat | −0.5 to +0.5 |
| `omega1`: 1st-order moment | moment expansion, only if `moment_beta_fixed` is set (now: −2.7) | flat | −1 to +1 |
| `omega2`: 2nd-order moment | moment expansion | flat | −1 to +1 |
| `omega3`: 3rd-order moment | moment expansion (order 3), when 3 is in `moment_orders` (on) | flat | −1 to +1 |

- **Matched curvature ranges:** the moment expansion's x² coefficient is ω₂/2, which spans
  −0.5 to +0.5, the same as C.
- **Prior width changes ln Z:** each extra parameter costs roughly
  ln(prior width / posterior width). A wider C range penalises the curved model more, a
  narrower one less. Part 3 of `notebooks/test_results.ipynb` Part C tests
  C ranges of ±0.5, ±0.2 and ±0.1.
- **The moment model can go negative inside its prior:** at LWA 35 MHz x² ≈ 6, so
  1 + ω₂x²/2 < 0 for ω₂ below about −0.33 (with ω₁ free the bracket is 1 + ω₁x + ω₂x²/2, and
  can still go negative). Those parameter sets get a very poor likelihood
  rather than breaking the fit, but part of the ω₂ prior is effectively wasted, which lowers
  the moment model's ln Z.
- **No zero-level offset parameter** is fitted. The `add_zhat_to_tris` switch only affects
  the map route.

## Assumptions and caveats

- **Units are treated as independent, so ln Z is summed.** Neighbouring ring samples are
  3° apart and the beam is about 20° wide, so they are strongly correlated. The sum
  overstates the evidence, in both directions.
- **The evidence depends on the prior ranges,** especially for the extra curvature
  parameter. The curved power law's C and the moment expansion's ω₂/2 are given the same
  range so that neither is favoured by its prior.
- **Curved power law and moment expansion are almost the same model to second order:**
  exp(C x²) ≈ 1 + C x². Expect similar evidences.
- **The moment expansion's slope comes from ω₁** while β is fixed (currently −2.7). With
  `moment_beta_fixed: null` it instead fixes ω₁ = 0 and fits β.
- **The pivot is 408 MHz** and A is the CMB-removed Rayleigh-Jeans temperature there, as in
  `sed_prototype`.
- **Survey errors:** every non-TRIS map uses its published scale,
  zero-level and noise errors in quadrature (kind `published` in `error_bars.csv`, the same
  values as `sed_prototype.PUBLISHED_ERRORS`; sources in `../DATA_SOURCES.md`). No point is
  marked as a placeholder. Results in `paper_outputs/model_comparison/results/` from 2026-10-20 used a 10% placeholder.
- **Errors are diagonal.** No inter-survey correlation, scale (gain) term or zero-level
  term is included for the surveys.
- **The TRIS points use the effective frequencies** (600.5 and 817.8 MHz), as
  `sed_prototype` does.
- **The sampler settings are UltraNest's defaults** (400 live points, dlogz 0.5).
