# Spectral model comparison (skeleton)

A Bayesian comparison of three spectral models with UltraNest, on the existing SED fit's
data:

- power law
- curved power law
- moment expansion

Written 2026-10-20 as a skeleton. The structure is in place, but the choices below are
first guesses, not decisions.

## Files

| File | What it is |
|---|---|
| `config.yaml` | every choice: data route, zero-level switch, models, **all priors in one block**, sampler |
| `error_bars.csv` | the error bar for each input, with an `is_placeholder` column. Kinds: `tris` (the code's own), `fraction` (`value` × \|T\|), `absolute_k` (`value` K), `published` (`scale_frac` × \|T\| ⊕ `zero_level_k` ⊕ `noise_k`) |
| `models.py` | the three models and the prior transform |
| `run_comparison.py` | loads the data, runs UltraNest per unit and model, prints and saves ln Z |
| `results/` | `evidence_per_unit_<route>.csv` and `evidence_summary_<route>.csv` |

## Running

```bash
python model_comparison/run_comparison.py [--config model_comparison/config.yaml]
```

Building the inputs takes about a minute. The fits take about 6 s per unit for all three
models, so the TOD route's 120 ring samples take about 13 minutes. The map route has far
more units.

## What it does

1. **Builds the SED data with the existing loader**, `moonrunes.sed_prototype.Inputs`.
   - **TOD route** (default): one SED per TRIS ring sample, made of the ring data plus every
     fit map observed through the TRIS beam. These are the same points as `tod_route`.
   - **Map route:** one SED per TRIS-stripe pixel, made of the stage 1 maps plus every
     beam-matched map. These are the same points as `map_route`.
2. **Replaces the non-TRIS error bars with the CSV's values.** TRIS keeps the code's own
   errors: `ring_sigma` for ring data, the stage 1 posterior for maps.
3. **Keeps a point only if T > 0, σ > 0 and both are finite,** the rule
   `sed_prototype.fit_power_law` uses. A unit needs at least `min_points` (4) usable points.
4. **Runs UltraNest for each unit and model** with a Gaussian likelihood in T. A model's
   total ln Z is the sum over units.
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

**The moment expansion's β is fixed at −2.7.** This is set by `models.moment_beta_fixed` in
`config.yaml`. With a number there, β̄ is fixed at that value and ω₁ is free (parameters A, ω₁,
ω₂). With `null`, β is free and ω₁ = 0 (parameters A, β, ω₂), the setup before 2026-10-22.
To test another value, change that one line and rerun. The code that uses it is
`build_models` in `models.py`. The value is printed in the run header and written to
`evidence_summary_*.csv`. Notebook 11 caches these fits under a `beta<value>` tag, so
different values don't overwrite each other.

Every prior is independent and flat between two limits, set in `config.yaml` (`priors:`).
None was tuned to the data.

| parameter | used by | prior | range |
|---|---|---|---|
| `A_K`: CMB-removed RJ temperature at 408 MHz | all three | flat in log(A) | 1 – 1000 K |
| `beta`: spectral index | power law, curved PL (moment expansion only if `moment_beta_fixed: null`) | flat | −3.5 to −1.5 |
| `C`: curvature | curved power law | flat | −0.5 to +0.5 |
| `omega1`: 1st-order moment | moment expansion, only if `moment_beta_fixed` is set (now: −2.7) | flat | −1 to +1 |
| `omega2`: 2nd-order moment | moment expansion | flat | −1 to +1 |
| `omega3`: 3rd-order moment | moment expansion, only if `moment_order: 3` (off) | flat | −1 to +1 |

- **Matched curvature ranges:** the moment expansion's x² coefficient is ω₂/2, which spans
  −0.5 to +0.5, the same as C.
- **Prior width changes ln Z:** each extra parameter costs roughly
  ln(prior width / posterior width). A wider C range penalises the curved model more, a
  narrower one less. Part 3 of `notebooks/11_model_comparison_diagnostics.ipynb` tests
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
- **Survey errors:** since 2026-10-21 every non-TRIS map uses its published scale,
  zero-level and noise errors in quadrature (kind `published` in `error_bars.csv`, the same
  values as `sed_prototype.PUBLISHED_ERRORS`; sources in `../DATA_SOURCES.md`). No point is
  marked as a placeholder any more. The results in `results/` from 2026-10-20 were made with
  the earlier 10% placeholder.
- **Errors are diagonal.** No inter-survey correlation, scale (gain) term or zero-level
  term is included for the surveys.
- **The TRIS points use the effective frequencies** (600.5 and 817.8 MHz), as
  `sed_prototype` does.
- **The sampler settings are UltraNest's defaults** (400 live points, dlogz 0.5).
