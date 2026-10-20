# Beam averaging with spectral-index blobs

**Question (from the supervisor):** when regions inside the TRIS beam have different
spectral indices, does beam averaging change the *shape* of antenna temperature vs LST
across frequencies, and not only its normalisation?

Run 2026-10-19/20. The recorded results use `HASLAM_DESTRIPED_ONLY.fits`. An earlier run on
`HASLAM.fits` is kept in `HASLAM/` for comparison.

## What the test does

1. **Load Haslam 408 MHz** with the repo's loader, `moonrunes.sky_maps.describe`
   (nside 512, Galactic, K, CMB included).
2. **Build a spectral-index map** β(pixel) on Haslam's own grid: a background β₀ = −2.7 plus
   Gaussian blobs centred on the TRIS stripe (Dec 42.43°, the site latitude).
   β = β₀ + Σᵢ (βᵢ − β₀) · exp(−r²/2σᵢ²), so βᵢ is the value at a blob's centre.
3. **Extrapolate:** T_ν = T₄₀₈ · (ν/408)^β(pixel) for ν = 408 MHz (reference) and the TRIS
   frequencies in `tris.frequencies_mhz` (600, 820 MHz).
4. **Observe along the ring** with the TOD route's own code,
   `sed_prototype.Inputs._native_tod` / `Inputs.observe`. This is limTOD's measured TRIS cut
   beam, zenith pointing, selfrot −7°, the horizon mask, nside 16, and the 120 ring LSTs
   (LST = RA). The methods are called unchanged on a stand-in object that carries only the
   attributes they read.
5. **Control:** the same steps with β = −2.7 everywhere.
6. **Compare.** Each curve is normalised by its mean over the ring, S = T / mean(T).
   - *Shape change* = 100 · (S_blob / S_control − 1) [%].
   - *Normalisation change* = mean(T_blob) / mean(T_control) − 1.
   - A *no-beam* version samples the unconvolved map at the ring points, to show how much
     the beam dilutes the effect.

The blob parameters are at the top of `beam_averaging_blobs.py`. They were fixed before
the first run and have not been tuned:

| blob | RA = LST | Dec | FWHM | β |
|---|---|---|---|---|
| 1 | 30° | 42.43° | 4° | −2.5 |
| 2 | 120° | 42.43° | 4° | −2.9 |
| 3 | 210° | 42.43° | 4° | −2.6 |
| 4 | 300° | 42.43° | 4° | −2.8 |

## Results

**Control:** with uniform β, the normalised curve is the same at every frequency to
7×10⁻¹⁶ (rounding). The method introduces no shape change of its own.

| ν | largest shape change | normalisation change | largest change without the beam |
|---|---|---|---|
| 600 MHz | 0.27% (LST 297°) | −0.044% | 8.07% |
| 820 MHz | 0.49% (LST 297°) | −0.077% | 15.07% |

- **The shape changes; the normalisation barely does.** Each blob leaves a bump (flatter
  spectrum than β₀) or a dip (steeper) about one beam wide. The effect is larger at 820 MHz,
  which is further from 408 MHz.
- **The beam dilutes the effect strongly:** a 4° blob fills a small part of the roughly
  19° × 23° beam.
- **Blob 4 sits about 1.7° from Cygnus A.** The destriped map keeps point sources, so this
  blob has the largest effect (−0.49% at 820 MHz). It was −0.23% on `HASLAM.fits`, which has
  the sources removed. With bright sources in the map, the result depends on what each blob
  happens to cover.
- On `HASLAM.fits` the largest shape change was 0.23% / 0.41% (600 / 820 MHz), at LST 30°.

## How to rerun

Use the `limtod` conda environment (it has healpy and limTOD). Both scripts take a few
seconds.

```bash
# the analysis: writes beam_averaging_blobs/<map stem>/
python beam_averaging_blobs/beam_averaging_blobs.py

# summary figures from the saved CSVs only (no analysis); RUN defaults to HASLAM_DESTRIPED_ONLY
python beam_averaging_blobs/replot_from_csv.py [RUN]
```

To change the input map, edit `HASLAM_FILE` at the top of `beam_averaging_blobs.py`. Any
full-sky map in `res/` works, and its results go to a new subfolder named after it. The
blob parameters are next to it.

## Files

| File | What it is |
|---|---|
| `beam_averaging_blobs.py` | the analysis (parameters at the top) |
| `replot_from_csv.py` | Plots A and B from the CSVs |
| `<run>/summary.csv` | per frequency: largest and rms shape change, normalisation change, T_ν/T₄₀₈ range, control check, largest no-beam change |
| `<run>/t_vs_lst.csv` | per LST and frequency: T (blob, control) [K], normalised S, shape change [%], T_blob − T_control [K] |
| `<run>/t_vs_lst_comparison.png` | per frequency: normalised T, T_ν/T₄₀₈ and shape change vs LST, blob vs control |
| `<run>/beta_map.png` | the β map in equatorial, with the ring, to check the blob positions |
| `<run>/plotA_shape_change_vs_lst.png` | shape change vs LST at 600 and 820 MHz, blobs marked, largest change annotated |
| `<run>/plotB_beam_dilution.png` | largest shape change without vs with the beam (bars) |

Only the largest no-beam change is saved, not the no-beam curve vs LST, so Plot B uses bars.

## Assumptions and caveats

- **The TRIS beam is the same at 600 and 820 MHz,** as everywhere in the pipeline.
- **Nominal frequencies** (600, 820 MHz) are used, not the rings' effective 600.5 and
  817.8 MHz.
- **The power law is applied to the map as given.** The CMB, the zero level and (in the
  destriped map) point sources are all scaled with β. If the CMB were subtracted first and
  added back, the control would no longer be exactly shape-invariant.
- **β is applied at nside 512, before any degrading,** so temperatures are averaged after the
  power law, not β before it.
- **Normalisation** divides by the mean over the 120 ring samples.
- **This compares the blob sky with a uniform sky.** That mixes the sky change itself with
  the beam's averaging; the no-beam numbers separate the two. It does *not* compare against
  a single beam-averaged β per LST.
- **The effect scales roughly with blob area and with βᵢ − β₀.** Larger or more numerous
  blobs give a proportionally larger effect.
- **The beam is pixelised at about 3.7°** (nside 16, `tris.nside_hires` 16). The integrated
  effect is sound, but the exact bump shapes depend on that pixelisation.
