# moonrunes — plot index

As of 2026-10-15. Online copy: <https://claude.ai/code/artifact/a5a7762b-b179-4813-8aa6-d3ea8d198d39>

## How to use this index

Every plot in the project lives in one of three places. Saved files are listed by name; plots that exist only inside a notebook are listed by notebook, cell number and section heading.

| Where | What is there | Made by |
| --- | --- | --- |
| `moonrunes_paper/figures/` (in `Papers_and_reports`) | the 5 paper figures, PDF + PNG, with captions in its `README.md` | `notebooks/10_paper_plots.ipynb` |
| `moonrunes/outputs/paper_plots/` | the working figure set (Figures 1–13b), PDF + PNG | `notebooks/08_update_plots.ipynb` |
| inside each notebook | diagnostic and exploration plots, not saved to files | notebooks 01–09 |
| `moonrunes/outputs/fwhm_scaling_test/` | the FWHM scaling test (added 2026-10-19) | `scripts/fwhm_scaling_test.py` |
| `moonrunes/beam_averaging_blobs/<map>/` | the spectral-index blob test (added 2026-10-19/20) | `beam_averaging_blobs/*.py` |

Cell numbers are 0-based, as `nbformat` counts them, and were taken on 2026-10-15. They shift if cells are added above. The section heading is the stable way to find a plot.

All SED results here exclude ARCADE 2 (`sed.exclude_maps`) and use the 23.366° circular beam-matching kernel (`beam_matching.kernel: circular_gaussian`).

## Paper figures — `moonrunes_paper/figures/`

Five figures, made by `notebooks/10_paper_plots.ipynb`, drawn at the MNRAS full width (174 mm). Full captions are in `figures/README.md`.

| File | Shows | Notebook 10 cell |
| --- | --- | --- |
| `tris_maps.pdf` | TRIS 600 and 820 MHz maps from stage 1 (limTOD), and the six 2.5 GHz samples vs RA | 3 |
| `beam_matching_stockert.pdf` | Stockert + Villa-Elisa 1420 MHz before and after matching to the 23.366° TRIS beam | 5 |
| `sed_fit_tod_route.pdf` | SED fits at six TRIS ring samples (TOD route), with residuals | 9 |
| `sed_fit_map_route.pdf` | SED fits at six TRIS-stripe pixels (map route), with residuals | 11 |
| `tods_all_maps.pdf` | every survey observed through the TRIS beam with limTOD, vs RA, with the TRIS ring data | 13 |

The README also holds three empty placeholder sections for planned figures: best-fit 408 MHz TOD, reforming the 408 MHz map, and extending to the RHINO band.

## Notebook 08 — update plots (`outputs/paper_plots/`)

`notebooks/08_update_plots.ipynb`, target = destriped-only Haslam. Every figure is saved as `<file>.pdf` and `.png`.

| Fig | File | Shows | Cell |
| --- | --- | --- | --- |
| 1 | `fig01_maps_used` | every map in the SED fit, native resolution, with the TRIS stripe | 4 |
| 2 | `fig02_target_stripe` | the destriped Haslam target in equatorial, with the stripe | 6 |
| 3 | `fig03_stripe_slices` | RA/Dec slices of every map along the TRIS stripe | 8 |
| 4 | `fig04_tris_maps_and_tods` | the stage 1 TRIS maps and the ring TODs | 10 |
| 5 | `fig05_sed_fits` | 12 pixel-by-pixel SED fits (map route) with residuals | 12 |
| 5b | `fig05b_sed_fits_tod` | 12 sample-by-sample SED fits (TOD route) with residuals | 14 |
| 6 | `fig06_fit_residual_maps` | best fit − map in K, every survey, over the stripe | 16 |
| 6b | `fig06b_tod_residuals` | best fit − data along the ring, every survey | 18 |
| 7 | `fig07_bestfit_haslam_both_routes` | best-fit Haslam from both routes, and each minus old Haslam | 20 |
| 7b | `fig07b_tod_route_ring` | old Haslam TOD vs the best-fit 408 MHz prediction along the ring | 22 |
| 7c | `fig07c_fractional_difference_hist` | histogram of (best fit − old)/old, both routes | 24 |
| 7d | `fig07d_scaled_tods_vs_lst` | every map's TOD scaled to 408 MHz by a power law fitted without TRIS, vs LST, with GSM2008, GSM2016 and LFSM; dashed lines mark the Haslam survey-zone boundaries along the ring | 26 |
| 8 | `fig08_residuals_by_survey` | pixel-by-pixel residuals (data − fit)/σ, per survey | 29 |
| 9 | `fig09a_chi2_maps`, `fig09b_chi2_tod` | χ²/dof maps (map route) and along the ring (TOD route) | 31 |
| 10 | `fig10a_beta_maps`, `fig10b_beta_tod` | β maps (map route) and along the ring (TOD route) | 33 |
| 11 | `fig11_corrected_full_resolution` | corrected Haslam at full resolution, both routes | 35 |
| 12 | `fig12a_tt_raw`, `fig12b_tt_beam_convolved` | T–T corner plots, raw and beam-convolved | 37 |
| 13 | `fig13_tt_beta` | β from every T–T pair, raw and beam-convolved, with bootstrap errors | 39 |
| 13b | `fig13b_beta_matrix_tt_vs_sed` | β matrix: T–T (below the diagonal) vs map-route SED fit (above), and their difference | 41 |

Also in the folder:

- `haslam_zones_tris_ring_check.png`: the TRIS ring traced over Fig. 1 of Nasirudin & Bull (2025), coloured by the Haslam zone it is read as. It is the check behind the zone lines in Figure 7d. Not made by a notebook cell.
- `fig07a_old_vs_bestfit_map_route` and `fig07b_old_vs_bestfit_tod_route` are leftovers from an earlier version of Figure 7. No current cell writes them.

## Notebook 07 — SED calibration

`notebooks/07_sed_calibration.ipynb`, target = the standard Haslam map. Plots stay in the notebook.

| Cell | Section | Shows |
| --- | --- | --- |
| 10 | 3. Calibration stage 1 — the TOD route | Haslam TOD vs predicted 408 MHz; prediction − Haslam; β and χ²/dof vs RA |
| 12 | β vs LST with different points left out | β along the ring with each survey left out, and with groups left out (no TRIS, no surveys below 100 MHz) |
| 18 | 5. The maps | 3 figures: old and TOD-route corrected Haslam; map-route corrected, correction and hybrid maps; map-route β, χ²/dof and points per fit |
| 20 | 6. The two routes against each other | map-route vs TOD-route correction, per stripe pixel |

## Notebook 03 — diagnostics

`notebooks/03_diagnostics.ipynb` has five parts. Plots stay in the notebook.

| Cell | Part / section | Shows |
| --- | --- | --- |
| 5 | Stage 1 · 2. Does the model describe the ring? | ring data vs model, residuals |
| 7 | Stage 1 · the uncertainties the solver uses | the per-sample error bars |
| 9 | Stage 1 · 3. The zero level | fitted zero level against the template deficit |
| 12 | Stage 1 · 4. Did the data move the map? | prior → posterior σ shrinkage, z-scores |
| 14 | Stage 1 · the band, drawn | the solved band maps |
| 16, 17 | Stage 1 · 5. Is this fit to hand to stage 3? | β between the solved maps; the degrade to nside 16 |
| 24, 26 | Stage 1 · prior-only runs | what the data moved: histograms and sky maps |
| 30, 32, 34 | Stage 1 · high-noise run | original vs high-noise run; per-pixel σ; per-pixel RMS error |
| 37 | Stage 1 · residuals with error bars | residuals with error bars, and the RMS maps |
| 42 | Elliptical beam · 1 | the three kernels on the sky |
| 44 | Elliptical beam · 2 | kernel profiles along E and H, and weight in the shoulders |
| 47 | Elliptical beam · 3 | which kernel each survey gets |
| 49 | Elliptical beam · 4 | a point source through the smoothing code; broad native beams |
| 52 | Elliptical beam · 5 | Haslam with the old vs new kernel |
| 55 | Elliptical beam · 6 | ARCADE 2's mask edge under both kernels |
| 58 | Elliptical beam · 7 | how the new kernel changes every map in the stripe |
| 64 | Harmonic checks · H1 | C_ℓ of every map, raw and after smoothing |
| 67 | Harmonic checks · H2 | the TRIS beam's C_ℓ at nside 16–128 |
| 69 | Harmonic checks · H3 | limTOD's pointed beam vs the exact beam |
| 71 | Harmonic checks · H4 | TODs: harmonic vs exact, at nside 16, 64, 128 |
| 73 | Harmonic checks · H5 | Haslam rotation round trip: ringing at Cas A |
| 75 | Harmonic checks · H6 | harmonic vs real-space Gaussian smoothing |
| 78 | Convolve → deconvolve | the two pointed beams and the operators' singular values |
| 80 | Convolve → deconvolve | Haslam: input, convolved, deconvolved, residual (both beams) |
| 83 | Convolve → deconvolve · the two beams | elliptical − circular, convolved and deconvolved |
| 85 | Convolve → deconvolve · which scales come back | transfer function, correlation by ℓ, error vs smoothing scale |
| 87 | Convolve → deconvolve · every map | rms error per survey, both beams |
| 89 | Convolve → deconvolve · along the TRIS ring | map / input − 1 vs LST, pipeline noise, both beams |
| 93 | Two-step vs one-step | the complicated survey beam: image and profiles |
| 94 | Two-step vs one-step | difference maps: survey beam × kernel |
| 96 | Two-step vs one-step | the difference along the TRIS ring, and a summary |

## Notebooks 01, 02, 04, 05, 06, 09

Plots stay in the notebooks.

| Notebook | Cell | Section | Shows |
| --- | --- | --- | --- |
| 01 explore data | 5 | the four TRIS products | the 600 / 820 MHz rings with errors |
| 01 | 9 | the temperature convention | β between the rings, with and without the CMB |
| 01 | 11 | the beam | archive E/H cuts vs a Gaussian |
| 01 | 15 | the prior template | GSM2008 through the TRIS beam vs the rings |
| 01 | 21, 24 | ARCADE 2 3.15 / 3.41 GHz | the maps and their masks |
| 01 | 30 | Haslam 408 MHz | the map |
| 01 | 37 | TRIS maps | the stage 1 maps |
| 01 | 47, 49, 51 | the annulus; nside 8 vs 16 | the TRIS footprint and its declination profile |
| 01 | 55, 56 | every map in `res/` | all maps side by side; slices along the TRIS stripe |
| 01 | 61 | the TRIS TOD each map predicts | one TOD panel per map, with coverage |
| 01 | 63, 65 | TODs against each other | normalised by the mean; scaled with β = −2.7 |
| 01 | 67 | spectral index against RA | β of each map vs TRIS 600 along the ring |
| 02 Berkhuijsen prior (superseded) | 4–13 | Haslam; region operators; prior circularity; independent prior | 6 figures |
| 04 coordinate frames | 5 | TRIS ring peak | does the ring peak at the plane crossing? |
| 04 | 14 | sanity check | known sources against literature positions |
| 04 | 17 | ARCADE 2 mask | does the zero mask survive rotation? |
| 05 beam matching (circular kernel) | 9, 10 | before and after | Haslam and ARCADE 2 smoothed to 23.366° |
| 05 | 17 | usable overlap | the four-way overlap with TRIS |
| 06 SED fit (TRIS + ARCADE 2 cross-check) | 10 | sanity checks | the four-point SEDs |
| 06 | 12 | full maps | best-fit Haslam, β and χ² |
| 09 beam FWHM tests | 11 | 4. archive-cut beam FWHM | cut beam vs Gaussian, in dB |
| 09 | 15 | 6. Test 2: round trip | FWHM error vs nside |

## Standalone tests (saved files)

Both are described, with rerun commands, under "Standalone tests" in `README.md`.

**FWHM scaling** — `outputs/fwhm_scaling_test/`, made by `scripts/fwhm_scaling_test.py`

| File | Shows |
| --- | --- |
| `fwhm_scaling_comparison.png` | two-step / one-step − 1 [%]: median, std (and std in the TRIS stripe), max \|Δ\|, for each FWHM scaling, one line per kernel |
| `fwhm_scaling_residual_maps.png` | the residual two-step − one-step [K] on the sky, 3 kernels × 4 scalings, one shared scale |

**Beam averaging with spectral-index blobs** — `beam_averaging_blobs/<map>/` (`HASLAM_DESTRIPED_ONLY/` is the recorded run; `HASLAM/` the earlier one)

| File | Shows | Made by |
| --- | --- | --- |
| `plotA_shape_change_vs_lst.png` | change in shape of T vs LST [%] at 600 and 820 MHz, blobs marked, largest change annotated | `replot_from_csv.py` |
| `plotB_beam_dilution.png` | largest change in shape without vs with the TRIS beam (bars, separate axes) | `replot_from_csv.py` |
| `t_vs_lst_comparison.png` | per frequency: normalised T, T_ν/T₄₀₈ and shape change vs LST, blobs vs uniform-β control | `beam_averaging_blobs.py` |
| `beta_map.png` | the β map in equatorial with the TRIS ring | `beam_averaging_blobs.py` |

## Find a plot by topic

| Question | Where to look |
| --- | --- |
| How low is Haslam against TRIS? | 08 Fig 7b; paper `sed_fit_tod_route`; 07 cell 10 |
| What do the SED fits look like? | paper `sed_fit_tod_route`, `sed_fit_map_route`; 08 Figs 5, 5b |
| Which surveys disagree with the power law? | 08 Figs 6, 6b, 8; 08 Fig 13b |
| How does β vary along the ring / over the sky? | 08 Fig 10; 07 cell 12 (β with surveys left out) |
| Is one power law a good fit? | 08 Fig 9 (χ²); 07 cell 12 |
| How do TRIS, Haslam and the sky models compare? | 08 Fig 7d |
| Where are the Haslam survey zones along the ring? | 08 Fig 7d (dashed lines); `outputs/paper_plots/haslam_zones_tris_ring_check.png` |
| β from T–T plots, and vs the SED fit | 08 Figs 12, 13, 13b |
| The corrected Haslam map | 08 Figs 7, 11; 07 cell 18 |
| What the TRIS maps look like | paper `tris_maps`; 08 Fig 4; 01 cell 37 |
| What beam matching does to a map | paper `beam_matching_stockert`; 05 cells 9–10 |
| Every survey through the TRIS beam | paper `tods_all_maps`; 01 cells 61–67 |
| Is the beam width / FWHM right? | 09 cells 11, 15 |
| Elliptical vs circular kernel | 03 cells 42–58 |
| Harmonic-space ringing, C_ℓ | 03 cells 64–75 |
| Convolve then deconvolve | 03 cells 78–89 |
| Two-step vs one-step beam matching | 03 cells 93–96 |
| Which FWHM to give the pipeline for an elliptical survey beam | `outputs/fwhm_scaling_test/` |
| Does beam averaging change the shape of T vs LST with varying β? | `beam_averaging_blobs/HASLAM_DESTRIPED_ONLY/plotA_*`, `plotB_*` |
| Stage 1 fit quality and residuals | 03 cells 5–37 |
| Coordinate frames | 04 cells 5–17 |
