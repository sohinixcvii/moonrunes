# Pipeline report

Written 2026-11-20T14:04:32+00:00 by `python run_pipeline.py --step report` from each step's `manifest.json`.

## Step 1 — Beam matching

Written 2026-11-16T19:31:14+00:00 (run_pipeline.py --step beam_match).

```json
{
  "kernel": "tris_cuts",
  "nside": 16,
  "maps": 14,
  "excluded_by_config": [
    "ARCADE2_315.fits",
    "ARCADE2_341.fits",
    "Juno_1248MHz.fits",
    "Juno_600MHz.fits"
  ],
  "per_map": {
    "LWA_35MHz.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "LWA1_38.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "LWA1_40.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "LWA1_45.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "MAIPU_45mhZ.fits": {
      "pixels_kept_in_tris_band": 1200,
      "of": 1272
    },
    "LWA1_50.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "LWA1_60.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "LWA1_70.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "LWA1_74.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "LWA_80MHz.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    },
    "CHIPASS_1400MHz.fits": {
      "pixels_kept_in_tris_band": 382,
      "of": 1272
    },
    "STOCKERT_VILLA_1400MHZ.fits": {
      "pixels_kept_in_tris_band": 1269,
      "of": 1272
    },
    "RHODES_HARTRAO_2326MHz.fits": {
      "pixels_kept_in_tris_band": 64,
      "of": 1272
    },
    "HASLAM_DESTRIPED_ONLY.fits": {
      "pixels_kept_in_tris_band": 1272,
      "of": 1272
    }
  }
}
```

## Step 2 — Ring TODs

Written 2026-11-20T14:03:24+00:00 (run_pipeline.py --step ring_tods).

```json
{
  "ring_samples": 120,
  "ra_deg": [
    0.0,
    357.0
  ],
  "data_sets": {
    "TRIS_600": {
      "freq_mhz": 600.5,
      "samples_with_a_point": 120,
      "median_T_K": 12.353404332430404,
      "median_sigma_K": 0.06675327707311456
    },
    "TRIS_820": {
      "freq_mhz": 817.8,
      "samples_with_a_point": 120,
      "median_T_K": 5.071597005680478,
      "median_sigma_K": 0.4301406746635338
    },
    "LWA_35MHz.fits": {
      "freq_mhz": 35.0,
      "samples_with_a_point": 120,
      "median_T_K": 14102.273066827263,
      "median_sigma_K": 2820.454613365453
    },
    "LWA1_38.fits": {
      "freq_mhz": 38.0,
      "samples_with_a_point": 120,
      "median_T_K": 12253.40745709792,
      "median_sigma_K": 2450.681491419584
    },
    "LWA1_40.fits": {
      "freq_mhz": 40.0,
      "samples_with_a_point": 120,
      "median_T_K": 10600.395926950001,
      "median_sigma_K": 2120.07918539
    },
    "LWA1_45.fits": {
      "freq_mhz": 45.0,
      "samples_with_a_point": 120,
      "median_T_K": 8059.664374674119,
      "median_sigma_K": 1611.932874934824
    },
    "MAIPU_45mhZ.fits": {
      "freq_mhz": 45.0,
      "samples_with_a_point": 120,
      "median_T_K": 7499.8886907003925,
      "median_sigma_K": 133.32954513470816
    },
    "LWA1_50.fits": {
      "freq_mhz": 50.0,
      "samples_with_a_point": 120,
      "median_T_K": 6251.03169123441,
      "median_sigma_K": 1250.2063382468823
    },
    "LWA1_60.fits": {
      "freq_mhz": 60.0,
      "samples_with_a_point": 120,
      "median_T_K": 4139.134743715687,
      "median_sigma_K": 827.8269487431375
    },
    "LWA1_70.fits": {
      "freq_mhz": 70.0,
      "samples_with_a_point": 120,
      "median_T_K": 2843.986257367901,
      "median_sigma_K": 568.7972514735802
    },
    "LWA1_74.fits": {
      "freq_mhz": 74.0,
      "samples_with_a_point": 120,
      "median_T_K": 2518.1048549991337,
      "median_sigma_K": 503.6209709998268
    },
    "LWA_80MHz.fits": {
      "freq_mhz": 80.0,
      "samples_with_a_point": 120,
      "median_T_K": 2135.893885375878,
      "median_sigma_K": 427.1787770751756
    },
    "CHIPASS_1400MHz.fits": {
      "freq_mhz": 1400.0,
      "samples_with_a_point": 0,
      "median_T_K": null,
      "median_sigma_K": null
    },
    "STOCKERT_VILLA_1400MHZ.fits": {
      "freq_mhz": 1420.0,
      "samples_with_a_point": 120,
      "median_T_K": 1.0680308887107668,
      "median_sigma_K": 0.5031309222941766
    },
    "RHODES_HARTRAO_2326MHz.fits": {
      "freq_mhz": 2326.0,
      "samples_with_a_point": 0,
      "median_T_K": null,
      "median_sigma_K": null
    }
  },
  "haslam": {
    "file": "HASLAM_DESTRIPED_ONLY.fits",
    "samples": 120,
    "median_T_K": 28.918145743180332,
    "median_sigma_K": 4.166844393221286
  },
  "sky_models_408": {
    "Haslam (uncorrected)": 28.918145743180546,
    "Haslam, Carter correction": 30.66677196973287,
    "Haslam, Irfan correction": 22.73384701037653,
    "GSM2008": 28.910596623252285,
    "GSM2016": 31.33275223029157,
    "GSM2016, GINAN correction": 32.79683607424917,
    "LFSM": 27.624243690194067
  },
  "irfan_split_tris": {
    "GSM2008 TRIS_600": {
      "median total (K)": -2.1888789956021455,
      "median offset part (K)": -2.2308536866878557,
      "median spectrum part (K)": 0.04289165329810185,
      "range (max \u2212 min) total (K)": 0.7523169975485611,
      "range (max \u2212 min) offset part (K)": 0.39161567489843696,
      "range (max \u2212 min) spectrum part (K)": 1.1171359490021384
    },
    "GSM2008 TRIS_820": {
      "median total (K)": -0.9198297221252976,
      "median offset part (K)": -0.9785856361912364,
      "median spectrum part (K)": 0.05998409191688836,
      "range (max \u2212 min) total (K)": 0.6684415493116325,
      "range (max \u2212 min) offset part (K)": 0.3637032502284101,
      "range (max \u2212 min) spectrum part (K)": 1.0079008536666567
    },
    "GSM2016 TRIS_600": {
      "median total (K)": -5.154496834865508,
      "median offset part (K)": -2.230853170676771,
      "median spectrum part (K)": -2.933623998249864,
      "range (max \u2212 min) total (K)": 1.9230693535847871,
      "range (max \u2212 min) offset part (K)": 0.39161455683465896,
      "range (max \u2212 min) spectrum part (K)": 2.2688030450328682
    },
    "GSM2016 TRIS_820": {
      "median total (K)": -4.113727186385985,
      "median offset part (K)": -0.97858520241182,
      "median spectrum part (K)": -3.1446695227827828,
      "range (max \u2212 min) total (K)": 1.498655485492769,
      "range (max \u2212 min) offset part (K)": 0.3637022652832034,
      "range (max \u2212 min) spectrum part (K)": 1.8184288659103722
    }
  }
}
```

## Step 3 — Model comparison (TRIS excluded)

Written 2026-11-16T19:55:12+00:00 (run_pipeline.py --step compare).

**Evidence of every model** (UltraNest ln Z summed over 120 ring samples; 1320 points in all; TRIS excluded):

| model | parameters | Σ ln Z | ± (error) | Δ ln Z vs best | samples where it has the highest ln Z |
|---|---|---|---|---|---|
| power law | 2 | -10424.45 | 2.153 | 0 | 98 of 120 |
| curved power law | 3 | -10460.72 | 2.294 | -36.27 | 22 of 120 |
| moment expansion (order 2) | 3 | -10565.16 | 2.367 | -140.72 | 0 of 120 |
| moment expansion (order 3) | 4 | -10573.89 | 2.368 | -149.44 | 0 of 120 |

```json
{
  "statistic": "evidence",
  "tris_excluded": true,
  "data_sets": [
    "LWA_35MHz.fits",
    "LWA1_38.fits",
    "LWA1_40.fits",
    "LWA1_45.fits",
    "MAIPU_45mhZ.fits",
    "LWA1_50.fits",
    "LWA1_60.fits",
    "LWA1_70.fits",
    "LWA1_74.fits",
    "LWA_80MHz.fits",
    "STOCKERT_VILLA_1400MHZ.fits"
  ],
  "units": 120,
  "skipped_units": 0,
  "points_per_unit": [
    11
  ],
  "sum_logZ": {
    "power law": -10424.446896774249,
    "curved power law": -10460.721760251216,
    "moment expansion (order 2)": -10565.162228608742,
    "moment expansion (order 3)": -10573.890676718751
  }
}
```

## Step 4 — Model selection

Written 2026-11-16T19:55:14+00:00 (run_pipeline.py --step select).

**Selection criteria** (rule: highest score; any model within tie_lnz = 1 of it with fewer parameters wins):

| criterion | setting / result |
|---|---|
| statistic compared | evidence (sum of ln Z over 120 ring samples) |
| 1. highest score | power law (score -10424.45) |
| 2. models within `select.tie_lnz` of the best | tie_lnz = 1.0: power law |
| 3. tie-break among those | fewest parameters, then highest score |
| chosen model | power law (parameters: A_K, beta) |

| model | parameters | score | score − best | within tie_lnz | samples won | chosen |
|---|---|---|---|---|---|---|
| power law | 2 | -10424.45 | 0 | yes | 98 | **yes** |
| curved power law | 3 | -10460.72 | -36.27 | no | 22 |  |
| moment expansion (order 2) | 3 | -10565.16 | -140.72 | no | 0 |  |
| moment expansion (order 3) | 4 | -10573.89 | -149.44 | no | 0 |  |

```json
{
  "model": "power law",
  "parameters": [
    "A_K",
    "beta"
  ],
  "statistic": "evidence (sum of ln Z over 120 ring samples)",
  "rule": "highest score; any model within tie_lnz = 1 of it with fewer parameters wins",
  "scores": {
    "power law": -10424.44689677425,
    "curved power law": -10460.721760251214,
    "moment expansion (order 2)": -10565.162228608744,
    "moment expansion (order 3)": -10573.89067671875
  },
  "score_minus_best": {
    "power law": 0.0,
    "curved power law": -36.27486347696322,
    "moment expansion (order 2)": -140.71533183449355,
    "moment expansion (order 3)": -149.4437799444986
  },
  "units_won": {
    "power law": 98,
    "curved power law": 22,
    "moment expansion (order 2)": 0,
    "moment expansion (order 3)": 0
  },
  "candidates_within_tie": [
    "power law"
  ]
}
```

## Step 5 — Our spectral model (TRIS included)

Written 2026-11-16T19:58:57+00:00 (run_pipeline.py --step fit).

```json
{
  "model": "power law",
  "parameters": [
    "A_K",
    "beta"
  ],
  "units_fitted": 120,
  "reused_from_cache": 0,
  "points_per_unit": [
    13
  ],
  "data_sets": [
    "LWA1_38.fits",
    "LWA1_40.fits",
    "LWA1_45.fits",
    "LWA1_50.fits",
    "LWA1_60.fits",
    "LWA1_70.fits",
    "LWA1_74.fits",
    "LWA_35MHz.fits",
    "LWA_80MHz.fits",
    "MAIPU_45mhZ.fits",
    "STOCKERT_VILLA_1400MHZ.fits",
    "TRIS_600",
    "TRIS_820"
  ],
  "median_parameters_maxlike": {
    "A_K": 32.14969680553959,
    "beta": -2.4865722154330734
  },
  "median_parameter_std": {
    "A_K": 0.16835815214236347,
    "beta": 0.007046582417976071
  },
  "median_chi2_per_dof": 0.5932041955893834
}
```

## Step 6 — Haslam at 408 MHz

Written 2026-11-16T19:58:58+00:00 (run_pipeline.py --step predict_408).

**Top 3 sky models against our prediction at 408 MHz** (model − prediction along the ring, TRIS beam; ranked by χ² with the prediction's σ; all 7 in `residuals_408.csv`):

| rank | sky model | χ²/N | median (K) | mean (K) | rms (K) | samples within 1σ | N |
|---|---|---|---|---|---|---|---|
| 1 | GSM2016 | 133.56 | -0.8351 | -1.28 | 2.212 | 5 | 120 |
| 2 | Haslam, Carter correction | 139.91 | -1.841 | -1.574 | 1.995 | 10 | 120 |
| 3 | LFSM | 387.20 | -3.001 | -3.221 | 3.359 | 0 | 120 |

```json
{
  "model": "power law",
  "median_prediction_K": 32.14969680553959,
  "median_prediction_sigma_K": 0.16835815214236347,
  "prediction_sigma_range_K": [
    0.15008889334210676,
    0.2202025725885234
  ],
  "haslam_minus_prediction": {
    "median_K": -3.250624809286233,
    "rms_K": 3.4940449082678873
  },
  "residuals_408": {
    "GSM2016": {
      "median_K": -0.8350627496653686,
      "rms_K": 2.2123694730928576,
      "chi2_per_N": 133.56172999552402,
      "rank": 1
    },
    "Haslam, Carter correction": {
      "median_K": -1.841278315081139,
      "rms_K": 1.995175139841652,
      "chi2_per_N": 139.90528243932118,
      "rank": 2
    },
    "LFSM": {
      "median_K": -3.001010143319788,
      "rms_K": 3.3585710472051593,
      "chi2_per_N": 387.2014838523264,
      "rank": 3
    },
    "Haslam (uncorrected)": {
      "median_K": -3.250624809286448,
      "rms_K": 3.4940449082678806,
      "chi2_per_N": 419.282348934845,
      "rank": 4
    },
    "GSM2008": {
      "median_K": -3.2625225288248965,
      "rms_K": 3.498956859769631,
      "chi2_per_N": 420.87331594123754,
      "rank": 5
    },
    "GSM2016, GINAN correction": {
      "median_K": 0.5043773608980491,
      "rms_K": 5.6729834721198396,
      "chi2_per_N": 1082.326402142309,
      "rank": 6
    },
    "Haslam, Irfan correction": {
      "median_K": -9.43492354209027,
      "rms_K": 9.558084030270308,
      "chi2_per_N": 3169.588096848587,
      "rank": 7
    }
  },
  "scale_offset_haslam": {
    "all": {
      "a": 1.016924281332967,
      "sigma_a": 0.038261520434066464,
      "b_K": -3.865317328780222,
      "sigma_b_K": 1.13766092251921,
      "N": 120
    },
    "excluding Cygnus A \u00b1 25\u00b0": {
      "a": 1.0552394103659433,
      "sigma_a": 0.05365482994963238,
      "b_K": -4.802001627299003,
      "sigma_b_K": 1.4671846482657804,
      "N": 103
    }
  },
  "cygnus_window_samples": 17
}
```

## Step 7 — RHINO prediction

Written 2026-11-20T14:04:32+00:00 (run_pipeline.py --step rhino).

**Top 3 sky models against our prediction, TRIS beam (ring at Dec +42.4°)** (model − this work; ranked by χ²/N with this work's σ over every sample and frequency; 6 models in `rhino_tods_tris_beam.csv`):

| rank | sky model | χ²/N | median 60 MHz (K) | median 65 MHz (K) | median 70 MHz (K) | median 75 MHz (K) | median 80 MHz (K) | rms, all (K) | N |
|---|---|---|---|---|---|---|---|---|---|
| 1 | GSM2016 | 43.5 | 194.66 | 121.44 | 72.6 | 42.2 | 24.97 | 302.13 | 600 |
| 2 | GSM2008, Carter correction | 43.92 | 77.18 | 79.86 | 73.25 | 65.92 | 58.67 | 289.73 | 600 |
| 3 | GSM2008 | 52.56 | -129.78 | -96.92 | -74.76 | -59.41 | -43.98 | 301.16 | 600 |

**Top 3 sky models against our prediction, RHINO beam (Jodrell Bank zenith)** (model − this work; ranked by χ²/N with this work's σ over every sample and frequency; 6 models in `rhino_tods_rhino_beam.csv`):

| rank | sky model | χ²/N | median 60 MHz (K) | median 65 MHz (K) | median 70 MHz (K) | median 75 MHz (K) | median 80 MHz (K) | rms, all (K) | N |
|---|---|---|---|---|---|---|---|---|---|
| 1 | GSM2016 | 154.96 | 114.28 | 72.42 | 40.65 | 21.95 | 9.306 | 149.55 | 600 |
| 2 | GSM2008, Carter correction | 170.43 | 23.41 | 55.56 | 68.31 | 76.69 | 78.95 | 179.24 | 600 |
| 3 | GSM2008 | 312.34 | -185.64 | -118.10 | -76.55 | -47.04 | -28.18 | 184.34 | 600 |

```json
{
  "model": "power law",
  "frequencies_mhz": [
    60.0,
    65.0,
    70.0,
    75.0,
    80.0
  ],
  "deconvolution_prior": "GSM2016",
  "this_work": {
    "TRIS": {
      "60": {
        "median_K": 3703.782941291103,
        "median_sigma_K": 57.13682366177778
      },
      "65": {
        "median_K": 3030.277892404568,
        "median_sigma_K": 45.162879814988415
      },
      "70": {
        "median_K": 2522.5910478928845,
        "median_sigma_K": 36.44964581970706
      },
      "75": {
        "median_K": 2129.1612480418216,
        "median_sigma_K": 29.82954877838586
      },
      "80": {
        "median_K": 1816.811753852573,
        "median_sigma_K": 24.693461628518392
      }
    },
    "RHINO": {
      "60": {
        "median_K": 3945.2631452175724,
        "median_sigma_K": 16.97872114924321
      },
      "65": {
        "median_K": 3204.138999693376,
        "median_sigma_K": 13.274530443496014
      },
      "70": {
        "median_K": 2658.9670363896917,
        "median_sigma_K": 10.694383516925804
      },
      "75": {
        "median_K": 2216.225489595982,
        "median_sigma_K": 9.020512746597442
      },
      "80": {
        "median_K": 1937.629235303768,
        "median_sigma_K": 6.622827277927949
      }
    }
  },
  "median_model_minus_this_work_fraction_rhino_beam": {
    "GSM2008": [
      -0.05492197499122678,
      -0.04149825094133662,
      -0.032431301685706404,
      -0.02385367745010694,
      -0.016425315174409854
    ],
    "GSM2016": [
      0.039573053691545526,
      0.031051228695665317,
      0.02102317540320677,
      0.013802578499438112,
      0.007137980853402745
    ],
    "LFSM": [
      0.1482505563540319,
      0.15449239291252742,
      0.15664800335436185,
      0.17013397198617924,
      0.18989208807669686
    ],
    "GSM2008, Carter correction": [
      0.007590903087830217,
      0.02048064669673322,
      0.03040392303175035,
      0.03860877126982835,
      0.045466492018685925
    ],
    "GSM2008, Irfan correction": [
      -0.2657976985706264,
      -0.25814727594667475,
      -0.24907277734520605,
      -0.24593979329969018,
      -0.21915850166038686
    ],
    "GSM2016, GINAN correction": [
      0.14435637754473307,
      0.16080662360740822,
      0.15821228246572402,
      0.14957015774072235,
      0.1431108141954305
    ],
    "GSM2008, Irfan spectrum only": [
      -0.060737202889393174,
      -0.05189662064771802,
      -0.04097425072363027,
      -0.0365148312490997,
      -0.015874875753490784
    ]
  },
  "max_closure_K": 230.1214422574203,
  "median_closure_K": {
    "60": 32.32172333324161,
    "65": 22.05220701387634,
    "70": 15.251351768009954,
    "75": 10.741791973685622,
    "80": 7.907151433325339
  },
  "irfan_offset_at_408_K": -6.1842987328039225,
  "irfan_split": {
    "TRIS beam": {
      "total": [
        -877.4295294714818,
        -721.7176628386425,
        -605.8045763759159,
        -511.1349307762763,
        -436.45719029192713
      ],
      "offset part": [
        -769.6577482840443,
        -631.7933187647881,
        -526.0397132147231,
        -443.3382085907115,
        -377.5820464115787
      ],
      "spectrum part": [
        -125.2850110777681,
        -95.96366352697737,
        -79.12632496183437,
        -66.46000278203678,
        -57.91230578297848
      ]
    },
    "RHINO beam": {
      "total": [
        -740.2838285195046,
        -630.7832307878449,
        -538.4404035736245,
        -462.5798088695325,
        -393.2467493451177
      ],
      "offset part": [
        -807.8896926868085,
        -657.1333083713935,
        -546.1800561284754,
        -459.1687020483229,
        -392.19421436025266
      ],
      "spectrum part": [
        91.0245438140214,
        43.96278245070357,
        17.305859677893068,
        3.6503432759826637,
        0.9507086575370067
      ]
    }
  }
}
```
