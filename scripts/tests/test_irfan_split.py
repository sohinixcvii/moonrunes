"""The Irfan-correction split (pipeline.irfan_split, used by steps 2 and 7).

total = corrected − base, offset part = corrected − spectrum only, spectrum part = spectrum only − base.  Every step
between the maps and the TODs is linear, so the parts observed through a beam add up to the total observed through it,
and the offset part is the 408 MHz offset times the beam-averaged spectral scaling.  Synthetic maps and a random
normalised beam operator; runs anywhere.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from moonrunes.pipeline import irfan_split  # noqa: E402


def _setup(seed=0, n=40, npix=300):
    rng = np.random.default_rng(seed)
    B = np.abs(rng.normal(size=(n, npix))) ** 3
    B /= B.sum(axis=1, keepdims=True)                     # a normalised beam operator, one row per ring sample
    base408 = rng.uniform(10.0, 100.0, npix)              # an uncorrected 408 MHz map [K]
    scale = rng.uniform(5.0, 40.0, npix)                  # a per-pixel S(nu) / S(408), as irfan_scaling returns
    base_nu = base408 * rng.uniform(4.0, 45.0, npix)      # the base model's own map at nu
    offset408 = -6.18                                     # -(8.9 K - T_CMB,RJ(408)) on a CMB-subtracted map
    corrected = (base408 + offset408) * scale             # the pipeline's Irfan-corrected model
    spectrum_only = base408 * scale                       # the same without the offset
    return B, scale, offset408, corrected, spectrum_only, base_nu


def test_offset_plus_spectrum_is_total_through_a_beam():
    B, scale, offset408, corrected, spectrum_only, base_nu = _setup()
    # each model observed separately, as the pipeline does
    total, offset, spectrum = irfan_split(B @ corrected, B @ spectrum_only, B @ base_nu)
    np.testing.assert_allclose(offset + spectrum, total, rtol=0, atol=1e-9 * np.abs(total).max())
    np.testing.assert_allclose(total, B @ (corrected - base_nu), rtol=1e-12)


def test_offset_part_is_the_offset_times_the_beam_averaged_scaling():
    B, scale, offset408, corrected, spectrum_only, base_nu = _setup(seed=1)
    _, offset, _ = irfan_split(B @ corrected, B @ spectrum_only, B @ base_nu)
    np.testing.assert_allclose(offset, offset408 * (B @ scale), rtol=1e-12)
    # uniform only where the scaling is: a constant scaling gives a constant offset part
    _, flat, _ = irfan_split(B @ ((spectrum_only / scale + offset408) * 7.0), B @ (spectrum_only / scale * 7.0),
                             B @ base_nu)
    np.testing.assert_allclose(flat, offset408 * 7.0, rtol=1e-12)


def test_nan_samples_stay_nan_and_a_broken_sum_raises():
    B, scale, offset408, corrected, spectrum_only, base_nu = _setup(seed=2)
    c, s, b = B @ corrected, B @ spectrum_only, B @ base_nu
    c[3] = s[3] = b[3] = np.nan                            # a sample below the coverage floor in all three
    total, offset, spectrum = irfan_split(c, s, b)
    assert np.isnan(total[3]) and np.isnan(offset[3]) and np.isnan(spectrum[3])
    with pytest.raises(AssertionError):                    # NaN in one model only: the parts no longer add up
        s2 = s.copy(); s2[5] = np.nan
        irfan_split(c, s2, b)


def test_step2_table_writes_one_row_per_model_and_band(tmp_path, monkeypatch):
    """Step 2's irfan_split_rows on synthetic TODs: every row has the same columns (_write_csv takes the
    columns from the first row) and the CSV writes."""
    import csv

    from moonrunes import pipeline

    healpy = pytest.importorskip("healpy")
    nside, n = 4, 30
    npix = healpy.nside2npix(nside)
    rng = np.random.default_rng(3)
    B = np.abs(rng.normal(size=(n, npix))) ** 3
    B /= B.sum(axis=1, keepdims=True)
    scales = {600.5: rng.uniform(0.3, 0.4, npix), 817.8: rng.uniform(0.14, 0.18, npix)}
    monkeypatch.setattr(pipeline, "irfan_scaling", lambda sky_cfg, ns, nu, from_mhz=408.0: scales[float(nu)])
    sky_cfg = {"irfan_monopole_k": 8.9}
    from moonrunes import sky_maps
    off408 = -(8.9 - float(sky_maps.cmb_rj_k(408.0)))
    nu_tris = {"TRIS_600": 600.5, "TRIS_820": 817.8}
    maps408, out = {}, {}
    for base in pipeline.IRFAN_SPLIT_BASES:
        m408 = rng.uniform(10.0, 60.0, npix)
        maps408[base] = (m408, "G", 1.0)
        for key_t, nu in nu_tris.items():
            out["modeltris::{}::{}".format(base, key_t)] = B @ (m408 * rng.uniform(0.2, 0.5, npix))
            out["modeltris::{}, Irfan correction::{}".format(base, key_t)] = B @ ((m408 + off408) * scales[nu])
            out["modeltris::{}::{}".format(pipeline.IRFAN_SPECTRUM_ONLY.format(base), key_t)] = B @ (m408 * scales[nu])
    observe = lambda sky, frame, native: (B @ np.asarray(sky, float), np.ones(n))
    rows = pipeline.irfan_split_rows(out, maps408, nu_tris, sky_cfg, observe, 0.5)
    assert len(rows) == len(pipeline.IRFAN_SPLIT_BASES) * len(nu_tris)
    assert all(list(r) == list(rows[0]) for r in rows)
    assert all(r["max |offset part − offset × scaling| (K)"] < 1e-9 for r in rows)
    pipeline._write_csv(tmp_path / "irfan_split_tris.csv", rows)
    with open(tmp_path / "irfan_split_tris.csv", newline="") as fh:
        back = list(csv.DictReader(fh))
    assert [r["base model"] for r in back] == [r["base model"] for r in rows]
