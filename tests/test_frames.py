"""Coordinate-frame tests.

The algebra tests use synthetic maps and run anywhere.  The FITS tests pin what
the downloaded files declare, and are skipped when the files are absent.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

healpy = pytest.importorskip("healpy")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from moonrunes import frames  # noqa: E402
from moonrunes.stage2_beam_matching import observed_mask  # noqa: E402

RES = REPO_ROOT / "res"


def _galactic_patch(nside=16, lat_max=25.0, lon_max=60.0):
    """A smooth ~3 K map observed only inside a Galactic lon/lat box, zero elsewhere."""
    lon, lat = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)), lonlat=True)
    lon = (lon + 180.0) % 360.0 - 180.0
    observed = (np.abs(lat) < lat_max) & (np.abs(lon) < lon_max)
    sky = 2.8 + 0.3 * np.exp(-(lat / 5.0) ** 2) * np.cos(np.radians(lon))
    return np.where(observed, sky, 0.0), observed


def test_masked_rotation_never_leaks_the_zero_fill():
    sky, observed = _galactic_patch()
    rotated, observed_rotated, _ = frames.rotate_masked(sky, observed)
    inside = rotated[observed_rotated]
    assert np.all(inside != 0.0)
    assert sky[observed].min() <= inside.min()
    assert inside.max() <= sky[observed].max()
    # Everything outside is the unambiguous sentinel, which stage 2 reads as
    # unobserved under the ARCADE convention.
    assert np.all(rotated[~observed_rotated] == healpy.UNSEEN)
    assert np.array_equal(observed_mask(rotated, "zeros"), observed_rotated)
    # A rotation preserves area; the weight floor may move the edge by a few pixels.
    assert abs(int(observed_rotated.sum()) - int(observed.sum())) <= 0.05 * observed.sum()


def test_plain_pixel_rotation_does_leak():
    """The failure rotate_masked exists to prevent, so the test above means something."""
    sky, observed = _galactic_patch()
    naive = healpy.Rotator(coord=["G", "C"]).rotate_map_pixel(sky)
    assert np.any((naive > 0) & (naive < sky[observed].min()))


def test_masked_rotation_moves_the_galactic_centre_to_sgr_a():
    nside = 64
    sky = np.zeros(healpy.nside2npix(nside))
    observed = np.ones(sky.size, dtype=bool)
    sky[healpy.ang2pix(nside, 0.0, 0.0, lonlat=True)] = 1.0
    rotated, _, _ = frames.rotate_masked(sky, observed)
    ra, dec = healpy.pix2ang(nside, int(np.argmax(rotated)), lonlat=True)
    offset = healpy.rotator.angdist([ra, dec], [266.41683, -29.00781], lonlat=True)
    assert np.degrees(offset[0]) < healpy.nside2resol(nside, arcmin=True) / 60.0


def test_full_sky_rotation_refuses_masked_maps():
    sky, _ = _galactic_patch()
    sky[sky == 0.0] = healpy.UNSEEN
    with pytest.raises(ValueError, match="rotate_masked"):
        frames.rotate_full_sky(sky)


@pytest.mark.parametrize("filename, keyword, ordering, nside", [
    ("HASLAM.fits", "COORDSYS", "RING", 512),
    ("ARCADE2_315.fits", "SKYCOORD", "NESTED", 16),
    ("ARCADE2_341.fits", "SKYCOORD", "NESTED", 16),
])
def test_fits_files_declare_galactic(filename, keyword, ordering, nside):
    path = RES / filename
    if not path.exists():
        pytest.skip("{} not downloaded".format(filename))
    declared = frames.declared_frame(path)
    assert declared["frame"] == "G"
    assert {s["keyword"] for s in declared["statements"]} == {keyword}
    assert declared["ordering"] == ordering
    assert declared["nside"] == nside
