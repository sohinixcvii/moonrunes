"""Coordinate frames -- put every map in equatorial (RA/Dec) before comparing.

What each input declares
------------------------
Read from the files, not from convention (``notebooks/04_coordinate_frames``):

=================  ======================  ========  =====  ==================
file               frame keyword           ordering  nside  mask on disk
=================  ======================  ========  =====  ==================
Haslam 408 MHz     ``COORDSYS='GALACTIC'`` RING      512    none (full sky)
ARCADE 2 3.15 GHz  ``SKYCOORD='Galactic'`` NESTED    16     exact 0.0
ARCADE 2 3.41 GHz  ``SKYCOORD='Galactic'`` NESTED    16     exact 0.0
TRIS stage 1       none -- ``.npz``        RING      8      NaN + band_mask
=================  ======================  ========  =====  ==================

ARCADE 2 *does* state its frame, under the non-standard ``SKYCOORD`` keyword
that healpy does not read -- which is why it looked unstated.  The TRIS
products carry no keyword at all; they are equatorial by construction, since
limTOD points the boresight at ``RA = LST, dec = latitude`` and builds the
operator on an equatorial pixel grid.

Which rotation, per map
-----------------------
healpy's own guidance is keyed on band-limitation and masking, not on nside:
``rotate_map_alms`` "is generally the best strategy", and pixel space is
better "for heavily masked maps where the spherical harmonics transform is not
well defined".  So:

* **Haslam** -- full sky, and band-limited at nside 512 (a 56 arcmin beam on
  6.9 arcmin pixels): :func:`rotate_full_sky`, harmonic space.  Its G->C->G
  round trip is ~6x more accurate than the pixel route, which smears compact
  sources through bilinear interpolation.
* **ARCADE 2** -- 93% of the sky is exact zero, a hard edge the harmonic
  transform would ring across: :func:`rotate_masked`, pixel space.

``hp.get_interp_val`` interpolates bilinearly over four neighbours, and it
knows nothing about masks: at the boundary it averages real data with the zero
fill and drags the edge pixels toward 0 K.  :func:`rotate_masked` interpolates
the data and the mask separately and divides, so every output value is a
weighted average of *observed* neighbours only.

Both functions take and return RING ordering.  ``hp.read_map`` already
reorders NESTED files to RING (it reads ``ORDERING``), and every consumer here
reads through it, so nothing converts back.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import numpy as np

#: The keywords a HEALPix FITS file uses to state its frame.  ``COORDSYS`` is
#: the HEALPix standard; ``SKYCOORD`` is what the ARCADE 2 archive writes.
FRAME_KEYWORDS = ("COORDSYS", "SKYCOORD")

_FRAME_CODES = {
    "G": "G", "GALACTIC": "G",
    "C": "C", "Q": "C", "EQUATORIAL": "C", "CELESTIAL": "C",
    "E": "E", "ECLIPTIC": "E",
}

#: Below this interpolated mask weight a rotated pixel is mask edge, not data.
#: Same meaning and value as ``beam_matching.weight_floor`` in stage 2.
DEFAULT_WEIGHT_FLOOR = 0.5


def declared_frame(path: Path) -> Dict[str, object]:
    """What a HEALPix FITS file states about itself, from every HDU.

    Returns the frame code healpy uses (``'G'``, ``'C'``, ``'E'``) or ``None``
    when no HDU states one, plus the keyword and HDU it came from, and the
    ordering and nside.  Raises if two HDUs disagree -- a file that contradicts
    itself is not something to pick a side of silently.
    """
    from astropy.io import fits

    found = []
    ordering = nside = None
    with fits.open(path) as hdus:
        for index, hdu in enumerate(hdus):
            header = hdu.header
            for keyword in FRAME_KEYWORDS:
                if keyword in header:
                    raw = str(header[keyword]).strip()
                    code = _FRAME_CODES.get(raw.upper())
                    if code is None:
                        raise ValueError(
                            "{} HDU {}: unrecognised {} = {!r}".format(
                                path, index, keyword, raw))
                    found.append((index, keyword, raw, code))
            ordering = header.get("ORDERING", ordering)
            nside = header.get("NSIDE", nside)

    codes = {code for *_, code in found}
    if len(codes) > 1:
        raise ValueError("{} states conflicting frames: {}".format(path, found))
    return {
        "frame": codes.pop() if codes else None,
        "statements": [
            {"hdu": i, "keyword": k, "value": v} for i, k, v, _ in found
        ],
        "ordering": None if ordering is None else str(ordering).strip(),
        "nside": None if nside is None else int(nside),
    }


def rotate_full_sky(
    sky: np.ndarray, frames: Tuple[str, str] = ("G", "C")
) -> np.ndarray:
    """Rotate a full-sky, band-limited RING map in harmonic space."""
    import healpy as hp

    sky = np.asarray(sky, dtype=float)
    if not np.all(np.isfinite(sky)) or np.any(sky == hp.UNSEEN):
        raise ValueError(
            "rotate_full_sky needs a map with no unobserved pixels; a masked "
            "map rings across its edge in harmonic space -- use rotate_masked"
        )
    return hp.Rotator(coord=list(frames)).rotate_map_alms(sky)


def rotate_masked(
    sky: np.ndarray,
    observed: np.ndarray,
    frames: Tuple[str, str] = ("G", "C"),
    weight_floor: float = DEFAULT_WEIGHT_FLOOR,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rotate a partial-sky RING map in pixel space without leaking the mask.

    Returns ``(rotated, observed_rotated, weight)``.  ``rotated`` is
    ``hp.UNSEEN`` wherever ``observed_rotated`` is False, i.e. wherever less
    than ``weight_floor`` of the bilinear interpolation weight fell on observed
    pixels.  Elsewhere it is a convex combination of observed input values, so
    it can never fall outside their range.
    """
    import healpy as hp

    sky = np.asarray(sky, dtype=float)
    observed = np.asarray(observed, dtype=bool)
    if sky.shape != observed.shape:
        raise ValueError("sky and observed must be the same shape")
    if not 0.0 < weight_floor <= 1.0:
        raise ValueError("weight_floor must be in (0, 1]")

    rotator = hp.Rotator(coord=list(frames))
    filled = np.where(observed, sky, 0.0)
    numerator = rotator.rotate_map_pixel(filled)
    weight = rotator.rotate_map_pixel(observed.astype(float))

    observed_rotated = weight >= weight_floor
    rotated = np.full(sky.shape, hp.UNSEEN)
    rotated[observed_rotated] = (
        numerator[observed_rotated] / weight[observed_rotated]
    )
    return rotated, observed_rotated, weight


def source_offsets_deg(
    sky: np.ndarray,
    sources: Dict[str, Tuple[float, float]],
    search_radius_deg: float,
) -> Dict[str, Dict[str, float]]:
    """Brightest pixel near each named source, and how far it is from it.

    ``sources`` maps a name to its literature ``(lon, lat)`` in degrees, in
    the map's own frame.  The peak is searched within ``search_radius_deg``;
    a large offset means the feature is not where the literature puts it.
    """
    import healpy as hp

    sky = np.asarray(sky, dtype=float)
    nside = hp.get_nside(sky)
    valid = np.isfinite(sky) & (sky != hp.UNSEEN)
    out: Dict[str, Dict[str, float]] = {}
    for name, (lon, lat) in sources.items():
        centre = hp.ang2vec(lon, lat, lonlat=True)
        disc = hp.query_disc(nside, centre, np.radians(search_radius_deg))
        disc = disc[valid[disc]]
        if disc.size == 0:
            out[name] = {"found": False}
            continue
        peak = disc[np.argmax(sky[disc])]
        peak_lon, peak_lat = hp.pix2ang(nside, peak, lonlat=True)
        offset = np.degrees(hp.rotator.angdist(
            hp.pix2vec(nside, peak), centre))
        out[name] = {
            "found": True,
            "peak_lon_deg": float(peak_lon),
            "peak_lat_deg": float(peak_lat),
            "peak_k": float(sky[peak]),
            "offset_deg": float(np.atleast_1d(offset)[0]),
            "pixel_deg": float(hp.nside2resol(nside, arcmin=True) / 60.0),
        }
    return out
