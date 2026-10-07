"""Read every map in ``res/``, whatever its layout, into one common description.

The surveys in ``res/`` arrive in four layouts and state their metadata in
different places (or not at all).  :func:`describe` reads one file and returns
the map as a RING float array with unobserved sky as NaN, in kelvin where the
unit is known, together with its frequency, frame, native beam FWHM and CMB
convention, and where each of those came from.  :func:`discover` does it for a
whole directory.  The same logic is shown, and checked against the sky, in
``notebooks/01_explore_data`` section 05.

Layouts handled:

* a standard HEALPix table -> ``hp.read_map``;
* an image in the FITS-standard HPX projection (CHIPASS) -> converted back to
  HEALPix pixel for pixel, with a one-to-one check;
* a HEALPix table without layout keywords (Rhodes/HartRAO) -> nside from the
  table size, RING vs NESTED chosen from the data.

Facts no file carries are in :data:`EXTERNAL_METADATA`, each with its source.
``cmb`` is ``"included"`` (subtract it before a spectral fit), ``"removed"``
(already gone) or ``"anisotropy"`` (a CMB-only map, no foreground).  ``usable``
is ``False`` for maps that cannot be a spectral-fit point, with the reason.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

FITS_NULL = -32768.0
SENTINEL_BELOW = -1e20            # hp.UNSEEN (-1.6375e30), CHIPASS's -1.0e30 / -1.1e30
ZERO_AS_MASK_FRACTION = 0.01      # exact zeros over >1% of the sky are a mask, not sky
UNIT_TO_KELVIN = {"K": 1.0, "KELVIN": 1.0, "KELVINS": 1.0, "MK": 1e-3}
FRAME_WORDS = {"G": "G", "GAL": "G", "GALACTIC": "G", "C": "C", "Q": "C", "EQ": "C",
               "EQUATORIAL": "C", "CELESTIAL": "C", "E": "E", "ECLIPTIC": "E"}

#: Facts the files do not carry themselves, keyed by filename, each with its source.
#: ``fwhm_deg`` is the native beam used for smoothing (the larger axis where the
#: beam is elliptical, so nothing is left under-smoothed).
EXTERNAL_METADATA: Dict[str, Dict[str, object]] = {
    "HASLAM.fits": dict(
        name="Haslam 408 MHz", cmb="included", role="target",
        source="CMB included: config haslam.remove_cmb_monopole"),
    "HASLAM_DESTRIPED_ONLY.fits": dict(
        name="Haslam 408 MHz, destriped only", cmb="included", role="target",
        source="Remazeilles et al. 2015, destriped but not desourced"),
    "MAIPU_45mhZ.fits": dict(
        name="Maipu + MU 45 MHz", unit="K", fwhm_deg=4.6, cmb="included",
        fwhm="4.6° × 2.4° (Maipu) / 3.6° (MU)",
        source="LAMBDA fg_maipu_info; Guzmán et al. 2011"),
    "LWA_35MHz.fits": dict(
        name="LWA1 35 MHz", frame="C", cmb="included",
        source="LAMBDA fg_lwa1_radio_maps_info: Equatorial; Dowell et al. 2017"),
    "LWA_80MHz.fits": dict(
        name="LWA1 80 MHz", frame="C", cmb="included",
        source="LAMBDA fg_lwa1_radio_maps_info: Equatorial; Dowell et al. 2017"),
    "STOCKERT_VILLA_1400MHZ.fits": dict(
        name="Stockert + Villa-Elisa 1420 MHz", unit="mK", freq_mhz=1420.0,
        fwhm_deg=35.4 / 60.0, fwhm="35.4′", cmb="included",
        source="LAMBDA fg_stockert_villa_info; Reich 1982, Testori et al. 2001"),
    "CHIPASS_1400MHz.fits": dict(
        name="CHIPASS 1400 MHz", unit="mK", fwhm_deg=14.4 / 60.0, fwhm="14.4′",
        cmb="included",
        source="LAMBDA fg_chipass_info; Calabretta et al. 2014. CMB measured as "
               "included (matches Stockert to +0.01 K)"),
    "RHODES_HARTRAO_2326MHz.fits": dict(
        name="Rhodes/HartRAO 2326 MHz", cmb="removed",
        source="Jonas et al. 1998, Platania et al. 2003 (header). CMB measured as "
               "removed (0.17 K at |b| > 40°)"),
    "ARCADE2_315.fits": dict(
        name="ARCADE 2 3.15 GHz", unit="K", thermodynamic=True, cmb="included",
        source="LAMBDA arcade_maps_info (K thermodynamic)"),
    "ARCADE2_341.fits": dict(
        name="ARCADE 2 3.41 GHz", unit="K", thermodynamic=True, cmb="included",
        source="LAMBDA arcade_maps_info (K thermodynamic)"),
    "WMAP_K.fits": dict(
        name="WMAP K band (header says V band)", frame="G", thermodynamic=True,
        cmb="included", usable=False,
        why="frequency contradicted (header FREQ 'V band', filename K) and the "
            "monopole is removed, so there is no absolute zero level",
        source="WMAP 7-yr band map, mK thermodynamic, monopole removed"),
    "WMAP_Ka.fits": dict(
        name="WMAP Ka band (header says V band)", frame="G", thermodynamic=True,
        cmb="included", usable=False,
        why="frequency contradicted (header FREQ 'V band', filename Ka) and the "
            "monopole is removed, so there is no absolute zero level",
        source="WMAP 7-yr band map, mK thermodynamic, monopole removed"),
    "WMAP_Q_BAND.fits": dict(
        name="WMAP Q band, foreground-cleaned", freq_mhz=40700.0, fwhm_deg=0.49,
        fwhm="0.49°", frame="G", thermodynamic=True, cmb="anisotropy", usable=False,
        why="foreground-cleaned CMB anisotropy: no foreground to fit",
        source="WMAP 9-yr: Q band 40.7 GHz, 0.49° beam, Galactic, mK thermodynamic"),
}


def _header_text(path: Path) -> Tuple[Dict[str, object], List[str]]:
    from astropy.io import fits

    keys: Dict[str, object] = {}
    lines: List[str] = []
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            for k in hdu.header:
                if k not in ("COMMENT", "HISTORY"):
                    keys[k] = hdu.header[k]
            lines += [str(c) for c in hdu.header.get("COMMENT", [])]
            lines += [str(c) for c in hdu.header.get("HISTORY", [])]
    return keys, lines


def _comment_value(lines: List[str], label: str) -> Optional[str]:
    for line in lines:
        hit = re.match(rf"\s*{label}\s*[-:=]\s*(.+)", line, flags=re.I)
        if hit:
            return hit.group(1).strip()
    return None


def _neighbour_corr(m: np.ndarray, nside: int, nest: bool) -> float:
    import healpy as hp

    nb = hp.get_all_neighbours(nside, np.arange(m.size), nest=nest)[0]
    ok = (nb >= 0) & np.isfinite(m) & np.isfinite(m[nb]) & (m != 0) & (m[nb] != 0)
    return float(np.corrcoef(m[ok], m[nb[ok]])[0, 1])


def read_sky(path: Path) -> Tuple[np.ndarray, int, str, Optional[str], Optional[str]]:
    """``(RING map as float, nside, how it was read, frame hint, hint source)``."""
    import healpy as hp
    from astropy.io import fits
    from astropy.wcs import WCS

    with fits.open(path, memmap=True) as hdul:
        h0 = hdul[0].header
        ctype = str(h0.get("CTYPE1", ""))
        if h0.get("NAXIS", 0) == 2 and ctype.endswith("-HPX"):
            img, wcs = np.asarray(hdul[0].data, dtype=float), WCS(h0)
            nside = img.shape[1] // 5                      # the HPX layout is 5 nside wide
            lon, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
            x, y = wcs.world_to_pixel_values(lon, lat)
            xi, yi = np.rint(x).astype(int), np.rint(y).astype(int)
            if np.unique(yi * img.shape[1] + xi).size != xi.size:
                raise ValueError("{}: HPX image is not a one-to-one HEALPix layout".format(path.name))
            hint = {"GLON": "G", "RA--": "C", "ELON": "E"}.get(ctype[:4])
            return img[yi, xi], nside, "HPX image -> HEALPix, 1:1", hint, "CTYPE1 " + ctype
        tables = [i for i, h in enumerate(hdul) if isinstance(h, fits.BinTableHDU)]
        if not tables:
            raise ValueError("{}: neither a HEALPix table nor an HPX image".format(path.name))
        t = tables[0]
        th = hdul[t].header
        if "NSIDE" in th and th.get("TTYPE1"):
            sky = hp.read_map(path, field=0, hdu=t, dtype=np.float64)
            return (sky, hp.get_nside(sky),
                    "HEALPix table, {}".format(str(th.get("ORDERING", "?")).strip()), None, None)
        if not str(th["TFORM1"]).strip().endswith("E"):
            raise ValueError("{}: unexpected TFORM1 {!r}".format(path.name, th["TFORM1"]))
        repeat = int(re.match(r"\s*(\d*)", str(th["TFORM1"])).group(1) or 1)
        npix = int(th["NAXIS2"]) * repeat
        offset = hdul.fileinfo(t)["datLoc"]
    raw = np.fromfile(path, dtype=">f4", count=npix, offset=offset).astype(float)
    nside = hp.npix2nside(npix)
    r, n = _neighbour_corr(raw, nside, False), _neighbour_corr(raw, nside, True)
    sky = raw if r >= n else hp.reorder(raw, n2r=True)
    how = "table w/o layout keys; {} from data ({:.2f} vs {:.2f})".format(
        "RING" if r >= n else "NESTED", max(r, n), min(r, n))
    return sky, nside, how, None, None


def _frequency_mhz(keys, lines, filename):
    for key, scale in (("FREQ", None), ("RESTFREQ", 1e-6), ("FREQUENC", None)):
        if key in keys:
            text = str(keys[key]).strip()
            hit = re.match(r"([\d.]+(?:[eE][+-]?\d+)?)", text)
            if hit:
                value = float(hit.group(1))
                if scale:
                    return value * scale, "header " + key
                return value * (1e3 if "GHZ" in text.upper() else 1.0), "header " + key
    text = _comment_value(lines, "Frequency")
    if text:
        hit = re.match(r"([\d.]+)\s*(MHZ|GHZ)", text.upper())
        if hit:
            return float(hit.group(1)) * (1e3 if hit.group(2) == "GHZ" else 1.0), "header comment"
    hit = re.search(r"(\d+(?:\.\d+)?)\s*(MHZ|GHZ)", filename.upper())
    if hit:
        return float(hit.group(1)) * (1e3 if hit.group(2) == "GHZ" else 1.0), "filename"
    return None, None


def _beam(keys, lines):
    """``(text, fwhm in degrees, source)`` -- the larger axis if elliptical."""
    if "BEAMSIZE" in keys:                                   # Haslam: arcmin
        v = float(keys["BEAMSIZE"]) / 60.0
        return "{:.0f}′".format(v * 60), v, "header BEAMSIZE"
    if "BEAMSZ" in keys:                                     # ARCADE 2: degrees
        v = float(keys["BEAMSZ"])
        return "{:.1f}°".format(v), v, "header BEAMSZ"
    if "BMAJ" in keys:                                       # LWA1: degrees
        a, b = float(keys["BMAJ"]), float(keys.get("BMIN", keys["BMAJ"]))
        text = "{:.2f}°".format(a) if a == b else "{:.2f}° × {:.2f}°".format(a, b)
        return text, max(a, b), "header BMAJ/BMIN"
    text = _comment_value(lines, "Resolution")              # Rhodes: "20 arcmin"
    if text:
        hit = re.match(r"([\d.]+)\s*(ARCMIN|DEG)", text.upper())
        if hit:
            v = float(hit.group(1)) / (60.0 if hit.group(2) == "ARCMIN" else 1.0)
            return text, v, "header comment"
    return None, None, None


def _unit(keys, lines):
    for key in ("TUNIT1", "BUNIT"):
        if str(keys.get(key, "")).strip():
            return str(keys[key]).strip(), "header " + key
    text = _comment_value(lines, "Pixel Units")
    return (text, "header comment") if text else (None, None)


def _frame(keys, lines, hint, hint_src):
    for key in ("COORDSYS", "SKYCOORD"):
        if key in keys:
            word = str(keys[key]).strip().upper()
            return FRAME_WORDS.get(word, FRAME_WORDS.get(word[:1])), "header " + key
    if hint:
        return hint, hint_src
    text = _comment_value(lines, "Coordinate system")
    if text and text.upper() in FRAME_WORDS:
        return FRAME_WORDS[text.upper()], "header comment"
    return None, None


def _unobserved(sky, keys):
    mask = ~np.isfinite(sky) | (sky < SENTINEL_BELOW) | (sky == FITS_NULL)
    if "TNULL1" in keys and float(keys["TNULL1"]) != 0.0:
        mask |= sky == float(keys["TNULL1"])
    zeros = sky == 0
    if zeros.mean() > ZERO_AS_MASK_FRACTION or keys.get("TNULL1") == 0:
        mask |= zeros
    return mask


def describe(path: Path) -> Dict[str, object]:
    """Read one map and everything known about it.  ``sky`` is in K where possible."""
    path = Path(path)
    keys, lines = _header_text(path)
    extra = EXTERNAL_METADATA.get(path.name, {})
    sky, nside, how, hint, hint_src = read_sky(path)

    freq, freq_src = _frequency_mhz(keys, lines, path.name)
    if freq is None and "freq_mhz" in extra:
        freq, freq_src = extra["freq_mhz"], "external"
    elif freq is not None and "freq_mhz" in extra and abs(float(extra["freq_mhz"]) - freq) > 1:
        freq, freq_src = extra["freq_mhz"], "external ({} says {:g})".format(freq_src, freq)
    unit, unit_src = _unit(keys, lines)
    if unit is None and "unit" in extra:
        unit, unit_src = extra["unit"], "external"
    fwhm, fwhm_deg, fwhm_src = _beam(keys, lines)
    if fwhm_deg is None and "fwhm_deg" in extra:
        fwhm, fwhm_deg, fwhm_src = extra.get("fwhm"), float(extra["fwhm_deg"]), "external"
    frame, frame_src = _frame(keys, lines, hint, hint_src)
    if frame is None and "frame" in extra:
        frame, frame_src = extra["frame"], "external"

    masked = _unobserved(sky, keys)
    scale = UNIT_TO_KELVIN.get(str(unit).upper()) if unit else None
    usable, why = bool(extra.get("usable", True)), extra.get("why")
    if usable and (freq is None or scale is None or frame is None or fwhm_deg is None):
        missing = [k for k, v in (("frequency", freq), ("unit", scale), ("frame", frame),
                                  ("beam", fwhm_deg)) if v is None]
        usable, why = False, "unknown " + ", ".join(missing)
    return dict(
        file=path.name, path=path, name=extra.get("name", path.stem), freq_mhz=freq,
        freq_src=freq_src, nside=nside, how=how, frame=frame, frame_src=frame_src,
        unit=unit, unit_src=unit_src, in_kelvin=scale is not None, fwhm=fwhm,
        fwhm_deg=fwhm_deg, fwhm_src=fwhm_src, cmb=extra.get("cmb", "included"),
        thermodynamic=bool(extra.get("thermodynamic", False)),
        role=extra.get("role", "fit"), usable=usable, why=why,
        source=extra.get("source"), n_unobserved=int(masked.sum()),
        sky=np.where(masked, np.nan, sky * (scale if scale else 1.0)),
    )


def discover(directory: Path) -> List[Dict[str, object]]:
    """:func:`describe` every ``*.fits`` in ``directory``, sorted by frequency."""
    maps = [describe(p) for p in sorted(Path(directory).glob("*.fits"))]
    return sorted(maps, key=lambda m: (m["freq_mhz"] is None, m["freq_mhz"] or 0.0))


# ---------------------------------------------------------------------------
# temperature conventions for a spectral fit
# ---------------------------------------------------------------------------
_H, _K, T_CMB = 6.62607015e-34, 1.380649e-23, 2.72548


def antenna_temperature(t_thermo_k, nu_mhz):
    """Thermodynamic -> Rayleigh-Jeans: T_A = T x / (e^x - 1), x = h nu / k T."""
    t = np.asarray(t_thermo_k, dtype=float)
    x = _H * nu_mhz * 1e6 / (_K * t)
    return t * x / np.expm1(x)


def cmb_rj_k(nu_mhz: float) -> float:
    """The CMB monopole as a Rayleigh-Jeans temperature at ``nu_mhz``."""
    return float(antenna_temperature(T_CMB, nu_mhz))


def rj_excess(t_k, m: Dict[str, object]):
    """RJ temperature with the CMB removed, per the map's own ``cmb`` flag."""
    nu = float(m["freq_mhz"])
    if m["cmb"] == "removed":
        return np.asarray(t_k, dtype=float)
    if m["thermodynamic"]:
        return antenna_temperature(t_k, nu) - antenna_temperature(T_CMB, nu)
    return np.asarray(t_k, dtype=float) - cmb_rj_k(nu)
