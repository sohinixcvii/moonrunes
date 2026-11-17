"""One palette, one colormap and one style file for every TRIS figure.

    from moonrunes.style import use_tris_style, COLORS
    use_tris_style()            # paper (one MNRAS column)
    use_tris_style(talk=True)   # slides: tris.mplstyle with tris_talk.mplstyle on top

Colormaps (registered with matplotlib by use_tris_style, or register_colormaps):
    "ring"      sequential; the default image.cmap of the style
    "ring_div"  diverging, for signed data: use it with CenteredNorm / TwoSlopeNorm(0)
                or symmetric vmin = -vmax
    and their reversed versions "ring_r", "ring_div_r".

COLORS names the colour of each role; use it wherever a colour must be set by hand.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

__all__ = ["COLORS", "TRIS10", "STYLE_FILE", "TALK_STYLE_FILE", "register_colormaps", "use_tris_style"]

_HERE = Path(__file__).resolve().parent
STYLE_FILE = _HERE / "tris.mplstyle"
TALK_STYLE_FILE = _HERE / "tris_talk.mplstyle"

# Same order as axes.prop_cycle in tris.mplstyle.
COLORS = {
    "data": "#E3A700",        # the TRIS data
    "survey": "#2E5C9A",
    "prediction": "#C8452D",
    "model": "#2A9D8F",
    "ink": "#1F2A44",         # text, axes, best-fit lines
    "group": "#B49AD6",
    "sky": "#7DB3E8",
    "rose": "#D0679A",
    "moss": "#7E8B2F",
    "umber": "#8B5A3C",
}

# 1-6: colour-vision safe; 7-10: need linestyle/marker as well
TRIS10 = list(COLORS.values())

_CMAPS = {
    "ring": ["#0B132B", "#5B3F8C", "#C8452D", "#FFD166"],
    "ring_div": ["#3E6FB0", "#F5F5F5", "#C8452D"],
}


def register_colormaps():
    """Register "ring" and "ring_div" (and their "_r" versions) with matplotlib (once)."""
    for name, colours in _CMAPS.items():
        cmap = LinearSegmentedColormap.from_list(name, colours)
        for c in (cmap, cmap.reversed()):
            if c.name not in mpl.colormaps:
                mpl.colormaps.register(c, name=c.name)


def use_tris_style(talk=False):
    """Register the colormaps, then apply the TRIS style (with the slide overlay if talk)."""
    register_colormaps()
    files = [STYLE_FILE, TALK_STYLE_FILE] if talk else [STYLE_FILE]
    plt.style.use([str(f) for f in files])
