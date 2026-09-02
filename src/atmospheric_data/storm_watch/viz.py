"""Storm-watch map: draw an NWS alert and the automatically-built simulation domain.

A light, dependency-free (matplotlib-only) situational-awareness figure -- the honest inputs the
automation acts on, not a rendered simulation.  It shows the alert polygon, the auto-domain box
(centred downstream along the storm motion), the storm-motion vector, and the nearest NEXRAD
radars + METAR stations the domain builder selected.  No cartopy: a plain lon/lat plot with an
aspect correction, so it runs anywhere the rest of the package does."""
from __future__ import annotations

import math
import os

from . import alerts as alerts_mod
from . import domain as dom
from .config import StormWatchConfig

_LEVEL_COLOR = {"confirmed": "#b8143c", "warning": "#e8663c", "watch": "#e3b23c", "info": "#5b8c5a"}


def plot_alert_domain(alert, cfg, out_path, level=None):
    """Render the alert + auto-domain map for ``alert`` to ``out_path`` (PNG). Returns the path."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Polygon as MplPoly

    if not alert.polygon:
        raise ValueError("alert has no polygon; nothing to map")
    d = dom.build_domain(alert, cfg)
    clat, clon = d["center_lat"], d["center_lon"]
    radars = dom.nearest_radars(clat, clon, cfg)
    metars = dom.nearest_metars(clat, clon)
    ex, ny = d["storm_motion_unit"]
    level = level or alerts_mod.classify_level(alert, cfg)
    col = _LEVEL_COLOR.get(level, "#e8663c")
    coslat = math.cos(math.radians(clat)) or 1.0

    poly = [(lon, lat) for lon, lat in alert.polygon]
    plons = [p[0] for p in poly]; plats = [p[1] for p in poly]
    hw = 0.5 * d["width_km"] / (111.0 * coslat); hh = 0.5 * d["width_km"] / 111.0   # square domain

    fig, ax = plt.subplots(figsize=(9, 8.4), facecolor="white")
    # auto-domain box (the simulation footprint)
    ax.add_patch(Rectangle((clon - hw, clat - hh), 2 * hw, 2 * hh, fill=False, ec="#2f6690",
                           lw=2.2, ls="--", zorder=4, label="auto-domain (%.0f km)" % d["width_km"]))
    # alert polygon (on top, bright -- it is the whole point; a county polygon is small vs the domain)
    ax.add_patch(MplPoly(poly, closed=True, facecolor=col, edgecolor=col, alpha=0.30, lw=2.6, zorder=7))
    ax.plot(plons + [plons[0]], plats + [plats[0]], color=col, lw=2.6, zorder=7, label="alert polygon")
    # storm-motion vector from the polygon centroid
    cc = dom.polygon_centroid(alert.polygon)
    if cc:
        L = 0.45 * d["width_km"] / 111.0
        ax.annotate("", xy=(cc[1] + L * ex / coslat, cc[0] + L * ny), xytext=(cc[1], cc[0]),
                    arrowprops=dict(arrowstyle="-|>", color="#333", lw=2), zorder=6)
        ax.text(cc[1] + 0.5 * L * ex / coslat, cc[0] + 0.5 * L * ny + 0.02, "storm motion",
                fontsize=8, color="#333", zorder=6)
    _bb = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75)
    # nearest NEXRAD radars (label offset alternates to avoid overlap at the storm cluster)
    for k, (sid, km) in enumerate(radars):
        if sid in dom.NEXRAD_STATIONS:
            la, lo = dom.NEXRAD_STATIONS[sid]
            ax.plot(lo, la, "v", color="#1b4965", ms=12, zorder=5,
                    label="NEXRAD" if k == 0 else None)
            dy = 0.06 if k % 2 == 0 else -0.06; va = "bottom" if dy > 0 else "top"
            ax.text(lo, la + dy, "%s (%.0f km)" % (sid, km), fontsize=8, ha="center",
                    va=va, color="#1b4965", zorder=8, bbox=_bb)
    # nearest METAR stations
    for k, (sid, km) in enumerate(metars):
        if sid in dom.METAR_STATIONS:
            la, lo, _ = dom.METAR_STATIONS[sid]
            ax.plot(lo, la, "o", color="#5b8c5a", ms=5, zorder=5,
                    label="METAR" if k == 0 else None)
    ax.plot(clon, clat, "P", color="#2f6690", ms=12, zorder=6, label="domain centre")

    # frame: the domain box plus a small margin
    ax.set_xlim(clon - 1.18 * hw, clon + 1.18 * hw); ax.set_ylim(clat - 1.18 * hh, clat + 1.18 * hh)
    ax.set_aspect(1.0 / coslat)                              # degrees lon are shorter at this lat
    ax.set_xlabel("longitude [deg]"); ax.set_ylabel("latitude [deg]")
    ax.grid(alpha=0.3, zorder=0)
    title = "%s  [%s]" % (alert.event or "alert", level.upper())
    sub = (alert.headline or alert.affected_area or "")[:96]
    ax.set_title(title + ("\n" + sub if sub else ""), fontsize=11, fontweight="bold", loc="left")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.text(0.5, -0.10, "storm-watch situational map — real NWS alert + the auto-built simulation domain "
            "(not a simulation)", transform=ax.transAxes, ha="center", fontsize=8, color="#666")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=140, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def map_alert_file(path, cfg=None, out_path=None):
    """Load the first alert with a polygon from ``path`` and map it. Returns the PNG path."""
    cfg = cfg or StormWatchConfig()
    als = [a for a in alerts_mod.load_alerts_file(path) if a.polygon]
    if not als:
        raise ValueError("no alert with a polygon in %s" % path)
    alert = als[0]
    import re
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", alert.alert_id or "alert")[-60:]
    out_path = out_path or os.path.join(getattr(cfg, "workdir", "outputs/storm_watch"),
                                        "map_%s.png" % safe)
    return plot_alert_domain(alert, cfg, out_path)
