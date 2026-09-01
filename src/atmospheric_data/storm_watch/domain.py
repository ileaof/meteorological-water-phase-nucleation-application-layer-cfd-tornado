"""Automatic domain + asset selection from an alert (storm-watch).

From the alert polygon: the geographic centre; a download domain enlarged by configurable
**upstream / downstream / lateral** margins (the upstream margin follows the storm-motion
vector); the nearest NEXRAD radars (distance + count limited); the nearest METAR/ASOS stations;
and the most recent HRRR cycle covering the period.
"""
from __future__ import annotations

import re
import math

# A compact NEXRAD WSR-88D station table (central-US emphasis; extend as needed). lat, lon.
NEXRAD_STATIONS = {
    "KTLX": (35.333, -97.278), "KFDR": (34.362, -98.976), "KINX": (36.175, -95.564),
    "KVNX": (36.741, -98.128), "KFWS": (32.573, -97.303), "KDYX": (32.538, -99.254),
    "KLBB": (33.654, -101.814), "KAMA": (35.233, -101.709), "KDDC": (37.761, -99.969),
    "KICT": (37.655, -97.443), "KGLD": (39.367, -101.700), "KTWX": (38.997, -96.232),
    "KEAX": (38.810, -94.264), "KSGF": (37.235, -93.401), "KSRX": (35.290, -94.362),
    "KOUN": (35.236, -97.462), "KLZK": (34.836, -92.262), "KSHV": (32.451, -93.841),
}
# a few ASOS/METAR stations (id: lat, lon, elev_m)
METAR_STATIONS = {
    "KOKC": (35.393, -97.601, 392), "KOUN": (35.181, -97.439, 362), "KTUL": (36.199, -95.888, 205),
    "KICT": (37.650, -97.433, 408), "KDFW": (32.897, -97.038, 171), "KAMA": (35.219, -101.706, 1099),
}


def _haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1); dl = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(min(1.0, math.sqrt(a)))


def polygon_centroid(polygon):
    """(lat, lon) centroid of a [[lon,lat],...] ring (simple vertex mean; robust for warnings)."""
    if not polygon:
        return None
    lons = [pt[0] for pt in polygon]; lats = [pt[1] for pt in polygon]
    return sum(lats) / len(lats), sum(lons) / len(lons)


_DIR8 = {"n": (0, 1), "north": (0, 1), "ne": (0.71, 0.71), "northeast": (0.71, 0.71),
         "e": (1, 0), "east": (1, 0), "se": (0.71, -0.71), "southeast": (0.71, -0.71),
         "s": (0, -1), "south": (0, -1), "sw": (-0.71, -0.71), "southwest": (-0.71, -0.71),
         "w": (-1, 0), "west": (-1, 0), "nw": (-0.71, 0.71), "northwest": (-0.71, 0.71)}


def storm_motion(alert, default=(0.6, 0.6)):
    """Estimate the storm-motion **unit** vector (east, north) from the alert text
    (e.g. "moving northeast at 45 mph"); default is a NE mover.  Returns (ex, ny) unit."""
    m = re.search(r"moving\s+([a-z]+)\s+at", alert.text())
    if m and m.group(1) in _DIR8:
        vx, vy = _DIR8[m.group(1)]
        n = math.hypot(vx, vy) or 1.0
        return vx / n, vy / n
    n = math.hypot(*default) or 1.0
    return default[0] / n, default[1] / n


def build_domain(alert, cfg):
    """Return a domain dict (center_lat/lon, width_km, height_km) enlarged by the margins,
    with the centre shifted downstream and the upstream margin along the inflow."""
    c = polygon_centroid(alert.polygon)
    if c is None:                                          # no polygon -> fall back to geocode/area
        raise ValueError("alert has no polygon; cannot build a domain")
    lat0, lon0 = c
    ad = cfg.automatic_domain
    ex, ny = storm_motion(alert)
    # extent along motion = upstream + downstream; across = 2*lateral
    along = ad.upstream_margin_km + ad.downstream_margin_km
    across = 2.0 * ad.lateral_margin_km
    width_km = max(along, across)
    # shift the centre downstream by (downstream-upstream)/2 along the motion unit vector
    shift_km = 0.5 * (ad.downstream_margin_km - ad.upstream_margin_km)
    dlat = (shift_km * ny) / 111.0
    dlon = (shift_km * ex) / (111.0 * math.cos(math.radians(lat0)) or 1.0)
    return {"center_lat": lat0 + dlat, "center_lon": lon0 + dlon,
            "width_km": width_km, "height_km": ad.vertical_extent_km,
            "storm_motion_unit": (ex, ny)}


def nearest_radars(lat, lon, cfg):
    """NEXRAD stations within ``maximum_distance_km``, closest first, up to ``maximum_radars``."""
    d = sorted(((sid, _haversine_km(lat, lon, la, lo)) for sid, (la, lo) in NEXRAD_STATIONS.items()),
               key=lambda t: t[1])
    d = [(s, round(km, 1)) for s, km in d if km <= cfg.radar.maximum_distance_km]
    if not cfg.radar.use_multiple_radars:
        return d[:1]
    return d[:max(1, cfg.radar.maximum_radars)]


def nearest_metars(lat, lon, max_km=200.0, n=4):
    d = sorted(((sid, _haversine_km(lat, lon, la, lo)) for sid, (la, lo, _) in METAR_STATIONS.items()),
               key=lambda t: t[1])
    return [(s, round(km, 1)) for s, km in d if km <= max_km][:n]


def latest_hrrr_run(alert_time_iso):
    """The HRRR cycle (YYYYMMDD, HH) at/just before the alert time — the most recent hourly run
    that can cover the period.  ``alert_time_iso`` like '2013-05-20T20:12:00Z'."""
    import datetime as dt
    try:
        t = dt.datetime.fromisoformat(alert_time_iso.replace("Z", "+00:00"))
    except Exception:
        t = dt.datetime.utcnow()
    t = t.replace(minute=0, second=0, microsecond=0)
    return t.strftime("%Y%m%d"), t.strftime("%H")
