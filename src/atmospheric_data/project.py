"""Horizontal map projection lat/lon <-> model x,y metres (ROADMAP §3a).

Uses ``pyproj`` (Lambert conformal, the HRRR/WRF standard) when installed; otherwise falls
back to a documented **local tangent-plane / equirectangular** approximation about the domain
centre -- adequate for a ~400 km storm domain, and the fallback is recorded in metadata so it
is never silent (task requirement 5/7).
"""
from __future__ import annotations

import importlib

import numpy as np

_R_EARTH = 6_371_000.0


def _has_pyproj():
    try:
        importlib.import_module("pyproj")
        return True
    except Exception:
        return False


class Projection:
    """Map between (lat, lon) [deg] and projected (x, y) [m] about a centre.

    ``kind``: 'lambert_conformal' (needs pyproj) or 'equirectangular' (fallback, always
    available).  ``.method`` records which was actually used, for provenance."""

    def __init__(self, center_lat, center_lon, kind="lambert_conformal"):
        self.lat0 = float(center_lat); self.lon0 = float(center_lon)
        self.kind = kind
        self._t = None
        if kind == "lambert_conformal" and _has_pyproj():
            import pyproj
            proj = pyproj.Proj(proj="lcc", lat_1=self.lat0 - 5, lat_2=self.lat0 + 5,
                               lat_0=self.lat0, lon_0=self.lon0, R=_R_EARTH)
            self._t = proj
            self.method = "pyproj_lambert_conformal"
        else:
            self.method = "equirectangular_fallback" if kind == "lambert_conformal" \
                else "equirectangular"

    def to_xy(self, lat, lon):
        lat = np.asarray(lat, float); lon = np.asarray(lon, float)
        if self._t is not None:
            return self._t(lon, lat)                       # (x, y) metres
        x = np.radians(lon - self.lon0) * _R_EARTH * np.cos(np.radians(self.lat0))
        y = np.radians(lat - self.lat0) * _R_EARTH
        return x, y

    def to_lonlat(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        if self._t is not None:
            lon, lat = self._t(x, y, inverse=True)
            return lat, lon
        lat = self.lat0 + np.degrees(y / _R_EARTH)
        lon = self.lon0 + np.degrees(x / (_R_EARTH * np.cos(np.radians(self.lat0))))
        return lat, lon
