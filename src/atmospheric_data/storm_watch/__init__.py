"""storm-watch: operational auto mode.

Continuously monitors official NWS/CAP alerts and, when a severe/tornadic storm is detected,
auto-selects a domain + radars, locates the HRRR run, downloads (or simulates the download of)
the data, validates it and generates the internal NetCDF — leaving the case READY_FOR_SIMULATION
(the CFD is auto-started only when explicitly enabled and the trigger level/severity are met).

Scientific safety: an alert triggers **data collection**, never the artificial insertion of a
tornado into the CFD — the tornado must emerge (or not) from the resolved dynamics. A Tornado
*Warning* is not a *confirmed* tornado; confirmation stages are preserved in the database.
"""
from __future__ import annotations

from .config import StormWatchConfig
from .db import WatchDB
from .monitor import StormWatchMonitor
from . import alerts, domain, casemachine

__all__ = ["StormWatchConfig", "WatchDB", "StormWatchMonitor", "alerts", "domain", "casemachine"]
