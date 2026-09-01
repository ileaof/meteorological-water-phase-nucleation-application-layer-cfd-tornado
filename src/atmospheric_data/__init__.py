"""Real atmospheric-data ingestion for the storm CFD model (ROADMAP §3a, "real_case" mode).

Downloads / reads real observations and analyses (HRRR, ERA5, NEXRAD Level II, radiosondes,
METAR/ASOS, Storm Events/SWDI), converts them into one unified CF-NetCDF internal format
(:mod:`atmospheric_data.internal`), and generates the model's initial + lateral boundary
conditions -- WITHOUT touching or replacing the idealized mode.  Heavy readers (GRIB, radar,
CDS) are optional: import them lazily and degrade gracefully when a library or the data is
absent (the idealized mode never depends on any of them).

CLI: ``python -m atmospheric_data <case-info|download|preprocess|validate-input|run-case|
compare-radar> config/<case>.yaml [--offline]``.
"""
from __future__ import annotations

from .config import CaseConfig
from .internal import AtmosphericState, STANDARD_VARS
from .cache import Cache
from . import thermo, units

__all__ = ["CaseConfig", "AtmosphericState", "STANDARD_VARS", "Cache", "thermo", "units"]
