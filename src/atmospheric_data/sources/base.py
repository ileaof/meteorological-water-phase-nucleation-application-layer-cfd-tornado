"""Data-source plumbing: optional-dependency guards + a registry (ROADMAP §3a).

Heavy readers (GRIB via cfgrib, radar via Py-ART/xradar, ERA5 download via cdsapi) must never
break the idealized mode.  Each source declares the modules it needs; :func:`optional_import`
returns ``None`` when absent and :class:`SourceUnavailable` carries a clear, actionable message.
"""
from __future__ import annotations

import importlib


class SourceUnavailable(RuntimeError):
    """Raised when a source's optional dependency (or its data) is missing."""


def optional_import(module_name):
    """Import ``module_name`` or return ``None`` (no exception) if it is not installed."""
    try:
        return importlib.import_module(module_name)
    except Exception:
        return None


def require(module_name, source, hint=""):
    """Import ``module_name`` or raise :class:`SourceUnavailable` with an install hint."""
    mod = optional_import(module_name)
    if mod is None:
        raise SourceUnavailable(
            "%s needs the optional package '%s' (pip install %s). %s"
            % (source, module_name, module_name, hint))
    return mod


def available(module_names):
    """True iff every module in ``module_names`` imports."""
    return all(optional_import(m) is not None for m in module_names)
