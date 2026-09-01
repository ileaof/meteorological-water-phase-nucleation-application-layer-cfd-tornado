"""Source registry + dispatch (ROADMAP §3a).

``load_atmosphere`` picks the configured atmospheric source, falling back to the fallback
source and finally to the synthetic sample if the primary is unavailable (offline / missing
dependency / no data) -- always with a clear log line, never a silent substitution.
"""
from __future__ import annotations

from . import synthetic
from .base import SourceUnavailable


def load_atmosphere(cfg, cache, logger=print):
    """Return an :class:`AtmosphericState` from the configured source, with fallback chain."""
    chain = [cfg.data.atmospheric_source, cfg.data.fallback_source, "synthetic"]
    seen = []
    last = None
    for name in chain:
        if name in seen:
            continue
        seen.append(name)
        try:
            if name == "hrrr":
                from . import hrrr; logger("[data] atmosphere: HRRR"); return hrrr.load(cfg, cache)
            if name == "era5":
                from . import era5; logger("[data] atmosphere: ERA5"); return era5.load(cfg, cache)
            if name == "synthetic":
                logger("[data] atmosphere: SYNTHETIC sample (no real data used)")
                return synthetic.synthetic_atmosphere(cfg)
            if name == "sounding":
                raise SourceUnavailable("'sounding' is a 1-D profile, not a gridded atmosphere; "
                                        "use it via use_sounding for the base state")
        except SourceUnavailable as e:
            logger("[data] %s unavailable (%s) -> trying next source" % (name, e))
            last = e
    raise SourceUnavailable("no atmospheric source available: %s" % last)


def load_radar(cfg, cache, logger=print):
    """Return a radar volume dict (real NEXRAD if available, else the synthetic sample)."""
    try:
        from . import nexrad
        if nexrad.available() and not cfg.offline:
            logger("[data] radar: NEXRAD Level II %s" % cfg.data.radar_station)
            return nexrad.load(cfg, cache)
    except SourceUnavailable as e:
        logger("[data] NEXRAD unavailable (%s) -> synthetic radar" % e)
    logger("[data] radar: SYNTHETIC sample")
    return synthetic.synthetic_radar(cfg)
