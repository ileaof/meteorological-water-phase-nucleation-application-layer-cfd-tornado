"""storm-watch configuration (operational auto mode).

Loads the ``storm_watch:`` block of a YAML into validated dataclasses.  Everything the task
lists is configurable: event filters, severity, trigger levels, the observational‑confirmation
phrases (so an *alert* is never mistaken for a *confirmed tornado*), automatic‑domain margins,
download window, resource limits, and the safe automation defaults (download/preprocess on,
**simulate off**).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict


@dataclass
class RegionCfg:
    country: str = "US"
    states: list = field(default_factory=lambda: ["OK", "TX", "KS"])


@dataclass
class ActionsCfg:
    auto_download: bool = True
    auto_preprocess: bool = True
    auto_simulate: bool = False              # SAFE default: never auto-run the CFD


@dataclass
class DataCfg:
    hrrr: bool = True
    hrrr_subhourly: bool = True
    nexrad: bool = True
    metar: bool = True
    sounding: bool = True


@dataclass
class RadarCfg:
    maximum_distance_km: float = 250.0
    use_multiple_radars: bool = True
    maximum_radars: int = 3


@dataclass
class TimeWindowCfg:
    before_alert_minutes: int = 120
    after_expiration_minutes: int = 60
    radar_update_seconds: int = 300
    alert_poll_seconds: int = 60


@dataclass
class AutomaticDomainCfg:
    upstream_margin_km: float = 250.0        # into the storm's inflow (uses the motion vector)
    downstream_margin_km: float = 350.0
    lateral_margin_km: float = 200.0
    vertical_extent_km: float = 20.0


@dataclass
class SimulationTriggerCfg:
    minimum_level: str = "confirmed"         # watch | warning | confirmed
    minimum_severity: str = "severe"


@dataclass
class ResourceLimitsCfg:
    maximum_active_cases: int = 2
    maximum_storage_gb: float = 500.0
    maximum_case_storage_gb: float = 100.0
    maximum_concurrent_downloads: int = 4
    delete_raw_data: bool = False            # never auto-delete scientific data
    compression: bool = True


@dataclass
class NotificationsCfg:
    desktop: bool = True
    email: bool = False


# observational-confirmation phrases: only these upgrade an alert to `confirmed`
_CONFIRMED_PHRASES = ["radar confirmed tornado", "observed tornado", "tornado debris signature",
                      "considerable damage threat", "catastrophic damage threat",
                      "confirmed tornado", "tornado emergency"]


@dataclass
class StormWatchConfig:
    enabled: bool = True
    region: RegionCfg = field(default_factory=RegionCfg)
    events: list = field(default_factory=lambda: ["Tornado Warning", "Severe Thunderstorm Warning",
                                                  "Tornado Watch"])
    minimum_severity: str = "severe"
    alert_poll_seconds: int = 60
    actions: ActionsCfg = field(default_factory=ActionsCfg)
    data: DataCfg = field(default_factory=DataCfg)
    radar: RadarCfg = field(default_factory=RadarCfg)
    time_window: TimeWindowCfg = field(default_factory=TimeWindowCfg)
    automatic_domain: AutomaticDomainCfg = field(default_factory=AutomaticDomainCfg)
    simulation_trigger: SimulationTriggerCfg = field(default_factory=SimulationTriggerCfg)
    resource_limits: ResourceLimitsCfg = field(default_factory=ResourceLimitsCfg)
    notifications: NotificationsCfg = field(default_factory=NotificationsCfg)
    confirmed_phrases: list = field(default_factory=lambda: list(_CONFIRMED_PHRASES))
    workdir: str = "outputs/storm_watch"

    @classmethod
    def from_yaml(cls, path):
        import yaml
        with open(path, "r", encoding="utf-8") as f:
            d = (yaml.safe_load(f) or {}).get("storm_watch", {})
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, d):
        d = dict(d or {})
        nested = {"region": RegionCfg, "actions": ActionsCfg, "data": DataCfg, "radar": RadarCfg,
                  "time_window": TimeWindowCfg, "automatic_domain": AutomaticDomainCfg,
                  "simulation_trigger": SimulationTriggerCfg, "resource_limits": ResourceLimitsCfg,
                  "notifications": NotificationsCfg}
        kw = {}
        for k, C in nested.items():
            if isinstance(d.get(k), dict):
                fields = {fl.name for fl in C.__dataclass_fields__.values()}
                kw[k] = C(**{f: v for f, v in d[k].items() if f in fields})
        for k in ("enabled", "events", "minimum_severity", "alert_poll_seconds",
                  "confirmed_phrases", "workdir"):
            if k in d:
                kw[k] = d[k]
        return cls(**kw)

    def to_dict(self):
        return asdict(self)
