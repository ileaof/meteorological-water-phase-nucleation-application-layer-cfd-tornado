"""Case configuration (YAML) for real-data runs (ROADMAP §3a).

Mirrors the task's YAML schema (``case / domain / data / model / processing / validation``) as
validated dataclasses with sensible defaults, so a real case is fully described by one file:

    python -m atmospheric_data preprocess config/moore_2013.yaml
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict


@dataclass
class CaseMeta:
    name: str = "unnamed_case"
    date: str = "2013-05-20"
    start_time_utc: str = "18:00"
    end_time_utc: str = "23:59"


@dataclass
class DomainCfg:
    center_lat: float = 35.34
    center_lon: float = -97.49
    width_km: float = 400.0
    height_km: float = 20.0
    projection: str = "lambert_conformal"


@dataclass
class DataCfg:
    atmospheric_source: str = "hrrr"        # hrrr | era5 | sounding | synthetic
    fallback_source: str = "era5"
    radar_station: str = "KTLX"
    use_sounding: bool = True
    use_metar: bool = True
    cache_directory: str = "data/cache"


@dataclass
class ModelCfg:
    execution_backend: str = "auto"         # auto | cpu | gpu
    input_mode: str = "real_case"           # idealized | real_case
    parent_dx_m: float = 1300.0
    nest_dx_m: float = 444.0
    fine_dx_m: float = 125.0
    moving_nest: bool = True


@dataclass
class ProcessingCfg:
    output_format: str = "netcdf4"
    interpolation: str = "conservative"     # conservative | linear
    hydrostatic_adjustment: bool = True
    temporal_interpolation: str = "linear"


@dataclass
class ValidationCfg:
    radar: bool = True
    surface_stations: bool = True
    tornado_track: bool = True


@dataclass
class CaseConfig:
    case: CaseMeta = field(default_factory=CaseMeta)
    domain: DomainCfg = field(default_factory=DomainCfg)
    data: DataCfg = field(default_factory=DataCfg)
    model: ModelCfg = field(default_factory=ModelCfg)
    processing: ProcessingCfg = field(default_factory=ProcessingCfg)
    validation: ValidationCfg = field(default_factory=ValidationCfg)
    offline: bool = False                    # set by --offline; no network access when True

    # ---- IO ----
    @classmethod
    def from_yaml(cls, path):
        import yaml
        with open(path, "r", encoding="utf-8") as f:
            d = yaml.safe_load(f) or {}
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, d):
        sub = lambda k, C: C(**{f: v for f, v in (d.get(k) or {}).items()
                                if f in {fl.name for fl in C.__dataclass_fields__.values()}})
        cfg = cls(case=sub("case", CaseMeta), domain=sub("domain", DomainCfg),
                  data=sub("data", DataCfg), model=sub("model", ModelCfg),
                  processing=sub("processing", ProcessingCfg),
                  validation=sub("validation", ValidationCfg))
        cfg.validate()
        return cfg

    def to_dict(self):
        return {k: asdict(getattr(self, k)) for k in
                ("case", "domain", "data", "model", "processing", "validation")}

    def validate(self):
        errs = []
        if not (-90 <= self.domain.center_lat <= 90):
            errs.append("domain.center_lat out of range")
        if not (-180 <= self.domain.center_lon <= 360):
            errs.append("domain.center_lon out of range")
        if self.domain.width_km <= 0 or self.domain.height_km <= 0:
            errs.append("domain width/height must be > 0")
        if self.model.input_mode not in ("idealized", "real_case"):
            errs.append("model.input_mode must be 'idealized' or 'real_case'")
        if self.model.execution_backend not in ("auto", "cpu", "gpu"):
            errs.append("model.execution_backend must be auto|cpu|gpu")
        if self.data.atmospheric_source not in ("hrrr", "era5", "sounding", "synthetic"):
            errs.append("data.atmospheric_source must be hrrr|era5|sounding|synthetic")
        if not (self.model.parent_dx_m > self.model.nest_dx_m > self.model.fine_dx_m > 0):
            errs.append("require parent_dx_m > nest_dx_m > fine_dx_m > 0")
        if errs:
            raise ValueError("invalid CaseConfig:\n  - " + "\n  - ".join(errs))
        return self
