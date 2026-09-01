"""NWS alert ingestion + classification (CAP v1.2 / GeoJSON).

Polls ``https://api.weather.gov/alerts/active`` (or reads a file for replay/offline), parses
every field the task lists, and classifies each alert into a trigger **level**
(``watch|warning|confirmed``).  Crucially, an alert only becomes ``confirmed`` when its text
contains an **observational‑confirmation phrase** (radar‑confirmed / observed / debris
signature / considerable|catastrophic damage) — a Tornado *Warning* alone is `warning`, never
`confirmed` (scientific-safety requirement).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

NWS_ACTIVE = "https://api.weather.gov/alerts/active"


@dataclass
class Alert:
    alert_id: str
    event: str = ""
    severity: str = "unknown"
    certainty: str = ""
    urgency: str = ""
    effective_time: str = ""
    onset_time: str = ""
    expiration_time: str = ""
    headline: str = ""
    description: str = ""
    instruction: str = ""
    affected_area: str = ""
    polygon: list = field(default_factory=list)     # [[lon,lat], ...]
    geocode: dict = field(default_factory=dict)
    issuing_office: str = ""

    def text(self):
        return " ".join([self.event, self.headline, self.description, self.instruction]).lower()


def parse_feature(feat):
    """Parse one GeoJSON alert feature into an :class:`Alert`."""
    p = feat.get("properties", {}) if isinstance(feat, dict) else {}
    geom = feat.get("geometry") or {}
    poly = []
    if geom.get("type") == "Polygon" and geom.get("coordinates"):
        poly = geom["coordinates"][0]
    return Alert(
        alert_id=p.get("id") or feat.get("id") or "",
        event=p.get("event", ""), severity=(p.get("severity") or "unknown").lower(),
        certainty=p.get("certainty", ""), urgency=p.get("urgency", ""),
        effective_time=p.get("effective", ""), onset_time=p.get("onset", ""),
        expiration_time=p.get("expires", ""), headline=p.get("headline", ""),
        description=p.get("description", ""), instruction=p.get("instruction") or "",
        affected_area=p.get("areaDesc", ""), polygon=poly,
        geocode=p.get("geocode", {}) or {}, issuing_office=p.get("senderName", ""))


def parse_feature_collection(data):
    feats = data.get("features", []) if isinstance(data, dict) else []
    return [parse_feature(f) for f in feats]


def load_alerts_file(path):
    """Read a saved GeoJSON alert (single feature or FeatureCollection) — used for replay."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if data.get("type") == "FeatureCollection":
        return parse_feature_collection(data)
    if data.get("type") == "Feature":
        return [parse_feature(data)]
    if "features" in data:
        return parse_feature_collection(data)
    return [parse_feature(data)]


def fetch_active_alerts(cfg, offline=False, timeout=30):
    """Fetch active alerts from the NWS API (filtered by region states).  Offline → []."""
    if offline:
        return []
    import requests
    params = {}
    if cfg.region.states:
        params["area"] = ",".join(cfg.region.states)
    r = requests.get(NWS_ACTIVE, params=params, timeout=timeout,
                     headers={"User-Agent": "met_h2o-storm-watch (research)", "Accept": "application/geo+json"})
    r.raise_for_status()
    return parse_feature_collection(r.json())


def _severity_ok(alert, minimum):
    rank = {"extreme": 4, "severe": 3, "moderate": 2, "minor": 1, "unknown": 0}
    return rank.get(alert.severity, 0) >= rank.get(str(minimum).lower(), 0)


def event_matches(alert, cfg):
    """True if the alert's event is in the configured filter list."""
    return any(e.lower() in alert.event.lower() for e in cfg.events)


def confirmation_stage(alert):
    """The scientific confirmation stage kept in the DB (never overstated)."""
    t = alert.text()
    if "damage" in t and ("catastrophic" in t or "considerable" in t):
        return "damage-confirmed"
    if any(k in t for k in ("radar confirmed tornado", "tornado debris signature")):
        return "radar-indicated"
    if "observed tornado" in t or "confirmed tornado" in t or "tornado emergency" in t:
        return "observed"
    return "alert-issued"


def classify_level(alert, cfg):
    """Return ``'watch' | 'warning' | 'confirmed'`` or ``None`` (filtered out)."""
    if not event_matches(alert, cfg) or not _severity_ok(alert, cfg.minimum_severity):
        return None
    t = alert.text()
    confirmed = any(ph.lower() in t for ph in cfg.confirmed_phrases)
    ev = alert.event.lower()
    if "warning" in ev:
        return "confirmed" if confirmed else "warning"
    if "watch" in ev:
        return "confirmed" if confirmed else "watch"       # rare, but honour explicit evidence
    return "warning" if confirmed else "watch"
