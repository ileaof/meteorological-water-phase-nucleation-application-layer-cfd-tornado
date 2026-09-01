# storm-watch — operational auto mode

Continuously monitors official **NWS/CAP alerts** and, when a severe/tornadic storm is
detected, automatically builds a domain, selects radars, locates the HRRR run, downloads (or
simulates) the data, validates it, and generates the internal NetCDF — leaving the case
**ready for simulation**. Part of `atmospheric_data`; the idealised mode is unaffected.

> **Scientific safety.** An alert triggers **data collection**, *never* the artificial
> insertion of a tornado into the CFD. The tornado must emerge (or not) from the resolved
> dynamics. A Tornado *Warning* is **not** a confirmed tornado — confirmation stages
> (`alert-issued → radar-indicated → observed → damage-confirmed`) are preserved in the DB.

## Quick start (offline replay — the full automation, no network)

```bash
python -m atmospheric_data storm-watch replay tests/data/sample_tornado_alert.json config/storm_watch.yaml
python -m atmospheric_data storm-watch cases  config/storm_watch.yaml
python -m atmospheric_data storm-watch alerts config/storm_watch.yaml
```

## Realtime

```bash
python -m atmospheric_data storm-watch start  config/storm_watch.yaml           # foreground loop
python -m atmospheric_data storm-watch start  config/storm_watch.yaml --max-iterations 10
python -m atmospheric_data storm-watch status config/storm_watch.yaml
python -m atmospheric_data storm-watch stop   config/storm_watch.yaml           # writes a STOP flag
python -m atmospheric_data storm-watch retry  CASE_ID config/storm_watch.yaml
```

## Alert monitoring

Polls `https://api.weather.gov/alerts/active` (GeoJSON/CAP v1.2), filtered by region states,
and extracts: `alert_id event severity certainty urgency effective/onset/expiration_time
headline description instruction affected_area polygon geocode issuing_office`. Duplicates are
never reprocessed (their id is stored in SQLite).

### Trigger levels

* **watch** — register + light environmental download; **no** full simulation.
* **warning** — create a case; download HRRR + radar + obs; preprocess.
* **confirmed** — prioritise; faster updates; consecutive radar volumes; auto-simulate allowed.

An alert becomes **confirmed only** on an observational-confirmation phrase (configurable):
`radar confirmed tornado`, `observed tornado`, `tornado debris signature`,
`considerable/catastrophic damage threat`.

## Automatic domain & assets

From the polygon: centroid → domain enlarged by `upstream/downstream/lateral` margins (the
upstream margin follows the **storm-motion vector**, parsed from "moving NE at N mph"); nearest
**NEXRAD** radars (distance + count limited); nearest **METAR/ASOS**; the most recent **HRRR**
cycle covering the period.

## Case state machine

```
DETECTED → DOWNLOADING → READY_FOR_PREPROCESSING → PREPROCESSING → READY_FOR_SIMULATION
        → (SIMULATING → VALIDATING → COMPLETED)   # only if auto_simulate AND the trigger is met
        → FAILED                                   # monitor retries with backoff
```

Every transition, failure and retry is logged in `outputs/storm_watch/storm_watch.sqlite`.

## Automation defaults (safe)

```yaml
actions: { auto_download: true, auto_preprocess: true, auto_simulate: false }
simulation_trigger: { minimum_level: confirmed, minimum_severity: severe }
```

Auto-simulation is **off** unless explicitly enabled *and* the level/severity gate is met.

## Resilience

Timeout + exponential backoff on API/HRRR/radar failure (a temporary outage never kills the
monitor); SQLite dedup + severity-priority queue that persists across restarts; resource limits
(`maximum_active_cases`, storage caps, concurrent downloads). If the current HRRR is not yet
posted: wait → try an earlier run → RAP → (ERA5 only for later post-processing, not as an
immediate operational substitute). **Scientific data is never auto-deleted** (`delete_raw_data:
false`); on hitting a storage cap the monitor suspends new downloads, preserves metadata, warns.

## Realtime vs historical

* **realtime** — monitor active alerts + current data.
* **historical / replay** — supply a date/time/region (or a saved alert file); no active-alert
  polling. `storm-watch replay FILE` runs the entire automation offline.

## Deployment (run continuously)

* **Linux (systemd):** `deploy/storm-watch.service` (`systemctl --user enable --now storm-watch`).
* **Windows (Task Scheduler / at logon):** `deploy/storm-watch-task.ps1` registers a task that
  starts the monitor at logon; or run it as a foreground/background process.
* **Docker:** any Python 3.11+ image with the repo installed; `CMD python -m atmospheric_data
  storm-watch start config/storm_watch.yaml`.

## Configuration

See `config/storm_watch.yaml` (the `storm_watch:` block) — events, severity, poll seconds,
actions, data toggles, radar limits, automatic-domain margins, time window, resource limits,
notifications, and the confirmation phrases.

## Tests

`tests/test_storm_watch.py` (9) covers: new/duplicate/filtered alerts, missing polygon, API
failure, two simultaneous alerts, resource-limit queueing, restart persistence, level
classification, auto-domain/radar/HRRR selection, and the offline replay reaching
`READY_FOR_SIMULATION` with the internal NetCDF written and logged.
