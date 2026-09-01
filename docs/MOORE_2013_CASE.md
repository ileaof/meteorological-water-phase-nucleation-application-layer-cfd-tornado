# Moore, Oklahoma — 20 May 2013 (EF5)

The documented historical case shipped in `config/moore_2013.yaml`. The architecture is
**general** — change `case`/`domain`/`data.radar_station` for any other event.

## The event (context, not model input)

* **Date/time:** 2013‑05‑20, tornado ~19:56–20:35 UTC (touchdown near Newcastle, OK).
* **Rating:** EF5; ~14 mile (≈22 km) track, up to ~1.1 mile wide, through Moore, OK.
* **Radar:** **KTLX** (Twin Lakes, OK) — Moore is well within range; **KOUN** soundings nearby.
* **Environment:** strongly unstable, sheared warm‑sector supercell environment (large CAPE,
  strong 0–6 km shear and 0–3 km SRH).

These facts (from Storm Events / damage surveys) are used to **select** the case and to
**validate** the simulated track/timing — never to build the flow fields.

## Configuration

```yaml
case:   { name: moore_2013, date: "2013-05-20", start_time_utc: "18:00", end_time_utc: "23:59" }
domain: { center_lat: 35.34, center_lon: -97.49, width_km: 400, height_km: 20 }
data:   { atmospheric_source: hrrr, fallback_source: era5, radar_station: KTLX }
model:  { input_mode: real_case, parent_dx_m: 1300, nest_dx_m: 444, fine_dx_m: 125, moving_nest: true }
```

## Run

```bash
python -m atmospheric_data case-info      config/moore_2013.yaml
python -m atmospheric_data download       config/moore_2013.yaml        # HRRR (+NEXRAD) from AWS
python -m atmospheric_data preprocess     config/moore_2013.yaml        # -> IC/BC/surface + QC
python -m atmospheric_data validate-input config/moore_2013.yaml
python -m atmospheric_data run-case       config/moore_2013.yaml --steps 200
python -m atmospheric_data compare-radar  config/moore_2013.yaml
```

Without `cfgrib`/`arm_pyart` installed, the commands run against a labelled **synthetic**
environment/radar (so the workflow is demonstrable offline); install those to use the real
HRRR/NEXRAD data.

## Expectations & honesty

The nested cascade (1.3 km → 444 m → ~125 m) resolves the **supercell and its low‑level
mesocyclone**; a resolved condensation funnel needs O(10 m). A right‑moving supercell with a
strengthening low‑level vortex — matched against KTLX radial velocity and the reported track —
is the realistic target. The tornado itself must **emerge** from the dynamics; it is not
imposed (limitation 5).
