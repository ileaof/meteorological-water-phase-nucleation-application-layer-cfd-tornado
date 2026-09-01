# Continuity — `atmospheric_data` (real_case ingestion)

Resume‑after‑session‑loss notes for the real‑data module (ROADMAP §3a). Built on branch
`feat/atmospheric-data`, additive, **idealised mode untouched**.

## What is built & verified (on synthetic samples, no downloads)

| Layer | Module | Status |
|---|---|---|
| Unified format | `internal.AtmosphericState` (CF‑NetCDF (time,z,y,x), standard names, provenance) | ✅ tested |
| Thermo (SI) | `thermo` (ρ, T_v, θ, hydrostatic base, humidity) reusing engine constants | ✅ tested |
| Units / projection | `units`, `project` (pyproj Lambert + equirectangular fallback) | ✅ tested |
| Config / cache | `config.CaseConfig` (YAML), `cache.Cache` (offline‑aware) | ✅ tested |
| Sources | `synthetic` (full backbone), `sounding`/`metar`/`storm_events` (CSV/JSON) | ✅ tested |
| Sources (dep‑gated) | `hrrr` (cfgrib+AWS), `era5` (netCDF+cdsapi), `nexrad` (pyart+AWS) | ⚠️ guarded, not run here (deps absent) |
| Interpolation | `interpolate` (bilinear + linear/conservative‑z, clamp + error log) | ✅ tested |
| Base state / IC / BC | `basestate`, `ic_bc` (initial_conditions + boundary_{5} + surface_forcing) | ✅ tested |
| Radial operator | `radial` (V_r = V·r̂), `validation` (RMSE/MAE/bias/corr/FSS/CSI/displacement) | ✅ tested |
| QC | `qc` (physical checklist → JSON + Markdown) | ✅ tested |
| Driver | `driver` (preprocess/build/run/compare, CPU/GPU auto+fallback) | ✅ tested (CPU) |
| CLI | `cli`/`__main__` (6 subcommands, `--offline`) | ✅ tested |

Tests: `tests/test_atmospheric_data.py` — **15 passed, 2 skipped** (HRRR/NEXRAD dep‑gated).
Configs: `config/moore_2013.yaml`, `config/local_case.yaml`. Docs: `docs/REAL_CASE_DATA.md`
(+ HRRR/ERA5/NEXRAD/MOORE_2013).

## How to resume

1. `pip install cfgrib eccodes arm_pyart cdsapi pyproj` to light up the real readers (each is
   optional; the module already runs without them via synthetic fallback).
2. Verify: `python -m pytest tests/test_atmospheric_data.py -q` (expect 15 pass, 2→run when deps
   present) and `python -m atmospheric_data validate-input config/local_case.yaml --offline`.
3. Real Moore run: `python -m atmospheric_data preprocess config/moore_2013.yaml` (downloads
   HRRR + NEXRAD from AWS), then `run-case` / `compare-radar`.

## Remaining refinements (honest)

* **HRRR/ERA5 pressure→height**: readers hand back a pressure‑proxy vertical; use the GRIB/NC
  geopotential‑height to make the vertical map fully height‑correct per column (structure is in
  place; `interpolate.regrid_to_model` consumes height‑z).
* **Extra HRRR groups**: CAPE/CIN/SRH/reflectivity/fluxes live in separate GRIB messages — add
  `filter_by_keys` groups to `hrrr._to_state` (isobaric group already mapped).
* **Time‑dependent LBC**: the boundary files carry a time axis and `limited_area` applies the
  Davies zone; wiring per‑step target interpolation between analysis times into the parent loop
  is a thin addition (`driver.run_case` currently nudges toward the t0 environment).
* **Reflectivity forward operator**: a calibrated Z(microphysics) for the reflectivity CSI/FSS.
* **Non‑periodic parent**: `real_case` currently runs the periodic storm grid with a relaxation
  sponge; a true limited‑area wall parent is a follow‑up.
* **Moving nest / AMR from real IC**: connect the ingested parent to the `run_multilevel_nest`
  cascade (1.3 km → 444 m → 125 m) with the Galilean storm‑motion frame.

## Next task queued by the user

After this: **TKE‑1.5 (Deardorff)** closure — `storm_dynamics/turbulence.py` currently raises
`NotImplementedError` for `tke15`.
