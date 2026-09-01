# Real atmospheric data — the `real_case` mode (`atmospheric_data`)

Initialise and validate the storm CFD with **real observations and analyses**, without
replacing the idealised mode. The `atmospheric_data` package downloads / reads data from six
sources, converts everything into one unified CF‑NetCDF format in SI units, builds the model's
initial + lateral boundary conditions, runs a quality‑control report, and validates the
simulated winds against Doppler radar in **radial‑velocity space**.

> The idealised mode (`storm_dynamics`, warm‑bubble supercell) is **untouched**. `real_case`
> is additive: `model.input_mode: idealized | real_case` in the case YAML.

## Sources

| # | Source | What it provides | Reader (optional dep) |
|---|---|---|---|
| A | **NOAA HRRR** | 3‑D T, p, q, wind, w, CAPE/CIN, reflectivity, terrain (GRIB2) | `cfgrib`+`eccodes`; download from AWS Open Data (no key) |
| B | **ERA5 / Copernicus** | global synoptic environment (NetCDF/GRIB) | `xarray` (read); `cdsapi`+`~/.cdsapirc` (download) |
| C | **NEXRAD Level II** | reflectivity, radial velocity, ZDR, ρhv, φdp (validation) | `arm_pyart` / `xradar`; download from AWS (no key) |
| D | **Radiosondes** | vertical profiles → base state + CAPE/shear/SRH | `pandas` (CSV/text) |
| E | **METAR/ASOS** | surface T, Td, p, wind (surface validation) | `pandas` (CSV) |
| F | **Storm Events / SWDI** | tornado track, EF rating, LSRs (case selection + track validation) | `pandas` (CSV/JSON); `geopandas` (shapefile) |

Only `numpy scipy xarray netCDF4 pandas pyyaml` are **required**; every heavy reader is
**optional** and degrades gracefully — if a library or the data is missing, the source is
skipped with a clear message and the pipeline falls back (finally to a labelled *synthetic*
sample), so the idealised mode never depends on any of them.

## Install

```bash
pip install numpy scipy xarray netCDF4 pandas pyyaml          # core (usually already present)
# optional, per source:
pip install cfgrib eccodes            # HRRR GRIB2      (conda-forge recommended for eccodes)
pip install cdsapi                    # ERA5 download   (+ ~/.cdsapirc credentials)
pip install arm_pyart                 # NEXRAD Level II
pip install pyproj                    # Lambert projection (else equirectangular fallback)
pip install geopandas shapely metpy   # Storm Events shapefile / extra diagnostics
```

## API credentials (never stored in the repo)

* **CDS API (ERA5):** create `~/.cdsapirc` per <https://cds.climate.copernicus.eu/how-to-api>.
* **HRRR / NEXRAD:** none — read from the AWS Open Data buckets `noaa-hrrr-bdp-pds` and
  `noaa-nexrad-level2`.

## Commands

```bash
python -m atmospheric_data case-info      config/moore_2013.yaml   # config + source availability
python -m atmospheric_data download       config/moore_2013.yaml   # fetch into the cache
python -m atmospheric_data preprocess     config/moore_2013.yaml   # -> IC/BC/surface + QC report
python -m atmospheric_data validate-input config/moore_2013.yaml   # QC only (exit != 0 on fail)
python -m atmospheric_data run-case       config/moore_2013.yaml   # start the CFD
python -m atmospheric_data compare-radar  config/moore_2013.yaml   # synthetic Vr vs NEXRAD
```

Add `--offline` to forbid all network access (uses only cached / local files), `--max-n N` to
cap grid points per axis (dev/tests; raise for production), `--steps N` for `run-case`.

### Offline / bring‑your‑own‑data

```bash
python -m atmospheric_data preprocess config/local_case.yaml --offline
```

`config/local_case.yaml` uses the synthetic sample (no network, no heavy deps). To run a real
case offline, set `atmospheric_source` to `hrrr`/`era5`, drop the downloaded files into
`data/cache/<source>/` (see [HRRR_IMPORT](HRRR_IMPORT.md) / [ERA5_IMPORT](ERA5_IMPORT.md) for
the expected filenames), then run `--offline`.

## Output files (per case, under `outputs/real_case/<name>/`)

```
initial_conditions.nc        full model-grid state at t0 (u,v,w,theta,qv,p)
boundary_{west,east,south,north,top}.nc   edge time series for the Davies relaxation zone
surface_forcing.nc           terrain (+ surface fluxes when available)
qc_report.{json,md}          quality-control report
radar_metrics.json           radial-velocity/reflectivity scores (compare-radar)
```

## Internal format, variables & units (SI)

One `xarray` dataset, dims `(time, z, y, x)`, standard names and per‑variable provenance
attributes (`units long_name source original_name interpolation_method valid_time projection
missing_value`):

`T` [K] · `p` [Pa] · `rho` [kg m⁻³] · `theta`,`theta_v` [K] · `qv,qc,qr,qi,qs,qg` [kg kg⁻¹] ·
`u,v,w` [m s⁻¹] · `terrain` [m] · `reflectivity` [dBZ]. Conversions applied to each field are
recorded (`interpolation_method`, e.g. `omega_to_w`, `theta(T,p)`).

## Model coupling (base state, IC, BC)

The anelastic core evolves perturbations about a hydrostatic base state. From the regridded
real fields we take the **horizontal‑mean** profile and re‑integrate `p0` hydrostatically
(`dp0/dz = -rho0 g`), then decompose `phi = phi0 + phi'`. Velocities are placed on the
staggered C‑grid; the first projection makes them anelastically divergence‑free. Lateral
boundaries use a **Davies relaxation zone** (`storm_dynamics.limited_area`) nudging the edges
toward the (optionally time‑dependent) environment, preventing reflections.

## Backends (CPU/GPU)

`model.execution_backend: auto` uses the GPU when present and **falls back to CPU
automatically**; `cpu`/`gpu` force a backend. CUDA is never required.

## Validation against radar (radial space)

A Doppler radar measures `V_r = V · r̂`, not the 3‑D vector. `compare-radar` interpolates the
CFD `(u,v,w)` to the radar gates and projects onto the beam, then scores against the observed
`V_r` (RMSE, MAE, bias, correlation) and reflectivity (CSI, FSS), and estimates the mesocyclone
displacement. See [NEXRAD_VALIDATION](NEXRAD_VALIDATION.md).

## Scientific limitations (must stay documented)

1. **HRRR does not resolve the tornado core** — it sets the storm‑scale environment.
2. **ERA5 is the synoptic environment**, not the tornadic vortex.
3. **A Doppler radar gives radial velocity**, not the full 3‑D wind (single‑radar inversion is
   under‑determined — hence validation in radial space).
4. **The reported track is observations / damage survey**, not a boundary condition.
5. **The tornado must emerge from the equations, resolution, supercell dynamics and physics** —
   it is not imposed.
6. **Δx ≈ 100–150 m** represents the tornadic circulation only in a limited way, not the full
   vortex sub‑structure.
7. **Direct radar assimilation** is a separate problem (3D/4D‑Var, EnKF) — this module does
   ingestion + a radial observation operator, not DA. Interpolation is not assimilation.

## References

* HRRR: <https://rapidrefresh.noaa.gov/hrrr/> · AWS: <https://registry.opendata.aws/noaa-hrrr-pds/>
* ERA5: <https://cds.climate.copernicus.eu/> · docs: <https://confluence.ecmwf.int/display/CKB/ERA5>
* NEXRAD Level II: <https://www.ncei.noaa.gov/products/radar/next-generation-weather-radar> ·
  AWS: <https://registry.opendata.aws/noaa-nexrad/> · Py‑ART: <https://arm-doe.github.io/pyart/>
* Storm Events / SWDI: <https://www.ncdc.noaa.gov/stormevents/> · <https://www.ncdc.noaa.gov/swdi/>
* Moore 2013: see [MOORE_2013_CASE](MOORE_2013_CASE.md).
