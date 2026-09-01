# ERA5 import (`atmospheric_data.sources.era5`)

ERA5 (ECMWF reanalysis, ~31 km) is the global **fallback** and the source of the **synoptic
environment**. It represents the large‑scale setting, **not** the tornadic vortex (limitation 2).

## Reading vs downloading

* **Reading** a NetCDF ERA5 file needs only `xarray` (always available).
* **Downloading** needs `cdsapi` and your CDS credentials in `~/.cdsapirc` (never in the repo):

```bash
pip install cdsapi
# ~/.cdsapirc  (from https://cds.climate.copernicus.eu/how-to-api)
url: https://cds.climate.copernicus.eu/api
key: <UID>:<APIKEY>
```

`download()` retrieves pressure‑level `temperature, specific_humidity, u/v/vertical_velocity,
geopotential` for the case time/area into `data/cache/era5/<date>_<HHMM>.nc`. Offline: place a
NetCDF there and run `--offline`.

## Variables read (SI)

`t`→`T`, `q`→`qv`, `u`,`v`; pressure from the level coordinate; `theta` derived. **Vertical
velocity** in ERA5 is pressure‑velocity ω [Pa s⁻¹]; converted to `w` [m s⁻¹] via
`w = -ω/(ρ g)` (recorded in `interpolation_method`). Geopotential `z` [m² s⁻²] → height via
`/g`.

## Notes

* The area is `[N, W, S, E]` around the domain centre, sized from `domain.width_km/height_km`.
* Latitude is usually descending in ERA5; the reader sorts the vertical to ascending pressure.
