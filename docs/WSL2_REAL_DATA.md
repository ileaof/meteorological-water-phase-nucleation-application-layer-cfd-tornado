# Real gridded data on WSL2 (HRRR/RAP GRIB + NEXRAD) — and ERA5 on Windows

The heavy readers (`cfgrib`/`eccodes` for GRIB, `arm_pyart` for NEXRAD) are **Linux-native**
and don't pip-install cleanly on Windows, but install fine on **WSL2** via conda-forge. **ERA5**
needs only `cdsapi` + a key and works on Windows directly.

## What each source needs (and where it runs)

| Source | Reader | Credentials | Windows | WSL2 |
|---|---|---|---|---|
| **ERA5** (2013 ✓) | `xarray` (read) + `cdsapi` (download) | CDS key `~/.cdsapirc` | ✅ | ✅ |
| **Radiosonde (IEM)** | built-in `iem_raob` | none | ✅ | ✅ |
| **HRRR / RAP** GRIB2 | `cfgrib`+`eccodes` | none (AWS) | ✗ (use WSL2) | ✅ |
| **NEXRAD Level II** | `arm_pyart` + `nexradaws` | none (AWS) | ✗ (use WSL2) | ✅ |

> **Note for Moore 2013:** HRRR did **not** exist yet (operational 2014‑09‑30). Use **ERA5**
> (covers 1940–present) for the gridded environment, and the real **KOUN sounding via IEM**
> (no credentials) as a proximity profile. RAP (2012+) is the GRIB alternative.

## ERA5 (works on Windows too)

1. Register + get a key: <https://cds.climate.copernicus.eu/how-to-api>.
2. **Accept the ERA5 licence** (downloads fail otherwise):
   <https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels> → *Terms of use*.
3. Create `~/.cdsapirc` (Windows: `C:\Users\<you>\.cdsapirc`):
   ```
   url: https://cds.climate.copernicus.eu/api
   key: <your-personal-access-token>
   ```
4. Run: `python -m atmospheric_data preprocess config/moore_2013_real.yaml` (downloads real ERA5).

## Full stack on WSL2

```bash
# inside WSL2 (Ubuntu), from the repo:
bash deploy/wsl2_setup.sh          # conda env 'met_real' with cfgrib/eccodes/pyart/nexradaws/...
conda activate met_real
export PYTHONPATH=$PWD/src
bash deploy/run_moore_real.sh      # ERA5 + NEXRAD KTLX + the full pipeline
```

`deploy/wsl2_setup.sh` builds a conda-forge env (`cfgrib eccodes pyproj cdsapi arm_pyart xradar
wradlib boto3 …` + `pip install nexradaws`).

## NEXRAD note

The `noaa-nexrad-level2` S3 bucket does **not** allow anonymous listing (GET-by-key only), so
the module uses **`nexradaws`** to query the archive and pick the scan nearest the case time
(no credentials). Reading the volume then uses **Py-ART**. Both are in the WSL2 env.

## Idealised mode is unaffected

None of this is required for the idealised mode or the offline synthetic pipeline — it only
lights up the real gridded/radar sources when present.
