# HRRR import (`atmospheric_data.sources.hrrr`)

NOAA's **High‑Resolution Rapid Refresh** (3 km CONUS) is the primary `atmospheric_source`. It
sets the storm‑scale **environment** — it does **not** resolve a tornado core (limitation 1).

## Dependencies

```bash
pip install cfgrib eccodes            # conda-forge strongly recommended for the ecCodes C lib
```

Absent → the source raises `SourceUnavailable` and the pipeline falls back (ERA5 → synthetic).

## Download (AWS Open Data, no credentials)

`download()` fetches the pressure‑level file for the case hour:

```
https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.<YYYYMMDD>/conus/hrrr.t<HH>z.wrfprsf00.grib2
```

into `data/cache/hrrr/<YYYYMMDD>_t<HH>z_wrfprs.grib2`. For **offline** use, place that file at
exactly that path and run with `--offline`.

## Variables read (SI, skipped if absent — never invented)

Isobaric group (`typeOfLevel=isobaricInhPa`): `t`→`T`, `q`→`qv`, `u`,`v`,`w`, plus pressure
from the level coordinate; `theta` is derived as `theta = T (P0/p)^κ`. HRRR also carries
surface/2 m/10 m groups and CAPE/CIN/SRH/reflectivity/terrain/fluxes in separate GRIB messages
(read via additional `filter_by_keys` — extend `_ISOBARIC`/add groups as needed).

## Vertical coordinate

HRRR is on pressure levels. The model vertical coordinate is **height**; the reader hands back
a native‑grid state whose vertical is a monotone pressure proxy, and `interpolate.regrid_to_model`
maps to the model heights. For a fully height‑correct map, convert pressure→height with the
GRIB geopotential‑height field (`gh`) before regridding (documented remaining refinement).

## Notes

* Different variable groups live in different GRIB level types — read each group with its own
  `filter_by_keys`; this keeps memory low and avoids cfgrib merge conflicts.
* The domain crop uses the projection (`pyproj` Lambert, else equirectangular fallback).
