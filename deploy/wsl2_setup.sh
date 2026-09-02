#!/usr/bin/env bash
# Set up the FULL real-data stack on WSL2 (Ubuntu) via conda-forge, where cfgrib/eccodes and
# Py-ART are native.  Run inside WSL2:  bash deploy/wsl2_setup.sh
set -e
if ! command -v conda >/dev/null 2>&1; then
  echo "Install Miniforge first:  https://github.com/conda-forge/miniforge  (then re-run)"; exit 1
fi
ENV=met_real
conda create -y -n "$ENV" -c conda-forge python=3.11 \
    numpy scipy xarray netcdf4 h5netcdf pandas pyyaml matplotlib requests \
    cfgrib eccodes pyproj cdsapi metpy dask boto3 \
    arm_pyart xradar wradlib shapely geopandas
# nexradaws (AWS NEXRAD archive query) is pip-only:
conda run -n "$ENV" pip install nexradaws
echo
echo "Done.  Use it with:"
echo "  conda activate $ENV"
echo "  export PYTHONPATH=\$PWD/src"
echo "  python -m atmospheric_data preprocess    config/moore_2013_real.yaml"
echo "  python -m atmospheric_data run-case       config/moore_2013_real.yaml --multilevel"
echo "  python -m atmospheric_data compare-radar  config/moore_2013_real.yaml"
echo
echo "ERA5 needs ~/.cdsapirc (https://cds.climate.copernicus.eu/how-to-api) and the ERA5 licence"
echo "accepted; NEXRAD/HRRR need no credentials."
