#!/usr/bin/env bash
# Real Moore 2013 case end-to-end (ERA5 environment + NEXRAD KTLX validation).
# Prereqs: deploy/wsl2_setup.sh done; ~/.cdsapirc set; ERA5 licence accepted.
set -e
export PYTHONPATH="$PWD/src"
CFG=config/moore_2013_real.yaml
python -m atmospheric_data case-info      "$CFG"
python -m atmospheric_data download       "$CFG"          # real ERA5 (+ NEXRAD if nexradaws)
python -m atmospheric_data validate-input "$CFG"
python -m atmospheric_data run-case       "$CFG" --multilevel
python -m atmospheric_data compare-radar  "$CFG"          # synthetic Vr / reflectivity vs real KTLX
echo "Artefacts in outputs/real_case/moore_2013_real/"
