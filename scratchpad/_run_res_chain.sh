#!/usr/bin/env bash
cd "c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
# wait for whatever uniform run is active
while pgrep -f "uniform_resolution_test_gpu.py" >/dev/null 2>&1; do sleep 30; done
# HORIZONTAL lever: dx 300 m at the default coarse near-surface grid
RES=240 TAG=dx300 python scratchpad/uniform_resolution_test_gpu.py >> outputs/unires_chain.log 2>&1
# VERTICAL lever (the dominant one): dz1 5.5 m at fixed dx 600 m
RES=120 NZ=64 Z_STRETCH=1.09 TAG=dz5 python scratchpad/uniform_resolution_test_gpu.py >> outputs/unires_chain.log 2>&1
echo "CHAIN COMPLETE" >> outputs/unires_chain.log
