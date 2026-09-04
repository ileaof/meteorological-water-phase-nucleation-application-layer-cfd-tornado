#!/usr/bin/env bash
cd "c:/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
for s in 1 2 3; do
  RES=120 SEED=$s TAG=A_seed$s python scratchpad/uniform_resolution_test_gpu.py >> outputs/seeds.log 2>&1
done
echo "SEEDS COMPLETE" >> outputs/seeds.log
