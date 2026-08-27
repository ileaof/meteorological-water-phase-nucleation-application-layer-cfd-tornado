#!/bin/bash
# Turnkey pyAMReX build for WSL2/Ubuntu -- the reproducible recipe for the M3
# "full AMR" framework path (see docs/amr_design.md).  Every blocker found while
# bringing this up on a fresh WSL is fixed here:
#
#   * pyAMReX >= 26.8 needs Python >= 3.11  -> use a conda (Miniforge) Python 3.12
#     env (Ubuntu 22.04's system Python 3.10 is too old; the Miniforge *base* is
#     3.14, too new -- 3.12 is the sweet spot).
#   * the pybind11 bindings are memory-hungry: a parallel build OOM-kills cc1plus
#     on a 7 GB WSL.  Build with few jobs (JOBS below), or give WSL more RAM via
#     %UserProfile%\.wslconfig  ([wsl2]\nmemory=12GB) then `wsl --shutdown`.
#
# Usage (inside WSL):  bash scripts/build_pyamrex_wsl.sh
set -euo pipefail

JOBS="${JOBS:-2}"                 # keep low on <=8 GB WSL to avoid OOM
PREFIX="$HOME/miniforge3"
ENV="amr312"

# 1. Miniforge (no sudo) if absent
if [ ! -x "$PREFIX/bin/conda" ]; then
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/mf.sh
  bash /tmp/mf.sh -b -p "$PREFIX"
fi

# 2. Python 3.12 env with the build tools
"$PREFIX/bin/conda" env list | grep -q "$ENV" || \
  "$PREFIX/bin/conda" create -y -n "$ENV" python=3.12 numpy cmake ninja
PY="$PREFIX/envs/$ENV/bin/python"; PIP="$PREFIX/envs/$ENV/bin/pip"
echo "Using $($PY --version)"

# 3. Source (official pyAMReX)
[ -d "$HOME/pyamrex_src" ] || \
  git clone --depth 1 --recursive https://github.com/AMReX-Codes/pyamrex.git "$HOME/pyamrex_src"
cd "$HOME/pyamrex_src"; rm -rf build _skbuild

# 4. Build 3D CPU, low parallelism to fit memory
echo "Building pyAMReX (3D, CPU) with JOBS=$JOBS -- slow but memory-safe"
CMAKE_BUILD_PARALLEL_LEVEL="$JOBS" "$PIP" install -v \
  --config-settings=cmake.define.AMReX_SPACEDIM=3 \
  --config-settings=cmake.define.AMReX_GPU_BACKEND=NONE \
  .

# 5. Verify
"$PY" -c "import amrex.space3d as amr; print('AMReX OK', amr.__version__)"
echo "Done.  Run with:  $PY  (or: conda activate $ENV)"
