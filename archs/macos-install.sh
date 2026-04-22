#!/usr/bin/env bash

# Runko installation script for a local macOS (Apple Silicon / M-series) machine.
# CPU build via the `macos-mpic++-vec-analysis` preset.
#
# Prerequisites (install with Homebrew first):
#   brew install open-mpi hdf5 python@3 cmake git
#
# This script runs end-to-end: it sets up submodules, downloads dependencies,
# creates a Python venv, installs Python packages, and builds runko.

set -ve

#--------------------------------------------------
# Change current directory to the runko repository root
cd "$(dirname "$0")"
cd ../

RUNKO_PATH=$(pwd)

#--------------------------------------------------
# Initialize git submodules and download rocThrust headers.
# rocThrust is needed for the CPU backend (tyvi uses the thrust API);
# CPU builds resolve it via the minimal cmake config (cmake/rocthrust-cpu/)
# that points to these downloaded headers.
git submodule update --init --recursive

ROCM_LIBS_DIR="$RUNKO_PATH/external/tyvi/rocm-libraries"
if [ ! -d "$ROCM_LIBS_DIR/projects/rocthrust/thrust" ]; then
    cd "$RUNKO_PATH/external/tyvi"
    git clone --no-checkout --depth=1 --filter=tree:0 https://github.com/ROCm/rocm-libraries.git
    cd rocm-libraries
    git sparse-checkout init --cone
    git sparse-checkout set projects/rocthrust
    git checkout develop
    cd "$RUNKO_PATH"
fi

#--------------------------------------------------
# Pre-download kokkos/mdspan into the repo tree; passed to CMake as
# -DCPM_mdspan_SOURCE so CPM does not re-fetch on every reconfigure.
MDSPAN_PATH="$RUNKO_PATH/external/mdspan"
if [ ! -d "$MDSPAN_PATH" ]; then
    git clone https://github.com/kokkos/mdspan.git "$MDSPAN_PATH"
fi

cd "$RUNKO_PATH"

#--------------------------------------------------
# Create Python virtual environment for runko CPU builds
python3 -m venv venv/runko-cpu

# Activate the virtual environment
source venv/runko-cpu/bin/activate

# Install necessary Python dependencies
pip install --upgrade pip
pip install h5py scipy matplotlib numpy

# Build mpi4py against the Homebrew MPI that `mpicc` resolves to.
# --force-reinstall --no-cache-dir rebuilds mpi4py cleanly when the script
# is re-run (e.g. after a Homebrew MPI upgrade).
MPICC=mpicc pip install --no-cache-dir --no-binary=mpi4py mpi4py

#--------------------------------------------------
# Patch the activate script with PYTHONPATH so pycorgi / runko_cpp_bindings
# are importable after `source venv/runko-cpu/bin/activate`.
P1="$RUNKO_PATH/"
P2="$RUNKO_PATH/external/corgi/lib"

cat >> venv/runko-cpu/bin/activate << EOL

# --- Runko environment (macOS CPU) ---
# Usage: source venv/runko-cpu/bin/activate
export PYTHONPATH="\$PYTHONPATH:${P1}:${P2}"
EOL

#--------------------------------------------------
# Configure and build with the macOS vectorization-analysis preset.
cmake --preset macos-mpic++-vec-analysis
cmake --build macos-mpic++-vec-analysis -j8

#--------------------------------------------------
# Run the multirank (MPI) test suite.
cd "$RUNKO_PATH/macos-mpic++-vec-analysis"
ctest -V

# Build complete: $RUNKO_PATH/macos-mpic++-vec-analysis/
# Re-enter the environment with: source venv/runko-cpu/bin/activate
