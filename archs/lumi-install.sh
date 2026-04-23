#!/usr/bin/env bash

# Runko installation script for on the LUMI supercomputer.
# Comment out unnecessary sections as you see fit.

set -v

# 1. Change current directory to the runko repository root
cd "$(dirname "$0")"
cd ../

RUNKO_PATH=$(pwd)

# 2. Initialize git submodules and download rocThrust headers.
# rocThrust is needed for both GPU and CPU backends (tyvi uses thrust API).
# GPU builds find rocthrust via the rocm module; CPU builds use a minimal
# cmake config (cmake/rocthrust-cpu/) that points to these downloaded headers.
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

# 3. Load standard prerequisite modules for runko:
module load LUMI/25.09
module load partition/G
module load PrgEnv-cray
module load rocm/6.4.4
module load cray-hdf5
module load craype-accel-amd-gfx90a
module load cray-mpich craype-network-ofi
module load buildtools
module load lumi-CrayPath


# 4. Create a Python virtual environment specially for runko,
#    located within the runko repository:
module load cray-python
python -m venv runko-venv

# 5. Load this runko virtual environment:
source runko-venv/bin/activate

# 6. Update the PYTHONPATH environment variable with required runko
#    and corgi-related paths in order to make our venv runko-aware::
P1="$RUNKO_PATH/"
P2="$RUNKO_PATH/external/corgi/lib"
export PYTHONPATH="$PYTHONPATH:$P1:$P2"

# 7. Install necessary Python dependencies
pip3 install h5py scipy matplotlib numpy

# 8. Build mpi4py against Cray MPICH (PrgEnv-cray).
# --force-reinstall ensures the venv gets its own copy instead of using the
# system cray-python mpi4py, which is built against GNU MPICH.
MPI4PY_BUILD_MPICC="cc -shared" pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py

# 9. Build runko (on login node; no GPU needed for compilation):
cmake --preset lumi-gpu-release
cmake --build lumi-gpu-release -j18

# 10. Patch the activate script with module loads, PYTHONPATH, and GTL preload
#     so that "source runko-venv/bin/activate" sets up the full environment:

# 11. Finally, we create a handy file which loads the runko virtual environment and
# necessary modules whenever called with "source runko-venv/bin/activate":
cat >> runko-venv/bin/activate << EOL
# Tool to load runko modules
# Usage: "source runko-venv/bin/activate"
# Load standard prerequisite modules for runko:
module load LUMI/25.09
module load partition/G
module load PrgEnv-cray
module load rocm/6.4.4
module load cray-hdf5
module load craype-accel-amd-gfx90a
module load cray-mpich craype-network-ofi
module load buildtools
module load lumi-CrayPath



# module load cray-python # not necessary as we are using a python virtual environment already

# GPU-aware MPI: preload GTL so it is available before mpi4py calls MPI_Init
export LD_PRELOAD=\${CRAY_MPICH_ROOTDIR}/gtl/lib/libmpi_gtl_hsa.so

export PYTHONPATH="\$PYTHONPATH:${P1}:${P2}"

EOL

#--------------------------------------------------
# 12. Building the CPU version (optional).
#
# The activate script loads craype-accel-amd-gfx90a and rocm for GPU builds.
# For CPU builds, these must be unloaded first — otherwise the Cray CC wrapper
# compiles everything as HIP/GPU code.
#
#   source runko-venv/bin/activate
#   module unload craype-accel-amd-gfx90a rocm
#   cmake --preset lumi-cpu-release
#   cmake --build lumi-cpu-release -j18
#

#--------------------------------------------------
# 13. Run tests on a GPU compute node.
#
# The login node has no GPUs, so tests must be run on a compute node.
# Get an interactive GPU allocation:
#
#   srun --account=<account> --partition=standard-g --gpus-per-node=1 --ntasks=1 --cpus-per-task=6 --mem=8G --time=00:30:00 --pty bash
#
# Then activate the environment and run tests:
#
#   source runko-venv/bin/activate
#   python -m unittest discover -s tests/ -v
#   cd lumi-gpu-release && ctest -V
#
