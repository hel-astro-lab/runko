#!/usr/bin/env bash

# Runko installation script for on the LUMI supercomputer.
# Comment out unnecessary sections as you see fit.

set -v

# 1. Change current directory to the runko repository root
cd "$(dirname "$0")"
cd ../

# 3. Load standard prerequisite modules for runko:
module load LUMI/24.03
module load partition/G
module load PrgEnv-cray
module load rocm
module load cray-hdf5
module load craype-accel-amd-gfx90a
module load cray-mpich craype-network-ofi
module load buildtools

# 4. Create a Python virtual environment specially for runko,
#    located within the runko repository:
module load cray-python
python -m venv runko-venv

# 5. Load this runko virtual environment:
source runko-venv/bin/activate

# 6. Update the PYTHONPATH environment variable with required runko
#    and corgi-related paths in order to make our venv runko-aware::
RUNKO_PATH=$(pwd)
P1="$RUNKO_PATH/"
P2="$RUNKO_PATH/lib"
P3="$RUNKO_PATH/external/corgi/lib"
export PYTHONPATH="$PYTHONPATH:$P1:$P2:$P3"

# 7. Install necessary Python dependencies
pip3 install h5py scipy matplotlib numpy

# 8. Build MPI4PY:
MPI4PY_BUILD_MPICC="cc -shared" python -m pip install --no-binary=mpi4py mpi4py

# 9. Build and unit test runko on 16 cores:
cmake -Bbuild -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=CC .
make -j16 -Cbuild

# 10. If all went well runko has now been compiled, and certain unit tests
#     should have automatically been executed and passed.

# 11. Finally, we create a handy file which loads the necessary runko modules
# whenever called with "source archs/runko-load-env":
cat > archs/lumi-load-runko-env << EOL
# Tool to load runko modules
# Usage: "source archs/runko-load-env"
# Load standard prerequisite modules for runko:
module load LUMI/24.03
module load partition/G
module load PrgEnv-cray
module load rocm
module load cray-hdf5
module load craype-accel-amd-gfx90a
module load cray-mpich craype-network-ofi
module load buildtools
module load cray-python
# Load runko virtual environment:
source ${RUNKO_PATH}/runko-venv/bin/activate
# Update the PYTHONPATH environment variable with runko path data
export PYTHONPATH="\$PYTHONPATH:${P1}:${P2}:${P3}"

EOL
