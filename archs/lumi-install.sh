#!/usr/bin/env bash

# Runko installation script for LUMI. Automates runko_lumi_v4.md.
#
# Usage:  ./archs/lumi-install.sh gpu
#         ./archs/lumi-install.sh cpu
#
# Run each from a fresh login shell: the two module stacks conflict.
# Both backends compile on the login node; tests need an allocation.

#--------------------------------------------------
# Settings

LUMI_PROJECT="project_462001358"
SCRATCH_ROOT="/pfs/lustrep3/scratch/${LUMI_PROJECT}/${USER}"

#--------------------------------------------------

# Refuse to be sourced: set -e would kill the calling shell.
if [ "${BASH_SOURCE[0]}" != "${0}" ]; then
    echo "error: run this script, do not source it: ./archs/lumi-install.sh {gpu|cpu}" >&2
    return 1
fi

set -e

usage() {
    echo "usage: $0 {gpu|cpu}" >&2
}

BACKEND="${1:-}"
case "$BACKEND" in
    gpu) PRESET="lumi-gpu"; VENV_NAME="runko-g"; MODULES="runko-modules-g.sh" ;;
    cpu) PRESET="lumi-cpu"; VENV_NAME="runko-c"; MODULES="runko-modules-c.sh" ;;
    *)   usage; exit 1 ;;
esac

# The clone this script lives in; build dir name matches the preset name.
RUNKO_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
VENV="$SCRATCH_ROOT/venvs/$VENV_NAME"
BUILD_DIR="$RUNKO_PATH/build/$PRESET"
ROCTHRUST_PREFIX="$SCRATCH_ROOT/rocm-libraries/projects/rocthrust/rocthrust-install"

if [ ! -f "$RUNKO_PATH/CMakePresets.json" ]; then
    echo "error: $RUNKO_PATH is not a runko checkout" >&2
    exit 1
fi

if [ "$RUNKO_PATH" != "$SCRATCH_ROOT/runko" ]; then
    echo "warning: repo is $RUNKO_PATH, not $SCRATCH_ROOT/runko" >&2
fi

set -v

#--------------------------------------------------
# Submodules

git -C "$RUNKO_PATH" submodule update --init --recursive

#--------------------------------------------------
# Module stacks

mkdir -p "$SCRATCH_ROOT/venvs"

cat > "$SCRATCH_ROOT/runko-modules-g.sh" << 'EOF'
module load LUMI/25.09
module load partition/G
module load PrgEnv-cray
module load rocm/6.4.4
module load craype-accel-amd-gfx90a
module load cray-mpich/9.0.1
module load craype-network-ofi
module load buildtools
module load cray-python
module load lumi-CrayPath
EOF

cat > "$SCRATCH_ROOT/runko-modules-c.sh" << 'EOF'
module load LUMI/25.09
module load partition/C
module load PrgEnv-cray
module load cray-mpich/9.0.1
module load craype-network-ofi
module load buildtools
module load cray-python
module load lumi-CrayPath
EOF

# Lmod returns nonzero in ordinary situations.
set +e
source "$SCRATCH_ROOT/$MODULES"
set -e

#--------------------------------------------------
# Virtual environment

if [ ! -d "$VENV" ]; then
    python -m venv "$VENV"
fi

# RUNKODIR and SKBUILD_BUILD_DIR give each backend its own build directory.
if ! grep -q "runko $BACKEND backend" "$VENV/bin/activate"; then
    cat >> "$VENV/bin/activate" << EOF

# runko $BACKEND backend
export RUNKODIR=$RUNKO_PATH
export SKBUILD_BUILD_DIR=\$RUNKODIR/build/$PRESET
EOF
    if [ "$BACKEND" = "cpu" ]; then
        echo "export ROCTHRUST_INSTALL_PREFIX=$ROCTHRUST_PREFIX" >> "$VENV/bin/activate"
    fi
fi

source "$VENV/bin/activate"

#--------------------------------------------------
# Dependencies

pip install --upgrade pip

# Link mpi4py against Cray MPICH; the PyPI wheel fails at MPI_Init.
MPI4PY_BUILD_MPICC=cc python -m pip install --no-cache-dir --no-binary=mpi4py mpi4py --force-reinstall

pip install pybind11 scikit-build-core numpy scipy matplotlib

#--------------------------------------------------
# rocThrust, CPU backend only

if [ "$BACKEND" = "cpu" ]; then
    if [ ! -d "$ROCTHRUST_PREFIX/lib/cmake/rocthrust" ]; then
        cd "$SCRATCH_ROOT"
        if [ ! -d rocm-libraries ]; then
            git clone --no-checkout --depth=1 --filter=tree:0 https://github.com/ROCm/rocm-libraries.git
            cd rocm-libraries
            git sparse-checkout init --cone
            git sparse-checkout set projects/rocthrust
            git checkout develop
            cd "$SCRATCH_ROOT"
        fi
        cd rocm-libraries/projects/rocthrust
        cmake -Bbuild -DTHRUST_DEVICE_SYSTEM=CPP -DLINK_HIP_DEVICE_LIBS=OFF -DCMAKE_INSTALL_PREFIX="$ROCTHRUST_PREFIX" .
        make -C build install
        cd "$SCRATCH_ROOT"
    fi

    # The lumi-cpu preset resolves rocthrust_DIR from here.
    if [ ! -d "$ROCTHRUST_PREFIX/lib/cmake/rocthrust" ]; then
        set +v
        echo "error: rocThrust config missing at $ROCTHRUST_PREFIX/lib/cmake/rocthrust" >&2
        exit 1
    fi
fi

#--------------------------------------------------
# Build and install runko

pip install --no-build-isolation --no-deps -v -e "$RUNKO_PATH" \
    --config-settings=cmake.args=--preset="$PRESET"

#--------------------------------------------------
# Next steps

set +v

if [ "$BACKEND" = "gpu" ]; then
    cat << EOF

Installed into $VENV; objects in $BUILD_DIR

Login nodes have no GPU. Test on a compute node:

  srun --account=$LUMI_PROJECT --partition=dev-g -G2 -c 32 --mem=64GB --time=00:30:00 --nodes=1 --pty bash
  source $VENV/bin/activate
  python -c 'import runko'
  ctest --test-dir $BUILD_DIR -j4
EOF
else
    python -c 'import runko'
    cat << EOF

Installed into $VENV; objects in $BUILD_DIR

Multi-rank tests need an allocation:

  srun --account=$LUMI_PROJECT --partition=small --ntasks=4 --cpus-per-task=1 --mem=16G --time=00:30:00 --pty bash
  ctest --test-dir $BUILD_DIR -j4
EOF
fi
