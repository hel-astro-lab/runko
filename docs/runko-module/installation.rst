Installation
############

Even if runko is a python module it is not (at least not yet) installable through pip or equivalent.
This page explains the steps required for compiling runko with gpu and gpu aware mpi support.


Requirements
============

Version numbers are the ones tested to work.
Newer ones should also work.

- LLVM 17 based C++ compiler that supports HIP (tested with GCC 13 libstdc++)
- Python 3.11.7
- Rocm 6.0.3
- HDF-5
- MPI implementation that supports GPU aware MPI
- CMake 3.23

Runko depends on following Python packages:

- h5py
- matplotlib
- numpy
- scipy
- palettable

Additionally mpi4py package is required.
However, the default installation with `python -m pip install mpi4py`
will most likely link to different mpi implementation that what runko is using.
This can lead to problems so we weed to install mpi4py manually:

.. code:: shell

   MPI4PY_BUILD_MPICC=<mpi-aware C compiler> python -m pip install --no-binary=mpi4py mpi4py


Obtaining Runko
===============

Runko depends on additional C++ libraries but they are vendored in as submodules.

.. code:: shell

   git clone https://github.com/hel-astro-lab/runko
   cd runko
   git checkout tyvifying-clang17
   git submodule update --init --recursive


.. note::

   As of 12.9.2025 `tyvifying-clang17` branch contains latest commits related to runko module.
   This includes refactors that remove use of features not implemented in clang 17
   and many workarounds for compiler bugs on LUMI.


.. tip ::

   To speed up `git submodule update --init --recursive` one can add e.g. `-j 4` flag
   to it, for `4` times parallelized git cloning.


Building
========

Runko is build with CMake. It is also used to locate dependencies.
Depending on the system configuration CMake may or may not find correct libraries automatically.
See CMake documentation on how they should be given.

Assuming dependencies are found automatically and `$CXX` is the compiler we want to use,
then runko can be build with:

.. code:: shell

   cmake -B<build-dir> -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=$CXX .
   make -j<num-of-threads> -C<build-dir>


.. note::

   Runko sets `CMAKE_CXX_FLAGS_{RELEASE,DEBUG}` in `CMakeLists.txt`.
   These can be adjusted for different optimization options.
   This might change in the future.


Building runko will install shared library representing `pyrunko` to `<path-to-runko>/lib`.
Similar pattern is used in corgi, so following paths have to be added to `$PYTHONPATH`
in order to use runko:

.. code:: shell

   RUNKO_PATH="<path-to-runko>"
   P1="$RUNKO_PATH/"
   P2="$RUNKO_PATH/lib"
   P3="$RUNKO_PATH/external/corgi/lib"

   # For example:
   export PYTHONPATH="$PYTHONPATH:$P1:$P2:$P3"


LUMI
----

As of 12.09.2025 following moudles are tested to fulfill the requirements on LUMI:

.. code:: shell

   module load LUMI/24.03
   module load partition/G
   module load PrgEnv-cray
   module load rocm
   module load cray-hdf5
   module load craype-accel-amd-gfx90a
   module load cray-mpich craype-network-ofi
   module load buildtools
   module load cray-python

   # Install Python dependencies...


Now runko can be build with:

.. code:: shell

   cmake -Bbuild -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=CC .
   make -j32 -Cbuild


.. warning::

   Remember not to compile on login nodes. Consider specific slurm job for building
   or using interactive jobs:
   `srun --account=<account> --partition=dev-g --time=00:30:00 --nodes=1 -c32 --pty bash`


To run runko here is a example slurm script (modified from CSC example):

.. code:: bash

   #!/bin/bash -l
   #SBATCH --account=<project-id>
   ##SBATCH --partition=standard-g
   #SBATCH --partition=dev-g
   #
   #SBATCH --job-name=decay
   #SBATCH --output=slurm-%x.out
   #SBATCH --error=slurm-%x.err
   #SBATCH --open-mode=truncate
   #
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=8
   #SBATCH --gpus-per-node=8
   #SBATCH --cpus-per-task=6
   #SBATCH --mem-pre-cpu=10G
   #SBATCH --time=0-00:60:00       # Run time (d-hh:mm:ss)

   # Loads correct modules and sets up PYTHONPATH.
   source <env-script>

   # Required to choose correct GPU for each task.
   cat << EOF > select_gpu
   #!/bin/bash
   export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID

   exec \$*
   EOF

   chmod +x ./select_gpu

   CPU_BIND="mask_cpu:7e000000000000,7e00000000000000"
   CPU_BIND="${CPU_BIND},7e0000,7e000000"
   CPU_BIND="${CPU_BIND},7e,7e00"
   CPU_BIND="${CPU_BIND},7e00000000,7e0000000000"

   export OMP_NUM_THREADS=6
   export MPICH_GPU_SUPPORT_ENABLED=1

   srun --cpu-bind=${CPU_BIND} ./select_gpu python ./runko/projects/pic2-turbulence/pic2-decay.py

   rm -f ./select_gpu


.. warning::

   If you run into a problem where python segfaults or crashes due to illegal instruction
   it can be worked around by adding `import matplotlib.pyplot` before `import runko`.

   I have no idea why this happens nor why importing it first fixes the problem.
