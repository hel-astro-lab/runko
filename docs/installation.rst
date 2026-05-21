Installation
############

This section contains installation instruction for users of runko.
See developer notes for development friendly installation.

Requirements
============

Listed software are not obtained automatically and should be made available.
Runko is build with CMake. It is also used to locate dependencies.
Depending on the system configuration CMake may or may not find correct libraries automatically.
See CMake documentation for instructions on how to make CMake find the dependencies.

.. hint::
   Install runko inside a python virtual environment.

Common requirements
(Version numbers are the ones tested to work.
Newer ones should also work.):

- Python >=3.11.7
- CMake >=3.23
- MPI

.. tabs::
   .. group-tab:: hip backend

      Consult your Linux distributions package manager for these packages.

      - LLVM 19 based C++ compiler with HIP support (tested with GCC 14 libstdc++)
      - Rocm 6.3.4
      - MPI implementation that supports GPU aware MPI

      .. include:: installation-mpi4py.rst

   .. group-tab:: cpu backend

      - C++26 compiler with OpenMP support (e.g. gcc)

      .. include:: installation-mpi4py.rst

      .. include:: installation-cpu-rocthrust.rst

   .. group-tab:: cpu backend on macos

      Install the toolchain and Xcode command-line tools:

      .. code-block:: bash

         brew install llvm cmake wget python3
         xcode-select --install

      Homebrew's ``open-mpi`` formula build against AppleClang and fails.
      You need to build OpenMPI manually from the latest stable source with 
      Homebrew's LLVM clang instead:

      .. code-block:: bash

         LLVM_PREFIX="$(brew --prefix llvm)"
         OMPI_PREFIX="$HOME/local/openmpi"
         OMPI_VERSION=5.0.10

         curl -fsSL -O "https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-$OMPI_VERSION.tar.bz2"
         tar -xjf "openmpi-$OMPI_VERSION.tar.bz2"
         cd "openmpi-$OMPI_VERSION"
         ./configure CC="$LLVM_PREFIX/bin/clang" CXX="$LLVM_PREFIX/bin/clang++" \
                     --prefix="$OMPI_PREFIX" \
                     --with-libevent=internal \
                     --with-hwloc=internal
         make -j "$(sysctl -n hw.ncpu)"
         make install

      Add the install ``bin/`` to ``PATH`` (e.g. append to your ``~/.zshrc``):

      .. code-block:: bash

         export PATH="$HOME/local/openmpi/bin:$PATH"

      Open a fresh shell and verify:

      .. code-block:: bash

         which mpic++         # -> $HOME/local/openmpi/bin/mpic++
         mpic++ --version     # -> Homebrew clang version <N>

      .. include:: installation-mpi4py.rst

      .. include:: installation-cpu-rocthrust.rst

   .. group-tab:: hip backend on LUMI

      Obtain required dependencies by loading modules as show:

      .. code-block:: shell

         module load LUMI/25.09 \
             partition/G \
             PrgEnv-cray \
             rocm/6.4.4 \
             craype-accel-amd-gfx90a \
             cray-mpich/9.0.1 \
             craype-network-ofi \
             buildtools \
             cray-python \
             lumi-CrayPath

      .. include:: installation-mpi4py.rst

      .. hint:: On LUMI ``<mpi-aware C compiler>`` is ``cc``.

.. _building:

Building
========

Runko ships a set of CMake configure presets in ``CMakePresets.json`` ---
one per (machine, backend) combination --- so the build is driven by a
single ``--config-settings=cmake.args=--preset=<NAME>`` argument. Available
presets:
``unix-cpu``, ``macos-cpu``, ``unix-gpu``, ``lumi-cpu``, ``lumi-gpu``, ``hile-cpu``, ``hile-gpu``.


.. tabs::
   .. group-tab:: hip backend

      Build and install git branch ``<branch>`` of runko with:

      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko@<branch> \
           --config-settings=cmake.args=--preset=unix-gpu

      .. note::

         ``unix-gpu`` will use ``hipcc`` as the compiler by default.


   .. group-tab:: cpu backend

      With ``$ROCTHRUST_INSTALL_PREFIX`` set to the cpu rocThrust install
      location (see Requirements above), build and install git branch
      ``<branch>`` of runko with:

      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko@<branch> \
           --config-settings=cmake.args=--preset=unix-cpu


   .. group-tab:: cpu backend on macos

      With ``$ROCTHRUST_INSTALL_PREFIX`` set per Requirements above and
      Homebrew's ``mpic++`` on PATH, build and install git branch
      ``<branch>`` of runko with:

      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko@<branch> \
           --config-settings=cmake.args=--preset=macos-cpu

   .. group-tab:: hip backend on LUMI

      After loading the modules given above,
      then git branch ``<branch>`` of runko can be build and installed with:

      .. hint::
         Remember to compile runko on a compute node
         and not to compile runko on a login node.
         For example with a interactive job:

         .. code:: shell

            srun --account=<account_id> --partition=dev-g -G1 -c 32 --mem=64GB --time=00:30:00 --nodes=1 --pty bash


      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko@<branch> \
           --config-settings=cmake.args=--preset=lumi-gpu


.. hint::

   Overriding preset variables:

   Provide ``--config-settings=cmake.define.<VAR>=<VALUE>`` to CMake.


.. hint::

   Explicitly choosing the compailer:

   Provide ``--config-settings=cmake.define.CMAKE_CXX_COMPILER=<compiler>`` to CMake.

Running
=======

.. tabs::
   .. group-tab:: hip backend

      Hip backend uses GPU-aware MPI which has to be enabled:

      .. code:: shell

         export MPICH_GPU_SUPPORT_ENABLED=1


   .. group-tab:: cpu backend

      Nothing special.

   .. group-tab:: hip backend on LUMI

      Hip backend uses GPU-aware MPI which has to be enabled:

      .. code:: shell

         export MPICH_GPU_SUPPORT_ENABLED=1


      Example slurm job script:

      .. code:: shell

         #!/bin/bash -l
         # source: https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
         #SBATCH --job-name=examplejob   # Job name
         #SBATCH --output=examplejob.o%j # Name of stdout output file
         #SBATCH --error=examplejob.e%j  # Name of stderr error file
         #SBATCH --partition=standard-g  # partition name
         ## SBATCH --exclusive           # Uncomment if not on standard-g partition.
         #SBATCH --nodes=2               # Total number of nodes
         #SBATCH --ntasks-per-node=8     # 8 MPI ranks per node, 16 total (2x8)
         #SBATCH --gpus-per-node=8       # Allocate one gpu per MPI rank
         #SBATCH --time=1-12:00:00       # Run time (d-hh:mm:ss)
         #SBATCH --account=project_<id>  # Project for billing

         # copy-paste module loads from requirements here:
         #
         #     module load ...
         #
         # or source a file that contains them:
         #
         #    source runko-modules.sh

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

         srun --cpu-bind=${CPU_BIND} ./select_gpu python <file>
         rm -rf ./select_gpu

      .. warning::

         - If runko crashes due to invalid GPU memory access,
           first try ``export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000``
           and if that does not fix it then try
           ``export MPICH_GPU_IPC_ENABLED=0``.
         - If you run into segfaults or crashes due to illegal instruction,
           it can be worked around by adding
           ``import matplotlib.pyplot`` before ``import runko``.

Now we can test if the installed package can be imported:

.. code:: shell

   python -c 'import runko'
