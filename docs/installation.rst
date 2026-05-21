Installation
############

This section contains installation instructions for users of runko.
See the developer notes for a development-friendly installation.

Requirements
============

The software listed below is not obtained automatically and must be made available.
Runko is built with CMake, which is also used to locate dependencies.
Depending on the system configuration, CMake may or may not find the correct libraries automatically.
See the CMake documentation for instructions on how to make CMake find the dependencies.

.. note::

   Runko offically supports installations only to a python virtual environment.
   Create and activate it before proceeding with the installation:

   .. code:: shell

      python -m venv venv
      source venv/bin/activate

   Virtual environment can be deactivated with ``deactivate`` command.

Common requirements (version numbers are the ones tested to work; newer versions should also work):

- Python >=3.11.7
- CMake >=3.23
- MPI

.. tabs::
   .. group-tab:: hip backend

      Consult your Linux distribution's package manager for these packages.

      - LLVM 19-based C++ compiler with HIP support (tested with GCC 14 libstdc++)
      - ROCm 6.3.4
      - MPI implementation that supports GPU-aware MPI

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

      Homebrew's ``open-mpi`` formula builds against AppleClang and fails.
      You need to build OpenMPI manually from the latest stable source with
      Homebrew's LLVM clang instead:

      .. code-block:: bash

         LLVM_PREFIX="$(brew --prefix llvm)"
         OMPI_PREFIX="$HOME/local/openmpi"
         OMPI_VERSION=5.0.10

         mkdir -p $OMPI_PREFIX
         cd $OMPI_PREFIX
         curl -fsSL -O "https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-$OMPI_VERSION.tar.bz2"
         tar -xjf "openmpi-$OMPI_VERSION.tar.bz2"
         cd "openmpi-$OMPI_VERSION"
         ./configure CC="$LLVM_PREFIX/bin/clang" CXX="$LLVM_PREFIX/bin/clang++" \
                     --prefix="$OMPI_PREFIX" \
                     --with-libevent=internal \
                     --with-hwloc=internal
         make -j 8
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

      Obtain the required dependencies by loading modules as shown:

      .. code-block:: shell

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

      .. include:: installation-mpi4py.rst

      .. hint:: On LUMI ``<mpi-aware C compiler>`` is ``cc``.

.. _building:

Building
========

Runko ships a set of CMake configure presets in ``CMakePresets.json`` ---
one per (machine, backend) combination --- so the build is driven by a
single ``--config-settings=cmake.args=--preset=<NAME>`` argument. The
available presets are:
``unix-cpu``, ``macos-cpu``, ``unix-hip``, ``lumi-cpu``, ``lumi-gpu``, ``hile-cpu``, ``hile-gpu``.


.. tabs::
   .. group-tab:: hip backend

      Build and install runko with:

      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko --config-settings=cmake.args=--preset=unix-hip

      .. note::

         ``unix-hip`` will use ``hipcc`` as the compiler by default.


   .. group-tab:: cpu backend

      With ``$ROCTHRUST_INSTALL_PREFIX`` set to the cpu rocThrust install 
      location (see Requirements above) and exported, build and install runko with:


      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko --config-settings=cmake.args=--preset=unix-cpu

      .. hint::

         ``ROCTHRUST_INSTALL_PREFIX`` has to be exported, meaning that it is defined as:

         .. code:: shell

            export ROCTHRUST_INSTALL_PREFIX=<...>

            # or

            ROCTHRUST_INSTALL_PREFIX=<...>
            export ROCTHRUST_INSTALL_PREFIX


   .. group-tab:: cpu backend on macos

      With ``$ROCTHRUST_INSTALL_PREFIX`` set per Requirements above and
      Homebrew's ``mpic++`` on PATH, build and install runko with:

      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko --config-settings=cmake.args=--preset=macos-cpu

   .. group-tab:: hip backend on LUMI

      After loading the modules given above, runko can be built and installed with:

      .. hint::
         Remember to compile runko on a compute node, not on a login node. For example, with an interactive job:

         .. code:: shell

            srun --account=<account_id> --partition=dev-g -G1 -c 32 --mem=64GB --time=00:30:00 --nodes=1 --pty bash


      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko --config-settings=cmake.args=--preset=lumi-gpu


.. hint::

   Overriding preset variables:

   Provide ``--config-settings=cmake.define.<VAR>=<VALUE>`` to CMake.


.. hint::

   Explicitly choosing the compiler:

   Provide ``--config-settings=cmake.define.CMAKE_CXX_COMPILER=<compiler>`` to CMake.

.. hint::

   Compiling a specific git <branch> can be done with ``pip install git+https://github.com/hel-astro-lab/runko@<branch>``


Running
=======

.. tabs::
   .. group-tab:: hip backend

      The HIP backend uses GPU-aware MPI, which has to be enabled:

      .. code:: shell

         export MPICH_GPU_SUPPORT_ENABLED=1


   .. group-tab:: cpu backend

      Nothing special.

   .. group-tab:: hip backend on LUMI

      The HIP backend uses GPU-aware MPI, which has to be enabled:

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
