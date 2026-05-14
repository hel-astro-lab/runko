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

      On MacOS these are easily installed with `homebrew <https://brew.sh/>`_ by running:

      .. code-block:: bash

         brew install gcc python3 cmake wget


      Make also sure that Xcode developer tools are installed by running

      .. code-block:: bash

         xcode-select --install


      MPI needs to be compiled separately because,
      by default, it uses the `AppleClang` compiler (instead of the `g++` that you just installed).

      You can compile OpenMPI via homebrew by first modifying your ``~/.zshrc`` to link the new compilers:

      .. code-block:: bash

         export HOMEBREW_CC=gcc-15
         export HOMEBREW_CXX=g++-15
         export OMPI_CC=gcc-15
         export OMPI_CXX=g++-15
         export CC=gcc-15  # NOTE: you need to change these after the installation
         export CXX=g++-15 # NOTE: you need to change these after the installation

      Then restart the terminal to reload the newly added environment variables.
      After restarting, install `OpenMPI` from source with

      .. code-block:: bash

          brew reinstall openmpi --build-from-source

      and then remove all the (possible) previous installations of mpi4py and re-install it using pip3

      .. code-block:: bash

          brew uninstall mpi4py
          pip3 uninstall mpi4py --break-system-packages
          pip3 install mpi4py --break-system-packages

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

.. tabs::
   .. group-tab:: hip backend

      Assuming dependencies are found automatically
      and ``<hip-compiler>`` is a compiler with HIP support,
      then git branch ``<branch>`` of runko can be build and installed with:

      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko@<branch> \
         --config-settings=cmake.define.tyvi_BACKEND=hip \
         --config-settings=cmake.define.CMAKE_CXX_COMPILER=<hip-compiler>

   .. group-tab:: cpu backend

      Assuming dependencies are found automatically
      and ``$ROCTHRUST_INSTALL_PREFIX`` is the cpu rocthrust install location
      (see Requirements from above),
      then git branch ``<branch>`` of runko can be build and installed with:

      .. code:: shell

         pip install git+https://github.com/hel-astro-lab/runko@<branch> \
         --config-settings=cmake.define.tyvi_BACKEND=cpu \
         --config-settings=cmake.define.rocthrust_DIR=$ROCTHRUST_INSTALL_PREFIX/lib/cmake/rocthrust/

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
         --config-settings=cmake.define.tyvi_BACKEND=hip \
         --config-settings=cmake.define.CMAKE_CXX_COMPILER=CC


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
