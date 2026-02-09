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
   git checkout dev-v5
   git submodule update --init --recursive


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


Building runko will install shared library representing `pyrunko` to repository root.
Similar pattern is used in corgi, so following paths have to be added to `$PYTHONPATH`
in order to use runko:

.. code:: shell

   RUNKO_PATH="<path-to-runko>"
   P1="$RUNKO_PATH/"
   P2="$RUNKO_PATH/external/corgi/lib"

   # For example:
   export PYTHONPATH="$PYTHONPATH:$P1:$P2"


LUMI
----

The easiest method to install runko on LUMI
is by using the handy script `archs/lumi-install.sh </archs/lumi-install.sh>`_ as follows:

.. code:: shell

   git clone git@github.com:hel-astro-lab/runko.git
   cd runko
   git checkout dev-v5
   git submodule update --init --recursive
   ./archs/lumi-install.sh


This will build runko in one go and append all necessary LUMI module loads to the virtual environment loader which can be used as follows:

.. code:: shell

   source runko-venv/bin/activate


to enable runko and all its required modules on a LUMI node.


**Warning: overloading login nodes**

Remember not to compile on login nodes. Consider using a specific slurm job instead for building or interactive jobs:


.. code:: shell

   srun --account=<account> --partition=dev-g --time=00:30:00 --nodes=1 -c32 --pty bash


**Running on LUMI**

For example run script on LUMI see: `projects/pic-turbulence/jobs/multi-node.lumi?`.


.. tip ::

   If runko crashes due to invalid GPU memory access,
   then adding `export MPICH_GPU_IPC_ENABLED=0` to the job script might fix the problem.


**Warning: LUMI segfaults**

If you run into a problem where python segfaults or crashes due to illegal instruction
it can be worked around by adding `import matplotlib.pyplot` before `import runko`.

I have no idea why this happens nor why importing it first fixes the problem.
