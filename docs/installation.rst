Installation
############

This section describes the basic installation of the framework.


Downloading/Cloning
===================

The framework relies on various small template libraries that are automatically obtained along the main code as submodules. Because of this, remember to always issue a recursive clone:

.. code-block:: bash

   git clone --recursive https://github.com/natj/runko.git

You also need to update all the submodules (so that other branches are also in sync) with

.. code-block:: bash

   git submodule update --recursive

External libraries
==================

Current external libraries shipped together with the framework are:

* `corgi <https://github.com/natj/corgi>`_ massively-parallel grid infrastructure
* `PyBind11 <https://github.com/pybind/pybind11>`_ binding library for seamless operability between C++ and Python
* `cppitertools <https://github.com/ryanhaining/cppitertools>`_ for various python-like C++ iterators etc.
* `ezh5 <https://github.com/natj/ezh5>`_ lazy template interface to HDF5 for easier usage (originally from [mileschen](https://github.com/mileschen360/ezh5))
* `mpi4cpp <https://github.com/natj/mpi4cpp>`_ high-level interface to MPI library


Requirements
============

Before proceeding to compilation, check that your system has these requirements installed:

* Modern C++ compiler (such as `Clang++ <https://clang.llvm.org/>`_, g++, ...)
* python3 
* `HDF5 <https://support.hdfgroup.org/HDF5/>`_ for I/O (serial mode is enough)
* `CMake <https://cmake.org/>`_ (>v3.0) for building and compiling
* `FFTW <http://www.fftw.org/>`_ for Fourier transforms
* MPI message passing interface


.. note::

    Note that g++-8 does not work because of a (known) compiler bug. Therefore, g++-9 is the recommended choice.


MacOS
-----

On MacOS these should be easily installed by using `homebrew <https://brew.sh/>`_ and running:

.. code-block:: bash

   brew install gcc@9 hdf5 python3 open-mpi cmake fftw


Linux
-----

On Linux (assuming Ubuntu) run:

.. code-block:: bash

   sudo -E apt-add-repository -y "ppa:ubuntu-toolchain-r/test"
   sudo apt-get install install clang-5.0 g++-9 hdf5-tools python3 python3-pip openmpi-bin libopenmpi-dev

Python libraries
================

All the python requirements can be installed via pip as

.. code-block:: bash

   pip3 install -r requirements.txt



Compiling
=========

After installing all the pre-requisites, you can proceed to compiling. First you need to configure the build. To use your (freshly installed) modern C++ compiler we need to export it as

.. code-block:: bash

   export CC=gcc-9
   export CXX=g++-9

You can also put this part into your `~/.bashrc` (or `~/.bash_profile` on MacOS) so correct compilers are automatically exported in the startup.

You should also add the python script directories into `PYTHONPATH` environment variable. Modify your `~/.bash_profile` (MacOS) or `~/.bashrc` (Linux) by appending `corgi` and `runko` libraries to the path by exporting

.. code-block:: bash

   export RUNKO=/path2repo/
   export PYTHONPATH=$PYTHONPATH:$RUNKO/corgi/lib
   export PYTHONPATH=$PYTHONPATH:$RUNKO/lib
   export PYTHONPATH=$PYTHONPATH:$RUNKO/python
   export PYTHONPATH=$PYTHONPATH:$RUNKO/analysis


where `path2repo` points to the location where you cloned the repository (i.e. path to `runko` directory).


Next we can proceed to compiling. Out-of-source builds are recommended so inside the repository run:

.. code-block:: bash

   mkdir build
   cd build
   cmake ..

And make sure that CMake finishes successfully. After that, you can try and compile the complete framework with

.. code-block:: bash

   make

When compiling and linking is finished, CMake runs few automated tests to check that everything is working. You should see a message *"XX tests finished succesfully"* in the end, if the build was successful.


.. note::

    Since the compiling can take quite a while, you can use the multi-core compilation by passing make the `-j8` option (or whatever number of tasks you want).


