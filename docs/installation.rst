Installation
############

This section describes the basic installation of the framework. See also the :doc:`clusters` page in case you are installing the code to some computing cluster.


Downloading/Cloning
===================

The framework relies on various small (template) libraries that are automatically obtained along the main code as submodules. Because of this, remember to always issue a recursive clone:

.. code-block:: bash

   git clone --recursive https://github.com/natj/runko.git

You also need to update all the submodules (so that other branches are also in sync) with

.. code-block:: bash

   git submodule update --recursive

External libraries
==================

External libraries installed together with the framework include:

* `corgi <https://github.com/natj/corgi>`_ massively-parallel grid infrastructure
* `PyBind11 <https://github.com/pybind/pybind11>`_ binding library for seamless operability between C++ and Python
* `cppitertools <https://github.com/ryanhaining/cppitertools>`_ for various python-like C++ iterators etc.
* `ezh5 <https://github.com/natj/ezh5>`_ lazy template interface to HDF5 for easier usage (originally from [mileschen](https://github.com/mileschen360/ezh5))
* `mpi4cpp <https://github.com/natj/mpi4cpp>`_ high-level interface to MPI library


Requirements
============

Before proceeding to compilation, check that your system has these requirements installed:

* Modern C++ compiler (such a`Clang++ <https://clang.llvm.org/>`_, g++, ...)
* python3 (with mpi4py)
* `HDF5 <https://support.hdfgroup.org/HDF5/>`_ for I/O (serial mode is enough)
* `CMake <https://cmake.org/>`_ (>v3.0) for building and compiling
* `FFTW <http://www.fftw.org/>`_ for Fourier transforms
* MPI message passing interface


.. note::

    Note that g++-8 does not work because of a (known) compiler bug. Therefore, g++-9 or newer is the current recommended choice.



MacOS
-----

On MacOS these are easily installed with `homebrew <https://brew.sh/>`_ by running:

.. code-block:: bash

   brew install gcc hdf5 python3 cmake wget


Make also sure that Xcode developer tools are installed by running

.. code-block:: bash

   xcode-select --install


MPI needs to be compiled separately because, by default, it uses the AppleClang compiler (instead of the latest `g++` just installed; currently `g++-12`). 

You can compile OpenMPI via homebrew by first modifying your `~/.bash_profile` to link the new compilers:

.. code-block:: bash

   export HOMEBREW_CC=gcc-12
   export HOMEBREW_CXX=g++-12

Then restart the terminal to reload the newly added environment variables. After restarting, install `openmpi` from source with

.. code-block:: bash

    brew reinstall openmpi --build-from-source



Alternatively, if you want even more control of the operation, you can compile it manually yourself by running:

.. code-block:: bash

   export MPI_IMPL=openmpi41
   mkdir -p $HOME/local/$MPI_IMPL/bin
   cd $HOME/local/$MPI_IMPL/bin
   mkdir -p openmpi && cd openmpi
   wget --no-check-certificate http://www.open-mpi.org/software/ompi/v4.1/downloads/openmpi-4.1.4.tar.bz2
   tar -xjf openmpi-4.1.4.tar.bz2
   cd openmpi-4.1.4
   export OMPI_CC=gcc-12
   export OMPI_CXX=g++-12
   ./configure CC=gcc-12 CXX=g++-12 --prefix=$HOME/local/$MPI_IMPL 
   make -j 4
   make install
   make clean
   cd ../../

   export PATH=$PATH:$HOME/local/$MPI_IMPL/bin
   export PATH=$PATH:$HOME/local/$MPI_IMPL/include
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/$MPI_IMPL/lib

This installs OpenMPI to `~/local` and exports the correct directories so that `mpic++` compiler wrapper becomes available. You should put the last 3 export commands to your `.bash_profile` for easier usage, in case you need to recompile Runko at some point.

Linux (Ubuntu)
-----

When compiling runko and running runko scripts, it is critical that you always use the same Python interpreter, C/C++ compiler and associated OpenMPI distribution, otherwise this can give several errors during installation. For this reason we recommend using vanilla Python and disabling anaconda (if you are using it) by commenting out its activation in your ``~/.bashrc`` file.

.. code-block:: bash

   # >>> conda initialize >>>
   # ...
   # <<< conda initialize <<<

You may find it necessary to delete folders containing older Python versionsthan your current one at `/usr/bin/python3.*`. In order to get a completely clean OpenMPI distribution first run:

.. code-block:: bash

   sudo apt-get remove mpich libopenmpi-dev openmpi-bin
   sudo apt-get update && sudo apt-get autoclean && sudo apt-get clean && sudo apt-get autoremove

Then run:

.. code-block:: bash

   sudo -E apt-add-repository -y "ppa:ubuntu-toolchain-r/test"
   sudo apt-get install libopenmpi-dev libhdf5-serial-dev hdf5-helpers openmpi-bin libblas-dev liblapack-dev python3 python3-pip

.. note::

   Recent Ubuntu (bionic) comes with gcc-7 which makes the installation easier. For previous versions you, additionally, need to install gcc-7 (or 9) and manually compile MPI similar to the MacOS discussed above.

You also need to export the HDF5 library location (since it is non-standard at least in Ubuntu) with

.. code-block:: bash

   export HDF5_INCLUDE_PATH=/usr/include/hdf5/serial

Finally, you can test that your Runko installation in Python is working properly by executing the test described in ``runko/projects/pic-shocks/README.md``.

Python libraries
================

All the python requirements can be installed via `pip` as

.. code-block:: bash

   pip3 install -r requirements.txt

.. note::

    If you had to manually install MPI in the previous section, then you need to remove mpi4py (`pip3 uninstall mpi4py`) and re-install it.



Compiling
=========

After installing all the pre-requisites, you can proceed to compiling. First you need to configure the build. To use your (freshly installed) modern C++ compiler we need to export them as

.. code-block:: bash

   export CC=mpicc
   export CXX=mpic++

Then make sure that everything works, check the output of

.. code-block:: bash

   $CC --version
   $CXX --version

This should indicate that the newly installed compilers are used.


You should also put this part into your `~/.bashrc` (or `~/.bash_profile` on MacOS) so correct compilers are automatically exported in the startup.

You should also add the python script directories into `PYTHONPATH` environment variable. Modify your `~/.bash_profile` (MacOS) or `~/.bashrc` (Linux) by appending `corgi` and `runko` libraries to the path by exporting

.. code-block:: bash

    export RUNKO=/path2repo
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/"
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/lib"
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/corgi/lib"
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/bindings/old"
    export PYTHONPATH

where `path2repo` points to the location where you cloned the repository (i.e. path to `runko` directory). Note that there is no trailing slash `/`. As an example, the path can be e.g., `/Users/natj/runko`.


Next we can proceed to compiling. Out-of-source builds are recommended so inside the repository make a new build directory, go into that and only then run the CMake. This can be done by running:

.. code-block:: bash

   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=$(which python3) ..

And make sure to check that `CMake` finishes successfully. After that, you are ready to compile the framework with

.. code-block:: bash

   make

When compiling and linking is finished, CMake runs few automated tests to check that everything is working. You should see a message *"XX tests finished succesfully"* in the end, if the build was successful.


.. note::

    Since the compiling can take quite a while, you can use the multi-core compilation by passing make the `-j8` option (or whatever number of tasks you want).


