Installation
############

This section describes the basic installation of the framework. See also the :doc:`clusters` page in case you are installing the code to some computing cluster.


Git + GitHub access with SSH
============================

Before we start, we need a developer level access to GitHub. Most importantly, you need to ensure that you can access git repositories via the SSH connection. Follow these tutorials to link your SSH key to your GitHub account:

* https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
* https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account


Downloading/Cloning
===================

The framework relies on various small (template) libraries that are automatically downloaded along the main code as git submodules. Because of this, remember to always issue a recursive clone when downloading the code from GitHub:

.. code-block:: bash

   git clone --recursive git@github.com:natj/runko.git

You also need to update all the submodules (so that various branches are in sync) with

.. code-block:: bash

   git submodule update --recursive

External libraries
==================

External libraries installed together with the framework include:

* `corgi <https://github.com/natj/corgi>`_ massively-parallel grid infrastructure
* `PyBind11 <https://github.com/pybind/pybind11>`_ binding library for seamless operability between C++ and Python
* `cppitertools <https://github.com/ryanhaining/cppitertools>`_ for various python-like C++ iterators etc.
* `ezh5 <https://github.com/natj/ezh5>`_ lazy template interface to HDF5 for easier usage (originally from `mileschen <https://github.com/mileschen360/ezh5>`_ )
* `mpi4cpp <https://github.com/natj/mpi4cpp>`_ high-level interface to MPI library


Requirements
============

Before proceeding with the code compilation, check that your system fulfills these requirements:

* Modern C++ compiler (such g++, `Clang++ <https://clang.llvm.org/>`_, etc.)
* python3 (with mpi4py, numpy, scipy, matplotlib)
* `HDF5 <https://support.hdfgroup.org/HDF5/>`_ for I/O (serial mode is enough)
* `CMake <https://cmake.org/>`_ (>v3.0) for building and compiling
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


MPI needs to be compiled separately because, by default, it uses the AppleClang compiler (instead of the `g++` that you just installed).

You can compile OpenMPI via homebrew by first modifying your `~/.bash_profile` to link the new compilers:

.. code-block:: bash

   export HOMEBREW_CC=gcc-13
   export HOMEBREW_CXX=g++-13
   export OMPI_CC=gcc-13
   export OMPI_CXX=g++-13

Then restart the terminal to reload the newly added environment variables. After restarting, install `openmpi` from source with

.. code-block:: bash

    brew reinstall openmpi -cc=gcc-13 --build-from-source

and 

.. code-block:: bash

    brew reinstall mpi4py



Linux (Ubuntu)
--------------

When compiling runko and running the scripts, it is critical that you always use the same Python interpreter, C/C++ compiler, and associated OpenMPI distribution, otherwise this can give several errors during the installation. For this reason we recommend using vanilla `python` and disabling anaconda (if you are using it) by commenting out its activation in your ``~/.bashrc`` file.

.. code-block:: bash

   # >>> conda initialize >>>
   # ...
   # <<< conda initialize <<<

You may find it also necessary to delete folders containing the older Python versions than your current one at `/usr/bin/python3.*`. In order to get a completely clean OpenMPI distribution first run:

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


Manual installation of OpenMPI (optional)
-----------------------------------------

Alternatively, if you want even more control of the operation, you can compile it manually yourself by running:

.. code-block:: bash

   export MPI_IMPL=openmpi41
   mkdir -p $HOME/local/$MPI_IMPL/bin/openmpi
   cd $HOME/local/$MPI_IMPL/bin/openmpi
   wget --no-check-certificate http://www.open-mpi.org/software/ompi/v4.1/downloads/openmpi-4.1.5.tar.bz2
   tar -xjf openmpi-4.1.5.tar.bz2
   cd openmpi-4.1.5
   export OMPI_CC=gcc-13
   export OMPI_CXX=g++-13
   ./configure CC=gcc-13 CXX=g++-13 --prefix=$HOME/bin/$MPI_IMPL 
   make -j 4
   make install
   make clean


This installs OpenMPI to `~/bin/` and exports the correct directories so that the `mpic++` compiler wrapper becomes available. You should then add to your `.bash_profile` (or `.zshrc` in latest MacOS) these exports (in case you need to re-compile the library):

.. code-block:: bash

   export OMPI_CC=gcc-13
   export OMPI_CXX=g++-13
   export MPI_IMPL=openmpi41
   export PATH=$PATH:$HOME/bin/$MPI_IMPL/bin
   export PATH=$PATH:$HOME/bin/$MPI_IMPL/include
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/bin/$MPI_IMPL/lib


After `openmpi` is installed we also need to re-install `mpi4py` because it uses the system-default mpi installation

.. code-block:: bash

   pip3 uninstall mpi4py --break-system-packages
   pip3 install mpi4py --break-system-packages


Note the additional ``--break-system-packages`` keyword that is needed for the latest python versions ``>3.12`` to install packages with pip and homebrew/apt-get.


Python libraries
================

All the python requirements can be installed via `pip` as

.. code-block:: bash

   pip install -r requirements.txt

.. note::

    If you had to manually install MPI in the previous section, then remember to re-install mpi4py.



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

You should also put this part into your `~/.bashrc` (or `~/.zshrc` in the latest MacOS) so that correct compilers are automatically exported during the startup.

You should also add the python script directories into `PYTHONPATH` environment variable. Modify your `~/.zshrc` (MacOS) or `~/.bashrc` (Linux) by appending `corgi` and `runko` libraries to the path by exporting

.. code-block:: bash

    export RUNKO=/path2repo
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/"
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/lib"
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/corgi/lib"
    PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$RUNKO/bindings/old"
    export PYTHONPATH

where `path2repo` points to the location where you cloned the repository (i.e. path to `runko` directory). Note that there is no trailing slash `/` in the commands. As an example, the path can be e.g., `/Users/natj/runko`.


Next we can proceed to compiling. Out-of-source builds are recommended: inside the repository directory, make a new `build` directory, go into that, and only then run the CMake configuration commands. This can be done by running (inside `runko` directory):

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



Testing of the new installation
-------------------------------

Runko comes with multiple tests (found in ``runko/tests/``) that are ran after every compilation. In general, if you see "All tests passes.", after the compilation, your installation should be succesfull.

The next step is to run an actual simulation with the code. For that, see the premade projects setups in ``runko/projects/``.
