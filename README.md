# Modern kinetic plasma simulation toolbox
<img align="top" src="notes/header.png">

[![Build Status](https://travis-ci.com/natj/plasmabox.svg?branch=master)](https://travis-ci.com/natj/plasmabox)

PlasmaBox is a collection of simulation modules written in modern C++14 to model astrophysical plasmas. The framework is build on top of massively parallel heterogeneous template library [corgi](https://github.com/natj/corgi) and relies on presenting the code as physical modules that can be either run individually or combined in a heterogeneous fashion.


## Installation

PlasmaBox relies on various small template libraries that are automatically obtained along the main code as submodules. Because of this, remember to always issue a recursive clone:
```
git clone --recursive https://github.com/natj/plasmabox.git
```
Currently the libraries shipped together with the code are:
- [corgi](https://github.com/natj/corgi) grid infrastructure
- [cppitertools](https://github.com/ryanhaining/cppitertools) for various python-like iterators etc.
- [ezh5](https://github.com/natj/ezh5) lazy template interface to HDF5 for easier usage (originally from [mileschen](https://github.com/mileschen360/ezh5))


### Requirements
Your system should have these installed:
- Modern C++ compiler (such as [Clang++](https://clang.llvm.org/), g++ [>v7.0], or similar)
- python3 for scripting
- [HDF5](https://support.hdfgroup.org/HDF5/) for I/O
- [CMake](https://cmake.org/) (>v3.0) for building and compiling
- [fmt](https://github.com/fmtlib/fmt) I/O formatting library
- [FFTW](http://www.fftw.org/) for Fourier transforms
- MPI message passing interface


### MacOS
On MacOS these should be (quite) easily obtained by using [homebrew](https://brew.sh/) and running:
- `brew install gcc@7 hdf5 python3 open-mpi cmake fftw fmt eigen`

### Linux
On Linux (assuming Ubuntu) run:
- `sudo -E apt-add-repository -y "ppa:ubuntu-toolchain-r/test"`
- `sudo apt-get install install clang-5.0 g++-7 hdf5-tools python3 libopenmpi-dev`


All the python requirements can be installed via pip as
- `pip install -r requirements.txt`
Note that on MacOS you should use `pip3` instead to get the correct python3 version of the libraries.



### Compiling

First you need to configure the build. To use your (freshly installed) modern C++ compiler export it as
```
export CC=gcc-7
export CXX=g++-7
```
You can also put this part into your `~/.bashrc` (or `~/.bash_profile` on MacOS) so correct compilers are automatically exported.

Next we can proceed to compiling. Out-of-source builds are recommended so inside the repository run:
```
mkdir build
cd build
cmake ..
```
And make sure that CMake finishes successfully. After that, you can try and compile the complete framework with
```
make
```

When compiling and linking is finished, the CMake runs few automated tests to check that everything is working. You should see a message *"XX tests finished succesfully"* in the end, if the build was successful.


After compilation, you should add the python modules into `PYTHONPATH` environment variable for easier access. Modify your `~/.bash_profile` (MacOS) or `~/.bashrc` (Linux) by appending `corgi` and `plasmabox` libraries to the path as
```
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/corgi/lib
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/lib
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/python
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/analysis
```
where `path2repo` points to the location where you cloned the repository.



## Code directory structure / modules
- `analysis`: analysis scripts for simulations
- `docs`: Doxygen documentation
- `jobs`: various cluster spesific job files and maintainance scripts
- `notes`: random (latex and markdown) notes 
- `projects`: projects files and tests
    - new simulation setups should go here
- `prototypes`: different prototypes and test versions of the code
- `python`: python3 bindings and scripts for running simulations using python
- `tests`: code unit tests
- `tools`: C++ libraries, bash scripts, git tools

### Physical modules include:
- `corgi`: corgi library
- `em-fields`: Maxwell field solver module
- `io`: input/output module
- `pic`: Particle-in-Cell module
- `vlasov`: Relativistic Vlasov module
- `radiation`: Radiation module 



