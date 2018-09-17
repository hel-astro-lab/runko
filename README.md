# Kinetic plasma simulation toolbox

Plasma-toolbox is a collection of simulation modules written in modern C++14/17 to model astrophysical plasmas. The framework is build on top of massively parallel heterogeneous template library [corgi](https://github.com/natj/corgi) and heavily relies on presenting the code as physical modules that can be combined or run individually.


## Installation

Plasmabox relies on various small template libraries that are automatically obtained along the main code as submodules. Because of this, one should always when downloading the repository issue a recursive clone as
```
git clone --recursive https://github.com/natj/plasma-toolbox.git
```
Currently these include:
- [cppitertools](https://github.com/ryanhaining/cppitertools) for various python-like iterators etc.
- [ezh5](https://github.com/natj/ezh5) lazy template interface to HDF5 for easier usage (originally from [mileschen](https://github.com/mileschen360/ezh5))
- [inih](https://github.com/benhoyt/inih) `.ini` file parser for configuration files.


### Requirements
Your system should have these installed:
- Modern C++ compiler (such as [Clang++](https://clang.llvm.org/), g++ [>v5.0], or similar)
- python3 for scripting
- [HDF5](https://support.hdfgroup.org/HDF5/) for I/O
- [CMake](https://cmake.org/) (>v3.0) for building and compiling
- [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra
- [FFTW](http://www.fftw.org/) for Fourier transforms
- [fmt](https://github.com/fmtlib/fmt) I/O formatting library

On MacOS these should be (quite) easily obtained by using [homebrew](https://brew.sh/) and running.
- `brew install gcc hdf5 python3 cmake eigen3 fftw fmt`

Additionally, especially for larger jobs and better performance, consider installing:
- MPI message passing library
- openMP (>v4.0) supported by the compiler 
    - NOTE: default clang++ in MacOS does not support this


### Compiling

First, you need to configure the build. Out-of-source builds are recommended so inside the repository run:
```
mkdir build
cd build
cmake ..
```
And make sure cmake finishes succesfully. After that, you can try and compile the framework as
```
make -j8
```

Different compilers can be selected by modifying the `CXX` and `CC` global environment variables.


After compilation, you should add the python modules into `PYTHONPATH` environment variable for easier access. Modify your `~/.bash_profile` (MacOS) or `~/.bashrc` (linux) by appending `corgi` and `plasmabox` libraries to the path as
```
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/corgi/lib
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/lib
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/python
export PYTHONPATH=$PYTHONPATH:/path2repo/plasmabox/analysis
```



## Code directory structure / modules
- `analysis`: analysis scripts for simulations
- `docs`: Doxygen documentation
- `jobs`: various cluster spesific job files and maintainance scripts
- `notes`: random (latex and markdown) notes 
- `projects`: projects files and tests
    - new simulation setups should go here
- `prototypes`: different prototypes and test versions of the code
- `python`: python3 bindings and scripts for running simulations using python
- `tests`: unittests
- `tools`: C++ libraries, bash scripts, git tools

### Physical modules include:
- `corgi`: corgi library
- `em-fields`: Maxwell field solver module
- `io`: input/output module
- `pic`: Particle-in-Cell module
- `vlasov`: Relativistic Vlasov module



