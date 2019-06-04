# Modern kinetic plasma simulation toolbox
<img align="top" src="notes/header.png">

[![Build Status](https://travis-ci.com/natj/plasmabox.svg?branch=master)](https://travis-ci.com/natj/plasmabox) [![MIT](https://badges.frapsoft.com/os/mit/mit.svg?v=102)](https://github.com/natj/plasmabox/LICENSE)


PlasmaBox is a collection of simulation codes written in modern C++14/Python3 to model astrophysical plasmas. The framework consists of various physical modules that can be run independently or combined to create multi-physics simulations. Low-level "kernels" are mainly implemented in modern C++ that allows to write modular and high performance code. By binding these fast low-level classes to Python objects it is also easy to use and extend them. This ensures efficient code, rapid prototyping, and ease of use.

Under the hood, the framework uses the massively parallel grid infrastructure library [corgi](https://github.com/natj/corgi) that relies on decomposing the grid to smaller subregions, called tiles, that can be operated on and updated independently of each other. [Corgi](https://github.com/natj/corgi) also automatically parallelizes the simulations and provides dynamic load-balancing capability. Therefore, small simulation setups can be tested locally on laptops and then extended for massively parallel supercomputer platforms (currently tested up to ~10k cores).


Current main physical simulation modules include:
- **Electromagnetic field solver** based on staggered Yee lattice (`fields/`)
- **2D3V Particle-In-Cell module** (`pic/`)
- **1D3V Relativistic Vlasov module** (`vlasov/`)

Additionally, modules under construction include:
- **Force-free MHD** module (`ffe/`)
- Non-linear Monte Carlo **radiation module** (`radiation/`)
- Full **3D3V Vlasov module**


## Quick getting started guide
1) Follow the installation instructions below to get PlasmaBox running on your laptop.
2) Add your own project repository under `projects/` directory.
	- This can be based on, for example, Python drivers in `projects/tests/`
3) Test & prototype your new simulation setups on your laptop/desktop.
4) Simulate on supercomputers!


## Installation

PlasmaBox relies on various small template libraries that are automatically obtained along the main code as submodules. Because of this, remember to always issue a recursive clone:
```
git clone --recursive https://github.com/natj/plasmabox.git
```
You also need to update all the submodules (so that other branches are also in sync) with
```
git submodule update --recursive
```
Current external libraries shipped together with the framework are:
- [corgi](https://github.com/natj/corgi) massively-parallel grid infrastructure
- [PyBind11](https://github.com/https://github.com/pybind/pybind11) binding library for seamless operability between C++ and Python
- [cppitertools](https://github.com/ryanhaining/cppitertools) for various python-like C++ iterators etc.
- [ezh5](https://github.com/natj/ezh5) lazy template interface to HDF5 for easier usage (originally from [mileschen](https://github.com/mileschen360/ezh5))
- [mpi4cpp](https://github.com/natj/mpi4cpp) high-level interface to MPI library


### Requirements
Before proceeding to compilation, check that your system has these requirements installed:
- Modern C++ compiler (such as [Clang++](https://clang.llvm.org/), g++ [>v6.0], or similar. Note that gcc8 does not work atm because of a compiler bug)
- python3 for scripting
- [HDF5](https://support.hdfgroup.org/HDF5/) for I/O (serial mode is enough)
- [CMake](https://cmake.org/) (>v3.0) for building and compiling
- [FFTW](http://www.fftw.org/) for Fourier transforms
- MPI message passing interface


### MacOS
On MacOS these should be easily installed by using [homebrew](https://brew.sh/) and running:
- `brew install gcc@7 hdf5 python3 open-mpi cmake fftw`

### Linux
On Linux (assuming Ubuntu) run:
- `sudo -E apt-add-repository -y "ppa:ubuntu-toolchain-r/test"`
- `sudo apt-get install install clang-5.0 g++-7 hdf5-tools python3 python3-pip openmpi-bin libopenmpi-dev`

### Python libraries
All the python requirements can be installed via pip as
- `pip3 install -r requirements.txt`



### Compiling
After installign all the pre-requisites, you can proceed to compiling. First you need to configure the build. To use your (freshly installed) modern C++ compiler we need to export it as
```
export CC=gcc-7
export CXX=g++-7
```
You can also put this part into your `~/.bashrc` (or `~/.bash_profile` on MacOS) so correct compilers are automatically exported in the startup.

You should also add the python script directories into `PYTHONPATH` environment variable. Modify your `~/.bash_profile` (MacOS) or `~/.bashrc` (Linux) by appending `corgi` and `plasmabox` libraries to the path by exporting
```
export PLASMABOXDIR=/path2repo/
export PYTHONPATH=$PYTHONPATH:$PLASMABOXDIR/corgi/lib
export PYTHONPATH=$PYTHONPATH:$PLASMABOXDIR/lib
export PYTHONPATH=$PYTHONPATH:$PLASMABOXDIR/python
export PYTHONPATH=$PYTHONPATH:$PLASMABOXDIR/analysis
```
where `path2repo` points to the location where you cloned the repository (i.e. path to plasmabox directory).

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

When compiling and linking is finished, CMake runs few automated tests to check that everything is working. You should see a message *"XX tests finished succesfully"* in the end, if the build was successful.





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

### Physical modules
- `corgi`: corgi library
- `em-fields`: Maxwell field solver module
- `io`: input/output module
- `pic`: Particle-in-Cell module
- `vlasov`: Relativistic Vlasov module
- `radiation`: Radiation module 



