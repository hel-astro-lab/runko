# Installing OpenMPI on MacOS

For convenience it might be better to compile MPI library yourself so that you get the MPI compiler wrappers together with new C++14 compliant compiler. Similarly, in Mojave this is needed because Apple changed header file location and so we get `no mpi.h header` errors.

Before beginning, uninstall old openmpi to avoid confusion of having multiple versions.
```
brew uninstall openmpi
```

## Compile and configure openmpi

Create a library folder into your home directory (~/) 
```
mkdir libs
cd libs
```

Download newest open-mpi (v.4.0 as of writing) here https://www.open-mpi.org/software/ompi/v4.0/ (.tar.gz) and unzip:
```
tar xf openmpi-4.0.0.tar.gz
```

Configure to use gcc-7 and to install into ~/libs
```
./configure CC=gcc-7 CXX=g++-7 --prefix=~/libs
```

Compile it
```
make all
```

Install
```
make install
```

Next, setup the new compilers and MPI in your `~/.bash_profile` as

Export the libs directory. Note the order of export.
```
export PATH=$HOME/libs/include:${PATH}
export PATH=$HOME/libs/bin:${PATH}
```

And set `mpic++` and `mpicc` wrappers as the default compilers:
```
export CXX=mpic++
export CC=mpicc
```


## reinstall mpi4py

Mpi4py needs to be reinstalled to use the newly compiled openmpi. This is non-trivial because pip tries to re-use old files as much as possible. Run:
```
pip3 uninstall mpi4py
pip3 install mpi4py --no-cache-dir
```

## Re-compile PlasmaBox

Now everything should be working. Normally, in build directory run
```
cmake ..
make -j4
```


