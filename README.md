# radpic

General relativistic radiative hybridly parallelized particle-in-cell code with adaptive mesh and dynamic load balancing.


### Physical notes (and also TODO)
    * (general) relativistic
    * radiation 


### Computational notes
    * Domain decomposed with MPI
    * Threaded using openMP
    * Vectorized with Eigen 
    * Adaptive mesh refinement with object oriented cell platform (dccrg)
    * Dynamically load balanced with hypergraph slicing (Zoltan)
    

## Installation

To compile one typically just needs to run `make`. However, be sure to first set up a proper environment in the `Makefile` (there are options for osx, osx-mpi, gnu...).

Dependencies are 
    * boost (c++ extensions)
    * Eigen (vector calculus)
    * dccrg (adaptive mesh)
    * Zoltan (dynamical load balancing)
    * hdf5 (optional)


### dccrg
dccrg is a header only library so one only needs to run

```
git clone https://github.com/fmihpc/dccrg.git
```

See also https://github.com/fmihpc/dccrg/wiki/Install


### Zoltan
http://www.cs.sandia.gov/zoltan/ug_html/ug_usage.html

```
CXX="mpic++" ../configure --prefix=/U
sers/natj/Documents/rndm/Zoltan --enable-mpi --with-mpi-compilers --with-id-typ
e=ullong
make
make install
``` 

### Boost


```
wget -c 'https://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.zip/download'
tar xf download
./boostrap.sh --prefix=/home/jatnat/lib
./b2 install
```

Include support for mpi:

```
echo "using mpi ;" >> ./tools/build/src/user-config.jam
./b2
./b2 --prefix=$HOME install
```

### Eigen 
Eigen is a header only and is only needed in the linking phase.
http://eigen.tuxfamily.org/

In case you want to compile (some parts) and copy (all):

```
mkdir build_dir
cd build_dir
cmake /home/jatnat/sources/eigen-... -DCMAKE_INSTALL_PREFIX=/home/jatnat/lib -DINCLUDE_DIR=/home/jatnat/lib
make install
```


### Additional info about MPI + OSX
```
./configure --prefix=/usr/local/openmpi
make
sudo make install
```

#MPI stuff (compiled from source)
```
export MPI_DIR=/usr/local/openmpi
export PATH=/usr/local/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH
```


### HDF5

```
http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.16.tar.gz
./configure --prefix=/home/jatnat/lib --enable-fortran -enable-cxx
#--with-zlib=/home/jatnat/lib,/home/jatnat/lib/include #obselete if szlib
--with-szlib=/home/jatnat/lib,/home/jatnat/lib/include
make
make check
make install
```





