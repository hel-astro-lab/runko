# Full installation with openMP + MPI

**Note:** Currently not working.



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

### Additional info about HPC in OSX
macOS does not have recent enough version of c++ compiler to be used with threading (openMP).
Options are to install `llvm` if you want to keep using `clang` or `gcc` if you prefer to use the gnu compiler.
Both options are valid.

To install you need `homebrew` (see http://brew.sh)

For `gcc`
```
brew install gcc --without-multilib
```
after which the compiler can be used with `g++-6`.

Or for `clang`
```
brew install llvm
```
where the new compiler is now in `/usr/local/opt/llvm/bin/clang++`.
Additionally, you have to append the library path with `-L/usr/local/opt/llvm/lib´ with `clang`.

Since 2016, both support threading out-of-the-box.

For message passing we need MPI compiled with aforementioned proper compiler with threading capability.
```
brew install openmpi --build-from-source --cc=gcc-6
```



```
./configure --prefix=/usr/local/openmpi
make
sudo make install
```

MPI stuff (compiled from source)
```
export MPI_DIR=/usr/local/openmpi
export PATH=/usr/local/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH
```




### vlsv+VisIt+silo
```
brew install homebrew/science/silo
git clone https://github.com/fmihpc/vlsv.git
```
Copy `tools/vlsv/Makefile.macos` to the vlsv folder.
then compile with
```
make "ARCH=macos"
```

### VisIt
Install visit.app bundle
go to /Applications/Visit/.../
update -L and -I paths to /Users/XXX/libs/vlsv



sudo xcode-select -switch /Library/Developer/CommandLineTools

build_visit2_12_0.sh:19270
    export VTK_VERSION=${VTK_VERSION:-"7.1.0"}
    export VTK_SHORT_VERSION=${VTK_SHORT_VERSION:-"7.1"}

sudo ln -s /usr/local/Cellar/qt5/5.7.1_1/mkspecs /usr/loc
al/mkspecs

sudo ln -s /usr/local/Cellar/qt5/5.7.1_1/plugins /usr/loc
al/plugins



