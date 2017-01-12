# radpic

General relativistic radiative hybridly parallelized particle-in-cell code with adaptive mesh and dynamic load balancing.

### TODO
    * typedefs and defines to define.h
    * vlsv file format
    * refactor
    * Doxygen

### Physical notes (and also TODO)
    * (general) relativistic
    * radiation 


### Computational notes
    * Domain decomposed with MPI
    * Threaded using openMP
    * Vectorized with Eigen 
    * Adaptive mesh refinement with object oriented cell platform (dccrg)
    * Dynamically load balanced with hypergraph slicing (Zoltan)
    
## Quick installation on macOS
For easy install you need `homebrew` (see http://brew.sh).
**Note:** this does not support threading but can be used to debug the code.

### Trivial dependencies
```
brew tap homebrew/science
brew install openmpi
brew install eigen
brew install boost
brew install zoltan
```

### dccrg mpi-grid
In addition create `libs` folder to your home-dir and inside there run
```
git clone https://github.com/fmihpc/dccrg.git
```
to get the dccrg grid library.

### SSE/AVX vectorization
C++ vector class library is similarly a header-only library and is obtainable from
```
http://www.agner.org/optimize/
```
It should be extracted to the lib dir into `vectorclass` folder.


### vlsv+VisIt+silo
**Work in progress**


