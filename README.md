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
    
## Quick installation on macOS
For easy install you need `homebrew` (see http://brew.sh).
**Note:** this does not support threading but can be used to debug the code.

```
brew tap homebrew/science
brew install openmpi
brew install eigen
brew install boost
brew install hdf5 --with-cxx
brew install zoltan
```

In addition create `libs` folder to your home-dir and
```
git clone https://github.com/fmihpc/dccrg.git
```

