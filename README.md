# Kinetic plasma project


## Downloading/cloning

The project depends on submodules so be sure to issue a recursive clone:
```
git clone --recursive https://github.com/natj/plasma.git
```

## Installation

First, you need to compile the libraries. Simple 
```
make
```
should do the trick. In some cases, you have to create your own `makefile` into `archs` directory.


After compilation, you should add the python modules into `PYTHONPATH` environment variable for easier accessing. Modify your `~/.bash_profile` (MacOS) or `~/.bashrc` (linux) by appending `corgi` and `plasma` as
```
export PYTHONPATH=$PYTHONPATH:/path2to/corgi/pycorgi
export PYTHONPATH=$PYTHONPATH:/path2to/plasma/python
```




## Code directory structure
- `notes`: random (latex) notes 
- `projects`: projects files and tests
- `prototypes`: different prototypes and test versions
- `tools`: unix etc. tools




