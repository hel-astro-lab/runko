<img align="top" src="docs/header.png">

# Modern C++14/Python3 toolbox for plasma simulations

[![Build Status](https://travis-ci.com/natj/runko.svg?branch=master)](https://travis-ci.com/natj/runko) [![Documentation Status](https://readthedocs.org/projects/runko/badge/?version=latest)](https://runko.readthedocs.io/en/latest/?badge=latest) [![MIT](https://badges.frapsoft.com/os/mit/mit.svg?v=102)](https://github.com/natj/runko/LICENSE) 

Runko is a collection of simulation codes written in modern C++14/Python3 to model astrophysical plasmas. The framework consists of various physical modules that can be run independently or combined to create multi-physics simulations. Low-level "kernels" are mainly implemented in modern C++ that allows to write modular and high performance code. By binding these fast low-level classes to Python objects it is also easy to use and extend them. This ensures efficient code, rapid prototyping, and ease of use.

Under the hood, the framework uses the massively parallel grid infrastructure library [corgi](https://github.com/natj/corgi) that relies on decomposing the grid to smaller subregions, called tiles, that can be operated on and updated independently of each other. [Corgi](https://github.com/natj/corgi) also automatically parallelizes the simulations and provides dynamic load-balancing capability. Therefore, small simulation setups can be tested locally on laptops and then extended for massively parallel supercomputer platforms (currently tested up to ~10k cores).

Documentation is available from [runko.readthedocs.io](https://runko.readthedocs.io/en/latest/?badge=latest) 


## Available modules
Current main physics simulation modules include:
- **FDTD electromagnetic field module** based on staggered Yee lattices (`fields/`)
- **2D3V Particle-In-Cell module** (`pic/`)
- **1D3V Relativistic Vlasov module** (`vlasov/`)

Additionally, modules under construction include:
- **Force-free MHD** module (`ffe/`)
- Non-linear Monte Carlo **radiation module** (`radiation/`)
- Full **3D3V Vlasov module**


## Quick getting started guide
1) Follow the [installation instructions](https://runko.readthedocs.io/en/latest/installation.html) to get Runko running on your laptop.
2) Add your own project repository under `projects/` directory.
	- This can be based on, for example, Python drivers in `projects/tests/`
3) Test & prototype your new simulation setups on your laptop/desktop.
4) Simulate on supercomputers!

## Showcase
<img align="right" src="https://cdn.jsdelivr.net/gh/natj/pb-utilities@master/movies/turb_small.gif">	

### Relativistic kinetic turbulence 	
PIC module has been used to simulate formation of turbulence in collisionless magnetically-dominated pair plasma.

We start by perturbing the initial magnetic field in an otherwise uniform box of plasma. Magnetic eddies are quickly seen to develop. When the eddies collide, thin current sheets are formed. These thin sheets are unstable for magnetic reconnection that starts to tear the sheets and produce plasmoid.

</br>
</br>
</br>

### Collisionless shocks
<img align="center" src="https://cdn.jsdelivr.net/gh/natj/pb-utilities@master/movies/shock.gif">

PIC module has also an initial **work-in-progress** collisionless shock setup. This is based on the common piston setup where plasma is reflected from the left simulation wall and made to collide with it self. A collisionless shock is quickly formed stopping the plasma.


### Plasma instabilities

<img align="right" src="https://cdn.jsdelivr.net/gh/natj/pb-utilities@master/movies/beam.gif">	

Vlasov module has been used to simulate the development of beam instability in a stratified medium. 

In this 1D1V simulation we track the development of stratified beam instability starting from a relative density of ~10^-5 as the incoming beam enters the box. At higher densities of ~1 the beam excites Langmuir wave turbulence and heats the background plasma. This is a density contrast of over 100000!

## How to cite?

You can use the following BibTeX template to cite Runko in any scientific discourse:
```
@ARTICLE{runko,
       author = {{N{\"a}ttil{\"a}}, J.},
        title = "{Runko: Modern multi-physics toolbox for simulating plasma}",
      journal = {arXiv e-prints},
     keywords = {Physics - Computational Physics, Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Plasma Physics},
         year = "2019",
        month = "Jun",
          eid = {arXiv:1906.06306},
        pages = {arXiv:1906.06306},
archivePrefix = {arXiv},
       eprint = {1906.06306},
 primaryClass = {physics.comp-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019arXiv190606306N},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

