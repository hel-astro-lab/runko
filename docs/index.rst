.. image:: header.png


Runko --- modern toolbox for plasma simulations 
==========================================================

Runko is a modern numerical toolkit for simulating astrophysical plasmas. It is written in C++17/Python3 and is designed to be highly modular. The name originates from the Finnish word *runko* meaning literally "a frame".  

To get started, follow these manual pages in order. This will get you up-to-date in the prerequisites, installation, and understanding the various algorithms. Notes and tips that focus on code development are listed at the end of the manual.


Physical modules
----------------

The framework consists of various physics modules that can be run independently or combined together to create multi-physics simulations. Different modules include:

* Finite difference time domain electromagnetic module 
* Particle-in-cell (PIC) module
* Force-free magnetohydrodynamics module 
* Relativistic Vlasov module
* Non-linear Monte Carlo Radiation module


Technical goodies
-----------------

* Uses *modern C++14/17* incl. ``std::vector``, ``auto`` keyword, etc.
* Low-level objects are binded to native *Python3* objects via `PyBind11 <https://pybind11.readthedocs.io/en/stable/>`_ providing seamless operability.
* Physics modules are based on *polymorphic classes*; child classes derive their functionality from the base class, occasionally overwriting the original functionality -- no more rewriting your algorithms!
* Same algorithm is specialized to a spesific space dimension during compile time with *template metaprogramming*.
* Algorithms are designed to have abstract *interface classes*; this makes adding new implementations a breeze - just add a new file!


About
-----

This project was originally developed by `Joonas Nättilä <http://natj.github.io/>`_ while in Nordic Institute for Theoretical Physics (NORDITA). 
Key contributors that provided additional features and/or improvements include

* Fredrik Robertsen (CSC)
* John Hope (Univ. Bath)
* Kristoffer Smedt (Univ. Leeds)
* Camilia Demidem (Nordita)
* Maarja Bussov (Univ. Helsinki)
* Alexanda Veledina (Univ. Turku)



.. raw:: html

    <style> .blue {color:blue;} </style>

.. role:: blue


:blue:`Recent changes`
----------------------

* :blue:`20/05/2021 Added a theory section for Vlasov-Maxwell systems and FDTD method`
* :blue:`19/05/2021 Added a unit conversion tool`
* :blue:`15/05/2021 Launched the improved web manual`




------

.. toctree::
   :maxdepth: 2

   installation

.. toctree::
   :caption: Theory:
   :maxdepth: 2

   theory
   units
   plasma


.. toctree::
   :caption: Algorithms:
   :maxdepth: 2

   fld
   pic
   vlv
   ffe
   qed


.. toctree::
   :caption: Tutorials:
   :maxdepth: 2

   shocks
   showcase
   user-defined modules


.. toctree::
   :caption: Developer notes:
   :maxdepth: 2

   unittool
   clusters
   vectorization 
   versions
   style
   debugging
   publications



