.. image:: header.png


Runko --- modern toolbox for plasma simulations 
==========================================================

Runko is a modern numerical toolkit for simulating astrophysical plasmas.
It is written in C++ & Python and is designed to be highly modular.
The name originates from the Finnish word *runko* meaning literally "a frame".

.. Warning::

   This documentation is about Runko v5 which contains major breaking changes compared to v4.
   Runko v5 comes with HIP based GPU backend which required significant refactors and API changes.
   Some features present in v4 are not (yet) ported for v5.


Physical modules
----------------

The framework consists of various physics modules that can be run independently
or combined together to create multi-physics simulations. Different modules include:

* Finite difference time domain electromagnetic module
* Particle-in-cell (PIC) module


About
-----

This project was originally developed by `Joonas Nättilä <http://natj.github.io/>`_ while in Nordic Institute for Theoretical Physics (NORDITA). 
Key contributors that provided additional features and/or improvements include

* Fredrik Robertsen (CSC)
* John Hope (Univ. Bath)
* Kristoffer Smedt (Univ. Leeds)
* Camilia Demidem (Nordita)
* Maarja Bussov (Univ. Helsinki)
* Alexandra Veledina (Univ. Turku)
* Miro Palmu (Univ. Helsinki)



.. raw:: html

    <style> .blue {color:blue;} </style>

.. role:: blue


------

.. toctree::
   :maxdepth: 2

   installation

.. toctree::
   :caption: Tutorials:
   :maxdepth: 2

   pic-tutorial

.. toctree::
   :caption: Theory:
   :maxdepth: 2

   theory
   units
   plasma
   emf
   pic


.. toctree::
   :caption: Developer notes:
   :maxdepth: 2

   dev-workflow
   unittool
   versions
   documentation
