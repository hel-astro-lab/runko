.. image:: header.png


Runko --- modern toolbox for plasma simulations 
==========================================================

Runko is a modern numerical toolkit for simulating astrophysical plasmas.
It is written in C++ & Python and is designed to be highly modular.
The name originates from the Finnish word *runko* meaning literally "a frame".

.. Warning::

   This documentation describes Runko v5, which contains major breaking changes compared to v4.
   Runko v5 comes with a HIP-based GPU backend that required significant refactoring and API changes.
   Some features present in v4 have not (yet) been ported to v5.


Physical modules
----------------

The framework consists of various physics modules that can be run independently
or combined to create multi-physics simulations. The modules currently include:

* Finite-difference time-domain electromagnetic module
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
   configuration-parameters

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

   contributing
   dev-workflow
   debugging
   unittool
   versions
   documentation
