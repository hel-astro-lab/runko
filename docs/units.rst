Units
==========

We typically use four fundamental quantities to characterize the simulations:

* numerical speed of light, :math:`\hat{c}` (``cfl`` or ``c``)
* particle number density, :math:`\hat{n_{s}}` (``ppc``), 
* numerical skin-depth resolution, :math:`\hat{R}` (``c_omp``), and
* plasma magnetization, :math:`\sigma` (``sigma``).


.. note::
    Technically we are using the cgs-Gaussian unit system in the following discussions. This is an unrationalized unit system meaning that factors of :math:`4\pi` appear in the Maxwell's equations.


Electromagnetic field
---------------------------

Faraday's and Ampere's laws are

.. math::
    \frac{\partial E}{\partial t} &=+c \nabla \times B - 4\pi J

    \frac{\partial B}{\partial t} &=-c \nabla \times E

that in code units corresponds to solving

.. math::
    \Delta_t \hat{E} &=+\hat{c} \nabla \times \hat{B} - \hat{c} \hat{J}

    \Delta_t \hat{B} &=-\hat{c} \nabla \times \hat{E}


The benefits are clear:
the equations appear symmetric and only have the numerical speed of light, :math:`\hat{c}` as a natural constant in them.



Lorentz force 
-------------------

The Lorentz force experienced by charged particles is


.. math::
    \frac{d (\gamma \beta) }{d t} =  \frac{q_s}{m_s} (E + \beta \times B)

that in code units corresponds to solving

.. math::
    \Delta_t (\gamma \beta) = \hat{c} \frac{\hat{q}}{\hat{m}} (\hat{E} + \beta \times \hat{B})

.. note::
    Only the charge-to-mass ratio enters the Lorentz force.


Current density 
---------------------

Charge current density and its numerical counterparts are

.. math::
    4\pi J = 4\pi q \beta c 
    \quad\mathrm{and}\quad
    \hat{J}\Delta t = \hat{q} \beta \hat{c}^2

Current density is technically a flux of charge through an area per time;
therefore, the current has numerical units of :math:`\hat{J} = \frac{v \hat{q}}{\Delta x^2 \Delta t}`.


.. note::
    Technically we store the current density per time step in the grid; hence the extra :math:`\hat{c} = \Delta t` factor in the expression of the numerical charge current density.


