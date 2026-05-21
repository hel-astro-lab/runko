.. default-role:: math

Vlasov-Maxwell equations
========================

.. note::

    This discussion is based on the `Nättilä (2019) <https://arxiv.org/abs/1906.06306>`_ paper, where an extended discussion can be found.


Let us discuss the special-relativistic formulation of the Vlasov/Boltzmann equation.
The spatial coordinate location vector is `\mathbf{x} \equiv (x,y,z)` and coordinate time is denoted by `t`.
The coordinate velocity (three-velocity) is `\mathbf{v} \equiv d_t \mathbf{x}`, and its individual Cartesian components are written as `\mathbf{v} = (v_x, v_y, v_z)`.
The proper time (measured with a co-moving clock), `\tau`, is related to the coordinate time through the Lorentz factor `\gamma \equiv d_{\tau} t`.
The proper velocity (the spatial components of the four-velocity) is `\mathbf{u} \equiv d_{\tau} \mathbf{x} = \gamma \mathbf{v}`.
The Lorentz factor and the velocities are related by the expression

.. math::

    \gamma^2 = 1 + (u/c)^2 = (1-(v/c)^2)^{-1},

where `c` is the speed of light, `u = |\mathbf{u}|` and `v = |\mathbf{v}|`.
Acceleration is denoted with `\mathbf{a} \equiv d_t \mathbf{u}`.

The six-dimensional phase-space density distribution for particle species s is `f_s \equiv f_s(\mathbf{x}, \mathbf{u}; t)`.
Thus, `f_s\, d^3 x \, d^3 u` is the number of particles in the six-dimensional differential phase-space volume between `\mathbf{x}`, `\mathbf{u}` and `\mathbf{x} + d\mathbf{x}`, `\mathbf{u} + d\mathbf{u}`.

The evolution of `f_s` is governed by the relativistic Boltzmann/Vlasov equation

.. math::

    \frac{\partial f_s}{\partial t} + \mathbf{v} \cdot \nabla_{\mathbf{x}} f_s + \mathbf{a}_{\mathrm{s}} \cdot \nabla_{\mathbf{u}} f_s  = C,

where `\nabla_{\mathbf{x}} = \frac{d}{d \mathbf{x}}` and `\nabla_{\mathbf{u}} = \frac{d}{d \mathbf{u}}` are the spatial and momentum parts of the differential operator `\nabla`, respectively.
The term on the right-hand side, defined as `C \equiv \partial_t f_s ~|_{\mathrm{coll}}`, is the collision operator.
For the Vlasov equation `C = 0`, i.e., the plasma is collisionless.

The acceleration of a charged particle, `\mathbf{a}_{\mathrm{s}}`, is governed by the Lorentz force

.. math::
    \mathbf{a}_{\mathrm{s}} \equiv \mathrm{d}_t \mathbf{u} = \frac{q_{\mathrm{s}} }{  m_{\mathrm{s}} } \left(\mathbf{E} + \frac{\mathbf{v}}{c} \times \mathbf{B} \right)
   = \frac{q_{\mathrm{s}} }{  m_{\mathrm{s}} } \left(\mathbf{E} + \frac{\mathbf{u}}{\gamma c} \times \mathbf{B} \right),

where `\mathbf{E}` and `\mathbf{B}` are the electric and magnetic fields, `q_{\mathrm{s}}` is the charge, and `m_{\mathrm{s}}` is the mass of the particle of species s.

Moments of the distribution function define macroscopic (bulk) quantities of the plasma.
The zeroth moment of the distribution function `f_s` defines the charge density of species s as

.. math::

    \rho_{\mathrm{s}} \equiv  q_{\mathrm{s}} \int f_s \, \mathrm{d} \mathbf{u}.

The total charge density is `\rho = \sum_{\mathrm{s}} \rho_{\mathrm{s}}`.
The first moment defines the current (charge flux) vector as

.. math::

    \mathbf{J}_{\mathrm{s}} \equiv q_{\mathrm{s}} \int f_s \mathbf{v} \, \mathrm{d} \mathbf{u} 
                      = q_{\mathrm{s}} \int f_s \, \frac{ \mathbf{u}}{\gamma} ~\mathrm{d} \mathbf{u}.

The total current is obtained by summing over all plasma species, `\mathbf{J} = \sum_{\mathrm{s}} \mathbf{J}_{\mathrm{s}}`.


Maxwell's equations
-------------------

The evolution of the electric field `\mathbf{E}` and magnetic field `\mathbf{B}` is governed by Maxwell's equations.
These are Gauss's law,

.. math::

    \nabla \cdot \mathbf{E} = 4\pi \rho,

Gauss's law for magnetism,

.. math::

    \nabla \cdot \mathbf{B} = 0,

Faraday's law,

.. math::

    \nabla \times \mathbf{E} = -\frac{1}{c}\frac{\partial \mathbf{B}}{\partial t},

and Ampère's law,

.. math::

    \nabla \times \mathbf{B} = \frac{4\pi}{c}\mathbf{J} +\frac{1}{c}\frac{\partial \mathbf{E}}{\partial t}.

Charge conservation follows from these by taking the divergence of Ampère's law and substituting Gauss's law to get

.. math::

    \frac{\partial \rho}{\partial t} + \nabla \cdot \mathbf{J} = 0.







