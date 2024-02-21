.. default-role:: math

Particle-in-cell 
================

Particle-in-cell (PIC) method can be used to solve the coupled Vlasov-Maxwell equations by sampling the phase space density `f_s(\mathbf{x}, \mathbf{u}; t)` with computational particles,

.. math::

    f_s(\mathbf{x}, \mathbf{u}; t) = \sum_p w_p ~S[\mathbf{x} - \mathbf{x}_p(t)] ~\delta[\mathbf{u} - \mathbf{u}_p(t)]

where `w_p` is the weight of the `p`:th particle,
`S[x]` is the particle's shape function (defining literally its spatial shape and extent),
and `\delta[x]` is the Dirac's delta function.
The computational particle is therefore assumed to be a monolithic collection of "real" charged particles that move in unison.

The computational particles interact with electromagnetic fields, that are defined on a grid, via the Lorentz force.
The electromagnetic fields are advanced self-consistently by the electric currents, `J = \sum_p w_p q_e \beta_p`, induced by the moving particles via the Ampere's and Faraday's laws. 



Initialization
--------------

In our typical PIC simulations the pair plasma particles, `e^-` and `e^+`, are initialized exactly on top of each other.
This means that the initial charge density is `\rho = 0`;
therefore, there is no need to solve the `\nabla \cdot \mathbf{\hat{E}} = 4\pi \rho` in any stage.

.. note::

    Cold mobile ions can be implemented into the simulation by having `\rho \ne 0` (i.e., not initializing particles such that their charge cancels);
    after this, the simulation behaves as if there are "ghost charges" (that do not evolve) on the grid.



Time advancement
----------------

The time-discretized variables are defined on a time step `n` such that for a general variable `Q` we write `Q^{(n)} \equiv Q(n\Delta t)`.
In the beginning of the time loop, 

.. math::

    \mathbf{\hat{x}} &\rightarrow \mathbf{\hat{x}}^{n}

    \mathbf{\hat{u}} &\rightarrow \mathbf{\hat{u}}^{n-\frac{1}{2}}

    \mathbf{\hat{E}} &\rightarrow \mathbf{\hat{E}}^{n}

    \mathbf{\hat{B}} &\rightarrow \mathbf{\hat{B}}^{n-\frac{1}{2}}

All of the updates proceed in a leapfrog-fashion such that a quantity `Q^n` is always updated with a derivative at `\Delta Q^{n+\frac{1}{2}}` making the time advancement second-order in time.


1.  Advance `\mathbf{\hat{B}}` half a time step forward according to Faraday's law.

.. math::

    \mathbf{\hat{B}}^{n} = \mathbf{\hat{B}}^{n-\frac{1}{2}} + \frac{1}{2} \hat{c} \Delta_t[\mathbf{\hat{B}}]


2. Interpolate the `\mathbf{\hat{E}}^n` and `\mathbf{\hat{B}}^n` to the particles' location `\mathbf{\hat{x}}^n`.

3. Compute Lorentz force on particles. This corresponds to updating the particle's velocity and location with a pusher, that computes the velocity in the middle of the step, `\mathbf{\hat{u}}^n`, so that 

.. math::

    \mathbf{\hat{u}}^{n+\frac{1}{2}} &= \mathbf{\hat{u}}^{n-\frac{1}{2}} + \hat{c} \frac{q_p}{m_p} 
    \left( \mathbf{\hat{E}}^n + \frac{ \mathbf{\hat{u}}^n }{\gamma^n} \times \mathbf{\hat{B}}^n \right)

    \mathbf{\hat{x}}^{n+1}           &= \mathbf{\hat{x}}^{n} + \hat{c} \frac{ \mathbf{\hat{u}}^n }{ \gamma^n}


4. Advance the second half of the `\mathbf{\hat{B}}` field according to the Faraday's law.

.. math::

    \mathbf{\hat{B}}^{n+\frac{1}{2}} = \mathbf{\hat{B}}^{n} + \frac{1}{2} \hat{c} \Delta_t[\mathbf{\hat{B}}]

5. Advance a full time step of the `\mathbf{\hat{E}}` field according to the Ampere's law.

.. math::

    \mathbf{\hat{E}}^{n+1} = \mathbf{\hat{E}}^{n} + \hat{c} \Delta_t[\mathbf{\hat{E}}]

6. Compute the electric current by the moving particles, `\mathbf{\hat{x}}^n \rightarrow \mathbf{\hat{x}}^{n+1}`, and deposit it into the middle of the trajectory, `\mathbf{\hat{x}}^{n+\frac{1}{2}}`, as

.. math::

    \mathbf{\hat{J}}^{n+\frac{1}{2}}( \mathbf{\hat{x}}^{n+\frac{1}{2}} ) = \sum_p q_e \mathbf{\hat{u}}^{n+\frac{1}{2}}


7. Optionally, smooth the currents with a smoothing operation, `S`,

.. math::

    \mathbf{\hat{J}}' = S[ \mathbf{\hat{J}} ]

8. Add the currents to the electric field,

.. math::

    \mathbf{\hat{E}}^{, n+1} = \mathbf{\hat{E}}^{n+1} - \hat{c} \mathbf{\hat{J}}^{n+\frac{1}{2}}


And we are finished! All the quantities are now defined at

.. math::

    \mathbf{\hat{x}} &\rightarrow \mathbf{\hat{x}}^{n+1}

    \mathbf{\hat{u}} &\rightarrow \mathbf{\hat{u}}^{n+\frac{1}{2}}

    \mathbf{\hat{E}} &\rightarrow \mathbf{\hat{E}}^{n+1}

    \mathbf{\hat{B}} &\rightarrow \mathbf{\hat{B}}^{n+\frac{1}{2}}

    \mathbf{\hat{J}} &\rightarrow \mathbf{\hat{J}}^{n+\frac{1}{2}}

and we can move on to the next time step.


Particle shapes
---------------

Electromagnetic field interpolation from the grid to the arbitrary particle location and the particle charge current projection from an arbitrary particle location to the electromagnetic grid are complementary operations.

Phase density `f` of `N_p` elementary particles (with location `\mathbf{x}` and coordinate-velocity `\mathbf{v}`) in PIC simulations is represented as 

.. math::

    f(\vec{x},\vec{v}) = N_p S(\vec{x}-\vec{x}_p) \cdot \delta(\vec{v} - \vec{v}_p).

Here `S` is the spatial shape function of the computational macro particle (representing `N_p` elementary particles).
The macro particle has a total charge `N_p q_s` and a total mass `N_p m_s`.

The shape function can have an extended spatial footpoint;
the computational particle can be thought of as a "cloud" covering many grid cells.
The fraction of a macro particle in a cell, `x, x+\Delta x`, is described by a weight function `W`

.. math::

    W(x) = \int_{x-\frac{1}{2}\Delta x}^{x+\frac{1}{2}\Delta x} S(x - x') dx'


The simplest shape function is the delta function `S_0(x) = \delta(x)`.
It, however, causes rapid bursts of electric current in the grid when a particle moves from one cell to another;
these current bursts drive unphysical high-frequency electromagnetic radiation.
Instead, we use higher-order B-splines as particle shapes.

.. note::
    
    Staggering of electromagnetic fields needs to be taken into account when interpolating and projecting. This somewhat complicates the above expressions because each field component and directions needs to be dealt independently.

Field interpolation
^^^^^^^^^^^^^^^^^^^

Electromagnetic fields are, in general, interpolated from the grid points, `\mathbf{x}_c`, to the particle's location, `\mathbf{x}_p`, as

.. math::

    \mathbf{E}_p &= \mathbf{E}(\mathbf{x}_p) = \int d\mathbf{x} ~ \mathbf{E}(\mathbf{x}) S(\mathbf{x}-\mathbf{x}_p)

    \mathbf{B}_p &= \mathbf{B}(\mathbf{x}_p) = \int d\mathbf{x} ~ \mathbf{B}(\mathbf{x}) S(\mathbf{x}-\mathbf{x}_p)

Using the connection between the shape and weight functions, this reduces to a summation operation for a regularly-spaced electromagnetic field lattice,

.. math::

    \mathbf{E}_p = \sum_{\mathbf{x}_c} \mathbf{E}(\mathbf{x}_c) W(\mathbf{x}_c - \mathbf{x}_p)

    \mathbf{B}_p = \sum_{\mathbf{x}_c} \mathbf{B}(\mathbf{x}_c) W(\mathbf{x}_c - \mathbf{x}_p)


Current projection
^^^^^^^^^^^^^^^^^^

The projection of a charge current into `\mathbf{\hat{x}}`, `\mathbf{\hat{y}}`, and `\mathbf{\hat{z}}` directions is an opposite operation involving both 
particle weight `W` to the direction of movement, and
particle shape `S` to the non-varying direction, 
expressed as

.. math::

    \hat{j}_x &= q \hat{v}_x W_x(x_c -x) ~ S_y (y_c - y_p ) ~ S(z_c - z)

    \hat{j}_y &= q \hat{v}_y S_x(x_c -x) ~ W_y (y_c - y_p ) ~ S(z_c - z)
    
    \hat{j}_z &= q \hat{v}_z S_x(x_c -x) ~ S_y (y_c - y_p ) ~ W(z_c - z)

Total current is a sum over all the particles `p`, and species `s`, `\hat{J}_i = \sum_s \sum_p \hat{j}_{i; ~s, p}`, and the total current vector is, `\mathbf{\hat{J}} = (\hat{J}_x, \hat{J}_y, \hat{J}_z)`.



Particle pushing
----------------

.. note:: 
    TODO






