.. default-role:: math

.. _emf-link:

Electromagnetic fields
=============================

Electromagnetic Maxwell's fields are propagated with a Finite-difference-time-domain (FDTD) method.

.. note::

    Due to the choice of our unit system, the time step `\Delta t` is basically equivalent to the numerical speed of light, `\hat{c}`, in the following discussion;
    therefore, when you see `\hat{c}` in the following formulas, you can mentally replace it by `\Delta t`.

Grid staggering
---------------

Yee lattice is used for `\mathbf{E}` and `\mathbf{B}` fields such that they are staggered both in space and in time,

.. math::

    \hat{\mathbf{E}} &= \left(
    \hat{E}_{x;\, i+\frac{1}{2}, j    , k    },\,
    \hat{E}_{y;\, i    , j+\frac{1}{2}, k    },\,
    \hat{E}_{z;\, i    , j    , k+\frac{1}{2}},
    \right)^{n} 

    \hat{\mathbf{B}} &= \left(
    \hat{B}_{x;\, i    , j+\frac{1}{2},  k+\frac{1}{2}},\,
    \hat{B}_{y;\, i+\frac{1}{2}, j    ,  k+\frac{1}{2}},\,
    \hat{B}_{z;\, i+\frac{1}{2}, j+\frac{1}{2},  k    },
    \right)^{n-\frac{1}{2}},


where `\hat{\mathbf{E}}` is located in the middle of the cell sides and `\hat{\mathbf{B}}` in the center of the cell faces.
This makes it easy to compute the curl of the fields in the equations below.
Note that there is at most one component at any grid location.
Therefore, we can drop the explicit component lables :math:`{x, y, z}` from the notation and for example simply write
:math:`B_{z;\, n+\frac12, m+\frac12, l} = B_{n+\frac12, m+\frac12, l}`
and :math:`E_{x;\, n+\frac12, m, l} = E_{n+\frac12, m, l}`.



The time staggering, in turn, increases the temporal order of the scheme since all the updates are leapfrog-like, `x^{n+1} = x^n + v^{n+\frac{1}{2}} \Delta t`.

FDTD method
-----------

In the time domain we update `\hat{\mathbf{E}}` and `\hat{\mathbf{B}}` fields with discrete forms of Maxwell's equations:

.. math::

    \Delta[ \hat{\mathbf{E}} ]_t = \hat{c} \Delta[ \hat{\mathbf{B}} ]_{\mathbf{x}} - \hat{c} \hat{J},

and

.. math::

    \Delta[ \hat{\mathbf{B}} ]_t =-\hat{c} \Delta[ \hat{\mathbf{E}} ]_{\mathbf{x}},

where `\Delta[Q]_{t,\mathbf{x}}` is the time differentiation or curl of a variable `Q` without the `\Delta x` or `\Delta t` denominator.

The only normalization factor entering these equations is the Courant velocity, `\hat{c}`, since everything else is absorbed into `B_0` and `J_0`.
There is no need to solve the divergence equations if a charge-conserving scheme is used together with the Yee staggering of the fields.


Finite-difference solvers
-------------------------

Different FDTD solvers can be obtained by defining different discrete curl-operators,  `\nabla \times \mathbf{Q} \rightarrow \Delta[ \hat{\mathbf{Q}} ]_{\mathbf{x}}`, i.e., different finite-difference operators for the calculation of the derivative.

FDTD2
.....

Second order accurate FDTD uses:

.. math::
   \left.\Delta[ \hat{\mathbf{Q}} ]_{\mathbf{x}}\right|_{n,m,l} =
   \begin{bmatrix}
     \left(Q_{n,m+\frac12,l} - Q_{n,m-\frac12,l}\right)
     - \left(Q_{n,m,l+\frac12} - Q_{n,m,l-\frac12}\right) \\
     \left(Q_{n,m,l+\frac12} - Q_{n,m,l-\frac12}\right)
     - \left(Q_{n+\frac12,m,l} - Q_{n-\frac12,m,l}\right) \\
     \left(Q_{n+\frac12,m,l} - Q_{n-\frac12,m,l}\right)
     - \left(Q_{n,m+\frac12,l} - Q_{n,m-\frac12,l}\right)
   \end{bmatrix}


For more detailed derivation, see
`Miro Palmu's master's theses <http://hdl.handle.net/10138/627079>`_.

Filtering
---------

In PIC, particles generate noisy current density to the grid.
This noise can be smoothed by applying digital filters.
Compared to using more particles or using higher order methods,
filtering is computationally cheaper.

.. admonition:: WIP

   Expand this section.
