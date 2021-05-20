.. default-role:: math

Electromagnetic fields
=============================

Electromagnetic Maxwell's fields are propagated with a Finite-difference-time-domain (FDTD) method.

Grid staggering
---------------

Yee lattice is used for `\mathbf{E}` and `\mathbf{B}` fields such that they are staggered both in space and in time,

.. math::

    \hat{\mathbf{E}} = \left(
    \hat{E}_{x;\, i+\frac{1}{2}, j    , k    },\,
    \hat{E}_{y;\, i    , j+\frac{1}{2}, k    },\,
    \hat{E}_{z;\, i    , j    , k+\frac{1}{2}},
    \right)^{n} 

    \hat{\mathbf{B}} = \left(
    \hat{B}_{x;\, i    , j+\frac{1}{2},  k+\frac{1}{2}},\,
    \hat{B}_{y;\, i+\frac{1}{2}, j    ,  k+\frac{1}{2}},\,
    \hat{B}_{z;\, i+\frac{1}{2}, j+\frac{1}{2},  k    },
    \right)^{n+\frac{1}{2}},


where `\hat{\mathbf{E}}` is located in the middle of the cell sides and `\hat{\mathbf{B}}` in the center of the cell faces. 
This makes it easy to calculate the curl of the fields in the following equations.


FDTD method
-----------

In the time domain we update `\hat{\mathbf{E}}` and `\hat{\mathbf{B}}` fields with discrete forms of Maxwell's equations:

.. math::

    \Delta[ \hat{\mathbf{E}} ]_t = \hat{c} \Delta[ \hat{\mathbf{B}} ]_{\mathbf{x}} - \hat{J},

and

.. math::

    \Delta[ \hat{\mathbf{B}} ]_t =-\hat{c} \Delta[ \hat{\mathbf{E}} ]_{\mathbf{x}},

where `\Delta[Q]_{t,\mathbf{x}}` is the time differentiation or curl of a variable `Q` without the `\Delta x` or `\Delta t` denominator.

Only normalization factor entering these equations is the Courant velocity, `\hat{c}`, since everything else is absorbed in `B_0` and `J_0`.
There is no need to solve divergence equations if charge-conserving scheme is used together with the Yee staggering of the fields.


Finite-difference solvers
-------------------------

.. note::

    TODO


Filtering
---------

.. note::

    TODO

