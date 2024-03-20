.. default-role:: math

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
This makes it easy to calculate the curl of the fields in the following equations.

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

Only normalization factor entering these equations is the Courant velocity, `\hat{c}`, since everything else is absorbed in `B_0` and `J_0`.
There is no need to solve divergence equations if charge-conserving scheme is used together with the Yee staggering of the fields.


Finite-difference solvers
-------------------------

Different FDTD solvers can be obtained by defining different discrete curl-operators,  `\nabla \times \mathbf{Q} \rightarrow \Delta[ \hat{\mathbf{Q}} ]_{\mathbf{x}}`, i.e., different finite-difference operators for the calculation of the derivative.

.. note::

    TODO

.. literalinclude:: ../em-fields/propagator/fdtd2.c++
   :linenos:
   :caption:
   :dedent:
   :language: c++
   :start-after: SPHINX emf docs pusher example start
   :end-before: SPHINX emf docs pusher example stop



Filtering
---------

.. note::

    TODO

