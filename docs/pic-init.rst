.. _pic-init-tutorial:

Initializing a PIC simulation
#############################

This tutorial walks through how to calculate initial values for a PIC simulation.
Prerequisites for this are the :ref:`units <units-theory>`
and :ref:`palsma parameters <plasma-params-theory>` theory sections.

Mass-to-charge ratio
====================

For particle species ``s`` we denote its numerical mass-to-charge ratio
as :math:`r_s = \frac{\hat m_s}{|\hat q_s|}`.
We can now write the numerical total plasma plasma frequency as:

.. math::

   \frac{\hat \omega_p^2}{\Delta t^2} = \frac{\sum_s \hat \omega_{p,s}^2}{\Delta t^2}
   = \sum_s \frac{4\pi \hat n_s \hat q_s^2 q_0^2}
   {\Delta x^3 \gamma_s r_s|\hat q_s| m_0}

Inserting the expression for :math:`\Delta x` from the units chapter
and then simplifyinig leads to:

.. math::

   \hat \omega_p^2 = \sum_s \frac{\hat n_s |\hat q_s|} {\gamma_s r_s }


Satisfying Gauss's law [for magnetism]
======================================

By using charge conserving current depositors,
if Guass's law [for magnetism] is satisfied at the beginning
then it is satisfied at all times.

.. note::
   Cold mobile ions can be implemented into the simulation by having `\rho \ne 0`
   (i.e., not initializing particles such that their charge cancels);
   after this, the simulation behaves as if there are "ghost charges"
   (that do not evolve) on the grid.

Pair-plasma
===========

We have two species labeled with ``e`` and ``i``.
We will assume following:

.. math:: |\hat q_e| = |\hat q_i| \equiv |\hat q|
          \quad \text{and} \quad \gamma_e = \gamma_i \equiv \gamma

Therefore :math:`\frac{\hat m_e}{r_e} = \frac{\hat m_i}{r_i}`.
In order for the plasma to be in quasineutrality we also have
:math:`\hat n_e = \hat n_i \equiv \hat n`.

Total plasma frequency can thus be written as:

.. math:: \hat \omega_p^2 = \frac{\hat n |\hat q|}{\gamma}\sum_s \frac1{r_s}

Now by choosing :math:`|\hat q_0|` and :math:`m_0` such that
:math:`r_e = \frac{\hat m_e }{| \hat q_e|} = 1` we get that
:math:`r_i = \frac{\hat m_i}{\hat m_e}` and totale plasma frequency becomes:

.. math:: \hat \omega_p^2 = \frac{\hat n |\hat q|}{\gamma}
          \left(1 + \frac{\hat m_e }{\hat m_i}\right)

which can be solved for :math:`|\hat q|`:

.. math:: |\hat q| = \frac{\hat \omega_p^2 \gamma}
          {\hat n \left(1 + \frac{\hat m_e }{\hat m_i}\right)}

This can be used to specify Runko parameters ``q0`` and ``q1``
when total plasma skin depth resolution :math:`R` is chosen
and substituted to :math:`\hat \omega_s = \frac{\hat c}{R}`.
