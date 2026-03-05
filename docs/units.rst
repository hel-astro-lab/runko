Units
=====

Numerical values
----------------

Runko does not actually use cgs-Gaussian units in the simulation.
Instead, values are normalized with a fiducial values:

.. math::

   t &= \hat t \Delta t

   x &= \hat x \Delta x


Here physical time and position are on the left hand side
and hat denotes a numerical (or code) value used in Runko.

For coordinate velocity:

.. math::

   v = \frac{dx}{dt} = \frac{d\hat x}{d\hat t} \frac{\Delta x}{\Delta t} = \hat v \frac{\Delta x}{\Delta t}


Applying this to speed of light we get numerical speed of light or CFL (Courant–Friedrichs–Lewy number):

.. math::

   \hat c = c \frac{\Delta t}{\Delta x}


and we can write :math:`v = \hat v \frac{c}{\hat c}`.

Rest of the physical values are written as:

.. math::

   B &= \hat B B_0

   E &= \hat E B_0

   q &= \hat q q_0

   m &= \hat m m_0

   J &= \hat J J_0

   A &= \hat A A_0

   n &= \frac{\hat n}{(\Delta x)^3}


Note the same fiducial value for magnetic and electric fields.

Choosing the fiducial values
----------------------------

Derivatives transform as:

.. math::

   \frac{\partial}{\partial t} & = \frac{1}{\Delta t} \frac{\partial}{\partial \hat t}

   \frac{\partial}{\partial x} & = \frac{1}{\Delta x} \frac{\partial}{\partial \hat x}

   \nabla &= \frac{1}{\Delta \hat x}\hat \nabla


Now we can write the dynamical Maxwell's equations as:

.. math::

   \frac{\partial \hat B}{\partial \hat t} &= -\hat c \hat \nabla \times \hat E

   \frac{\partial \hat E}{\partial \hat t} &= +\hat c \hat \nabla \times \hat B
   - \frac{4\pi J_0}{B_0} \Delta t\hat J

and acceleration due to Lorentz force as:

.. math::

   \frac{\partial \hat u}{\partial \hat t} = \frac{\hat c \hat q q_0 B_0}{c \hat m m_0} \Delta t
   \left(\hat E + \frac{\hat v}{\hat c} \times \hat B\right)


By choosing :math:`B_0 = \frac{m_0c}{q_0 \hat c \Delta t}` and :math:`J_0 = \frac{B_0}{4\pi \Delta t}`
the equations above can be written as:

.. math::

   \frac{\partial \hat B}{\partial \hat t} &= -\hat c \hat \nabla \times \hat E

   \frac{\partial \hat E}{\partial \hat t} &= +\hat c \hat \nabla \times \hat B - \hat J

   \frac{\partial \hat u}{\partial \hat t} &= \frac{\hat q}{\hat m}
   \left(\hat E + \frac{\hat v}{\hat c} \times \hat B\right)


Current density is :math:`J = qnv = \frac{q_0}{\Delta t (\Delta x)^2} \hat q \hat n \hat v`.
By demanding that numerical current density can be written as :math:`\hat J = \hat n \hat q \hat v`
we can deduce that

.. math::
   J_0 = \frac{q_0}{\Delta t (\Delta x)^2}


Equating this with :math:`J_0` defined earlier, we can solve for:

.. math::

   \Delta x = 4 \pi \frac{q_0^2}{m_0}\left(\frac{\hat c}{c}\right)


Magnetic field in terms of vector potential is :math:`B = \nabla \times A = \frac{A_0}{\Delta x} \hat \nabla \times \hat A`.
By demanding that numerical magnetic field can be written as :math:`\hat B = \hat \nabla \times \hat A`
we can deduce that:

.. math::

   A_0 = \Delta x B_0
