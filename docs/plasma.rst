Plasma parameters
=================

Here we list a few typical plasma parameters that come up often when running simulations.
For each quantity we give both the physical form and its numerical "code counterpart".
Code variables are marked with a hat, :math:`\hat{x}`.


We keep the grid spacing and time step, :math:`\Delta x` and :math:`\Delta t`, explicitly written in the equations even though they are set equal to :math:`1` in the actual numerical calculations;
this allows converting code quantities to physical quantities when needed, by plugging in the spatial and temporal scales.

Plasma skin depth
-----------------

Plasma perturbations that move at the speed of light :math:`c` and oscillate at the plasma/Langmuir frequency
define a length scale known as the plasma skin depth :math:`d_s`.
PIC simulation resolution is usually specified in terms of how many :math:`\Delta x` there are in one :math:`d_s`:

.. math:: d_s = \frac{c}{\omega_{p,s}} = R_s \Delta x

Plasma frequency
----------------

The plasma oscillation frequency, also known as the Langmuir frequency :math:`\omega_{p,s}`,
for a particle species :math:`s` with number density :math:`n_s` and charge :math:`q_s` is

.. math:: \omega_{p,s} = \sqrt{\frac{4\pi n_s q_s^2}{m_s}}

So called total plasma frequency is defined to be:

.. math:: \omega^2_p = \sum_s \omega^2_{p,s}

For relativistic bulk flows with Lorentz factor :math:`\gamma`
or relativistically hot plasmas with mean Lorentz factor :math:`\langle \gamma \rangle \gtrsim 1`
the plasma frequency becomes :math:`\omega_{p,s} \rightarrow \omega_{p,s}/\sqrt{\gamma}`.

Inserting :math:`\hat \omega = \omega \Delta t` to the plasma skin depth
we express the plasma frequency in terms of the resolution:

.. math:: \hat{\omega}_{p,s} = \frac{\hat c}{R_s}

Plasma magnetization
--------------------

A finite magnetic field introduces another set of degrees of freedom to the plasma.
The ratio of the magnetic field line tension, :math:`(B \cdot \nabla) B/4\pi \propto B^2/4\pi`, to the plasma rest mass (relativistic enthalpy density), :math:`\langle \gamma \rangle n m c^2`, is known as the plasma magnetization,

.. math::
    \sigma_s = \frac{B^2}{4\pi n_s \langle \gamma \rangle  m_s c^2}



The magnetization can be used to express the code magnetic field as

.. math::
    \hat{B} = \sqrt{\sigma_s}\hat\omega_{p,s}\hat{c}\langle\gamma\rangle
    \frac{\hat m_s}{|\hat q_s|}
    = \sqrt{\sigma_s \langle \gamma \rangle \hat m_s \hat n_s \hat c^2}

.. admonition:: Derivation
   :collapsible: closed

   .. math::
      \sigma_s = \frac{B^2 q_s^2}
      {\omega_{p,s}^2 \langle\gamma\rangle^2 m_s^2 c^2 }
      = \frac{\Delta t ^2 B_0^2\hat{B}^2 \hat q_s^2 q_0^2}
      {\hat\omega_{p,s}^2 \langle\gamma\rangle^2 \hat{m}_s^2 m_0^2 c^2 } \\

      \therefore \hat{B} = \sqrt{\sigma_s}\hat\omega_{p,s}
      \frac{q_0 \hat c \Delta t}{m_0c}
      \frac {\langle\gamma\rangle \hat{m}_s m_0 c } {\Delta t |\hat q_s| q_0}
      = \sqrt{\sigma_s}\hat\omega_{p,s}
      \frac{\hat{c}\langle\gamma\rangle \hat m_s}{|\hat q_s|} \\

      \omega_{p,s}^2 = \frac{\hat\omega_{p,s}^2}{\Delta t^2}
      = \frac{4\pi n_s q_s^2}{\gamma m_s}
      = \frac{\hat n_s q_s^2}{\Delta x^3 \gamma \hat m_s} \frac{4\pi q_0^2}{m_0}
      = \frac{\hat n_s q_s^2}{\Delta x^3 \gamma \hat m_s} \frac{c^2\Delta x}{\hat c^2} \\

      \therefore \hat B = \sqrt{\sigma_s
      \frac{\Delta t^2 \hat n_s q_s^2}{\Delta x^2 \gamma \hat m_s} \frac{c^2}{\hat c^2}}
      \hat c \gamma \frac{\hat m_s}{|\hat q_s|}
      = \sqrt{\sigma_s \gamma \hat m_s \hat n_s \hat c^2}


.. note::
    The magnetization can also be expressed as a ratio of the magnetic energy density, :math:`U_B = B^2/8\pi` to the plasma rest mass, so that the ratio is :math:`\sigma = B^2/8\pi n m c^2`, i.e., there is a difference of a factor :math:`4\pi` vs :math:`8\pi`.



Gyrofrequency
-------------

The angular frequency of a charged particle gyrating around a magnetic field :math:`B` is known as the gyrofrequency,

.. math::
    \omega_{B,s} = \frac{|q_s| B}{\gamma m_s c}
    \quad\mathrm{and}\quad
    \hat \omega_{B,s} = \sqrt{\sigma_s}\hat \omega_{p,s}

.. admonition:: Derivation
   :collapsible: closed

   .. math::
      \omega_{B,s} = \frac{ \hat \omega_{B,s}}{\Delta t}
      = \frac{|\hat q_s|q_0 \hat B B_0}{\gamma \hat m_s m_0 c}
      = \frac{|\hat q_s|q_0}{\gamma \hat m_s m_0 c}
      \sqrt{\sigma_s}\hat\omega_{p,s}
      \frac{\hat{c}\langle\gamma\rangle \hat m_s}{|\hat q_s|}
      \frac{m_0c}{q_0 \hat c \Delta t} \\

      \therefore \hat \omega_{B,s} = \sqrt \sigma_s\hat \omega_{p,s}



Larmor radius
-------------

A charged particle gyrating in a magnetic field forms a "ring" around the magnetic field line with a radius known as the gyroradius, or Larmor radius,

.. math::
    r_{L,s} = \frac{c}{\omega_{B,s}}
    \quad\mathrm{and}\quad
    \hat r_{L,s} = \frac{\hat c}{\sqrt {\sigma_s}\hat \omega_{p,s}}

.. admonition:: Derivation
   :collapsible: closed

   .. math::
      r_{L,s} = \hat r_{L,s}\Delta x
      = \frac{c\Delta t}{\hat \omega_{B,s}} \\

      \therefore \hat r_{L,s} =
      \frac{c}{\sqrt {\sigma_s}\hat \omega_{p,s}} \frac{\hat c}{c}
      = \frac{\hat c}{\sqrt {\sigma_s}\hat \omega_{p,s}}




Additional ratios of scales
---------------------------

These definitions also allow a slightly different way of expressing the magnetization, as a ratio of the gyrofrequency to the plasma frequency, or as a ratio of the Larmor radius to the skin depth,

.. math::
    \sigma = \left( \frac{\omega_{B,s}}{\omega_{p,s}} \right)^2 = \left( \frac{d_s}{r_{L,s}} \right)^2

Note that a high magnetization means that the gyrofrequency increases, :math:`\omega_B \propto \sqrt{\sigma}` and Larmor radius decreases, :math:`r_L \propto \gamma/\sqrt{\sigma}`;
therefore, we need to be careful that the particle gyrations are still resolved, :math:`\Delta x < r_L`.




