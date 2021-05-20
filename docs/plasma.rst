Plasma parameters
=================

Here we list a few typical plasma parameters that are often encountered when dealing with the simulations.
We list both the physical parameters and the numerical "code counterparts" for each quantity. 
Code variables are marked with a hat, :math:`\hat{x}`.


We keep the grid spacing and time step, :math:`\Delta x` and, :math:`\Delta t`, explicitly written in the equations even though they are set equal to :math:`1` in the actual numerical calculations;
this allows converting the code quantities to physical quantities, when needed, by plugging in the spatial and temporal scales.




Plasma frequency
----------------

The plasma oscillation frequency, also known as the Langmuir frequency, :math:`\omega_{p,s}` for a particle species :math:`s` with number density :math:`n_s` and charge :math:`q_s` is

.. math::
    \omega_{p,s} = \sqrt{\frac{4\pi n_s q_s^2}{m_s}}
    \quad\mathrm{and}\quad
    \hat{\omega_p} = \frac{\hat{c}}{\hat{\mathcal{R}}} \frac{1}{\Delta t} 

For relativistic bulk flows with Lorentz factor :math:`\Gamma` or relativistically hot plasmas with mean Lorentz factor :math:`\langle \gamma \rangle \gtrsim 1` the plasma frequency becomes :math:`\omega_p \rightarrow \omega_p/\sqrt{\gamma}`.


Simulation laps are typically in units of total plasma frequency so that one time step :math:`\Delta t = \hat{\omega_p}^{-1}` in physical units.


Plasma skin depth
-----------------

Plasma perturbations that move with the speed of light :math:`c` and oscillates with the Langmuir frequency define a length scale known as the plasma skin depth,

.. math::
    d_s = \frac{c}{\omega_{p,s}}
    \quad\mathrm{and}\quad
    \hat{d}_e = \hat{\mathcal{R}} \Delta x

One skin depth in the code is :math:`\hat{\mathcal{R}}` grid cells.


Plasma magnetization
--------------------

A finite magnetic field sets another set of degrees of freedom to the plasma. 
The ratio of the magnetic field line tension, :math:`(B \cdot \nabla) B/4\pi \propto B^2/4\pi` to the plasma rest mass (relativistic enthalphy density), :math:`\langle \gamma \rangle n m c^2` is known as the plasma magnetization,

.. math::
    \sigma_s = \frac{B^2}{4\pi n_s \langle \gamma \rangle  m_s c^2}



The magnetization can be used to express the code magnetic field as

.. math::
    \hat{B} = \frac{ \hat{c}^2 \sqrt{\sigma} }{\hat{\mathcal{R}}} \frac{\Delta x}{\Delta t^2}


.. note::
    The magnetization can also be expressed as a ratio of the magnetic energy density, :math:`U_B = B^2/8\pi` to the plasma rest mass, so that the ratio is :math:`\sigma = B^2/8\pi n m c^2`, i.e., there is a difference of a factor :math:`4\pi` vs :math:`8\pi`.



Gyrofrequency
-------------

Angular frequency of a charged particle gyrating around a magnetic field :math:`B` is known as the gyrofrequency,

.. math::
    \omega_B = \frac{|q_s| B}{\gamma m_s c}
    \quad\mathrm{and}\quad
    \hat{\omega_B} = \frac{\hat{c} \sqrt{\sigma} }{\hat{\mathcal{R}}} \frac{1}{\Delta t}


Larmor radius
-------------

The charged particle gyrating in the magnetic field forms a "ring" around the magnetic field line with a radius known as the gyroradius or a Larmor radius,

.. math::
    r_L = \frac{c}{\omega_B} = \frac{\gamma \beta m_s c^2}{|q_s| B}
    \quad\mathrm{and}\quad
    \hat{r}_L = \frac{\hat{c}}{ \sqrt{\sigma}} \Delta x



Additional ratios of scales
---------------------------

These definitions also allows a slightly different way of expressing the magnetization as a ratio of the gyrofrequency to the plasma frequency, or a ratio of the Larmor radius to the skin depth,

.. math::
    \sigma = \left( \frac{\omega_B}{\omega_p} \right)^2 = \left( \frac{d_e}{r_L} \right)^2

Note that a high magnetization means that the gyrofrequency increases, :math:`\omega_B \propto \sqrt{\sigma}` and Larmor radius decreases, :math:`r_L \propto \gamma/\sqrt{\sigma}`;
therefore, we need to be careful that the particle gyrations are still resolved, :math:`\Delta x < r_L = \sqrt{\sigma}/\hat{\mathcal{R}} \hat{c}`.




