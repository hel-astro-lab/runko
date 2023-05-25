.. default-role:: math

Overview of algorithms 
======================

.. raw:: html

    <style> .blue {color:blue;} </style>
    <style> .red {color:red;} </style>

.. role:: blue

.. role:: blue


.. list-table:: Runko algorithms
    :widths: 20 30 50 100 50
    :header-rows: 1

    * - Module
      - Solver
      - Instance
      - Description
      - Reference
    * - `fld <fld.rst>`_
      - 
      - 
      - Maxwell's field equation module
      -
    * -
      - Propagator 
      -
      - Field solver; update `\mathbf{\hat{E}}` and `\mathbf{\hat{B}}`
      -
    * - 
      - 
      - `\texttt{FDTD2}`
      - 2nd-order FDTD solver
      - [Yee1966]_
    * - 
      - 
      - `\texttt{FDTD4}`
      - 4th-order FDTD solver
      - [Greenwood2004]_
    * - 
      - 
      - `\texttt{FDTDGen}`
      - FDTD solver with free coefficients
      - [Blinne2018]_
    * -
      - Filter 
      -
      - Current smoothing
      -
    * - 
      - 
      - `\texttt{Binomial2}`
      - 3-point Binomial digital filter 
      - [Birdsall1985]_
    * - 
      - 
      - `\texttt{Compensator2}`
      - 3-point digital compensator filter 
      - [Birdsall1985]_

    * - `pic <pic.rst>`_
      - 
      - 
      - Particle-in-cell module
      -

    * -
      - Pusher 
      -
      - Particle `\mathbf{\hat{x}}` and `\mathbf{\hat{u}}` update
      -
    * - 
      - 
      - base class 
      - Velocity verlet propagator
      - [Verlet1967]_
    * - 
      - 
      - `\texttt{BorisPusher}`
      - Relativistic Boris pusher
      - [Boris1970]_
    * - 
      - 
      - `\texttt{VayPusher}`
      - Vay pusher
      - [Vay2008]_
    * - 
      - 
      - `\texttt{HigueraCaryPusher}`
      - Higuera-Cary pusher
      - [HigueraCary2017]_
    * - 
      - 
      - `\texttt{rGCAPusher}`
      - Reduced guiding-center pusher
      - [Bacchini2020]_

    * -
      - Interpolator 
      -
      - **Field interpolation to particle's location**
      -
    * - 
      - 
      - `\texttt{LinearInterpolator}`
      - Linear 1st-order interpolator
      - 
    * - 
      - 
      - `\texttt{QuadraticInterpolator}`
      - Quadratic 2nd-order interpolator
      - 
    * - 
      - 
      - `\texttt{CubicInterpolator}`
      - Cubic 3rd-order interpolator
      - 
    * - 
      - 
      - `\texttt{QuarticInterpolator}`
      - Quartic 4th-order interpolator
      - 

    * -
      - Depositer 
      -
      - **Current deposition**
      -
    * - 
      - 
      - `\texttt{ZigZag}`
      - 1st-order ZigZag scheme
      - [Umeda2003]_
    * - 
      - 
      - `\texttt{ZigZag\_2nd}`
      - 2nd-order ZigZag scheme
      - [Umeda2005]_
    * - 
      - 
      - `\texttt{Esikerpov\_2nd}`
      - 2nd-order Esikerpov scheme
      - [Esikerpov2001]_
    * - 
      - 
      - `\texttt{Esikerpov\_4th}`
      - 4th-order Esikerpov scheme
      - [Esikerpov2001]_




-------

References
----------

.. [Yee1966] Yee (1966)

.. [Greenwood2004] Greenwood, Cartwright, Luginsland, Baca (2004)

.. [Blinne2018] Blinne, Schinkel, Kuschel, et al. (2018)

.. [Birdsall1985] Birdsall & Langdon (1985)
    Plasma physics via computer simulations

.. [Verlet1967] Verlet (1967)

.. [Boris1970] Boris (1970)

.. [Vay2008] Vay (2008)

.. [HigueraCary2017] Higuera & Cary (2017)

.. [Bacchini2020] Bacchini, Ripperda, Philippov, Parfrey (2020)

.. [Umeda2003] Umeda, Omura, Tominaga, Matsumoto (2003)

.. [Umeda2005] Umeda, Omura (2005)

.. [Esikerpov2001] Esirkepov (2001)




