.. default-role:: math

.. _qed-link:

QED Monte Carlo
===============

The ``qed`` module makes the particle species of a tile react through the elementary QED
processes: Compton scattering, pair annihilation, photon–photon pair production and synchrotron
self-absorption between pairs of particles, and quantum synchrotron emission and magnetic
(Breit–Wheeler) pair production of single particles in the tile's magnetic field. It is a
Monte Carlo of weighted particles: every computational particle stands for a number of real ones,
and reactions are drawn so that the *expected* numbers of every process are exact while the
particle count stays adaptive.

Species and weights
-------------------

A QED tile is a tile of the simulation that carries a ``qed::tile_state`` (the layouts of its
species and the sampler) next to its Yee lattice and particle containers.
``runko.qed.threeD.Tile(idx, conf)`` builds one in a context of its own, as
``runko.pic.threeD.Tile`` does, and a ``TileGrid`` takes it like any other tile. Its species play
three roles found from their mass and charge: a massless species (``m = 0``) holds the photons, a species of unit mass the leptons —
electrons by a negative charge, positrons by a positive one — and any heavier species
(``m > 1``: protons, ions) takes no part in the reactions::

    q0: -1
    m0: 1
    w0: 1e-3     # lepton weight, in Thomson depths across the tile
    q1: 1
    m1: 1
    w1: 1e-3
    q2: 0        # photons: massless, a weight per particle
    m2: 0
    n_payload2: 1
    particle_pusher: ["boris", "boris", "ballistic"]

Leptons carry one weight per species, ``w{i}`` (1 if absent), equal for both signs; photons
carry a weight per particle, which the resolution function of the sampler adapts. It is the
first float of the per-particle payload of the species, so the photon species needs
``n_payload{i}: 1`` (a tile without it, or with a payload on a lepton species, is refused).
``tile.get_weights(i)`` returns the weights of the live particles of species ``i`` either way.
Photons deposit no current, and take the ``ballistic`` pusher, which moves them at the speed of
light along their momentum; a neutral species under any other pusher is refused when the pushers
are loaded.

Normalization
-------------

Weight is measured in units of the Thomson depth a particle contributes across the tile,
`W_0 = H^2/\sigma_\mathrm{T}`, so that the pair depth `\tau_\mathrm{T}` is the summed lepton
weight, and time in light crossings of the tile, `\mathrm{d}t = \mathrm{cfl}/N_\mathrm{cells}`.
The tile edge in Compton wavelengths, ``qed_H_lamC``, fixes the one-body clock
`\mathcal{N}_1 = (\lambda_C/H)/(\alpha\,\mathrm{d}t)`; the two-body clock is
`\mathcal{N}_2 = 1/\mathrm{d}t`. The magnetic field enters as the magnetic compactness
``qed_ell_B`` of the box edge and the floor of the synchrotron band as ``qed_x_min``; without a
field the one-body channels and self-absorption are off. In a 0-D box the field is uniform and
perpendicular to every particle. Inside PIC the channels read the tile's own lattice: with
``qed_B0``, the code-unit field strength that ``qed_ell_B`` stands for, the quantum parameter of a
particle is the invariant of the field interpolated at its position,
`\chi = \gamma\, (b_0/B_0)\, |(\mathbf{E} + \boldsymbol{\beta}\times\mathbf{B})^2 -
(\boldsymbol{\beta}\cdot\mathbf{E})^2|^{1/2}` (for a photon `\gamma \to x` and
`\boldsymbol{\beta}` its direction), so that emission and magnetic pair production follow the
local field. Self-absorption, whose envelope needs one field per tile, takes the tile's median
`|B|` along the lepton's pitch angle to the tile's mean field — within 10–20 % of the
local-field opacity in developed turbulence, where an external `B_0` is off by 3–10× because
the absorption of the cool tail is exponential in the field.

Configuration
-------------

On top of the ``pic`` keys (the tile must be cubic):

===========================  ==========================================================  =========
key                          meaning                                                     default
===========================  ==========================================================  =========
``qed_H_lamC``               tile edge in Compton wavelengths                            required
``qed_ell_B``                magnetic compactness `\sigma_\mathrm{T} H U_B/(m_e c^2)`     0
``qed_B0``                   the field of ``qed_ell_B`` in code units; 0: uniform field     0
``qed_x_min``                floor of the synchrotron band                               0
``qed_channels``             list of channels to run, e.g. ``["compton-", "pairann"]``   all
``qed_nbin``                 energy bins per reaction cell and species                   16
``qed_n_reaction_cells``     reaction sub-cells per tile edge                            1
``qed_nsort``                steps between padded sorts                                  10
``qed_slack``                slack slots per key, as a fraction of its live particles    0.5
``qed_smin``                 floor on that count                                         4
``qed_a_res``, ``qed_x_ref`` photon resolution function `(x/x_\mathrm{ref})^a`            0.7, 1e-5
``qed_n_cap``                largest copy number of one reaction                         32
``qed_acc_syn``              coefficient of the emission accumulation, 0 disables        3.0
``qed_kappa_acc``            merged exchange `\le \kappa(\gamma-1)`, 0 disables           0.1
``qed_nsub_max``             sub-steps of an additive channel before rows are clamped    8
``qed_seed``                 seed of the counter-based random streams                    1
===========================  ==========================================================  =========

The channels are ``compton-``, ``compton+``, ``pairann``, ``photann``, ``ssa-``, ``ssa+``
(two-body) and ``syn-``, ``syn+``, ``bw`` (one-body).

The step
--------

Every ``qed_nsort`` steps, and whenever a species has lost its layout, the sampler lays each
species out by a counting sort on its key — the reaction cell, the energy bin and, for photons,
a random bit — with slack behind every key for the products of the coming steps; this layout is
the sort of a QED tile, and the sampler makes it itself (``tile.layout()`` makes it at once).
The channels then run in a random order. A two-body channel proposes, per pair of energy bins, an envelope of
events from the bin tables alone, pairs particles by two-level stratified sampling, thins each
proposal against the true rate with a single uniform, and realizes the accepted events with
per-cell systematic rounding; a one-body channel evaluates the rate of every particle and fires
on it. Products land in the slack of their own key, reactive at once, or in the tail, inert
until the next sort; a parent that a kick took out of its energy bin is likewise inert until
the next sort, so that the bounds of a bin stay those of the bin.

In a lap the driver calls ``x.prtcl_interact()`` for the reactions and, for an open box,
``x.grid_heat_diffusive(ell)``, ``x.grid_inject_pairs_plaw(ell, p, gmin, gmax)``,
``x.grid_inject_photons(ell, theta, w)`` and ``x.grid_escape_photons(tau_T, dt)`` — the slab
escape against the pair depth ``tau_T`` and the step ``dt`` in crossings of the slab, by default
the tile's own. Each is an action (``interact``, ``heat_diffusive``, ...) over every local QED
tile; the arguments are floats. Every tile keeps its budgets, the escaping spectrum and its
angular distribution, the channel counts and a kernel report, which the ``ProxyTile`` of the tile
reads (``tile.budget()``, ``tile.escaped_spectrum()``, ``tile.escaped_mu_spectrum()``,
``tile.counts()``, ``tile.moved()``, ``tile.report()``, ``tile.totals("e-")``, ...). The same
``ProxyTile`` methods (``tile.interact()``, ``tile.escape_photons()``, ``tile.set_particles(s, u, w)``,
...) act on that one tile only; the species are named ``"ph"``, ``"e-"`` and ``"e+"``. The
``ProxyTile`` methods that runko's pic and emf layers provide (``tile.sort_particles()``,
``tile.push_particles()``, ``tile.inject(...)``, ``tile.set_EBJ(...)``, ...) instead evaluate an
action over every local tile of the context the proxy belongs to: the tile alone for a tile of
``runko.qed.threeD.Tile``, but all local tiles for one of ``simulation.local_tiles()``.

Inside PIC
----------

QED tiles take part in an otherwise ordinary PIC lap. The reactions run first, when every particle is inside its tile (the last thing that moved them was
the exchange of the previous lap), so that the products, born at their parents' positions, take
the push and the exchange of the lap like any particle. The sampler lays its species out
itself, so the lap needs no ``prtcl_sort``; calling it is harmless, since a pic sort only drops
the layout and the next ``interact`` lays the species out again::

    x.prtcl_interact()
    x.grid_escape_photons(tau_box, dt_box)    # tau_box: the lepton weight of all ranks / n_t^2
    x.grid_push_half_b(); comm B
    x.prtcl_push(); comm particles
    x.prtcl_deposit_current(); comm and filter J
    x.grid_push_half_b(); comm B
    x.grid_push_e(); x.grid_deposit_antenna_current(); x.grid_add_current(); comm E

The layout survives the particle exchange and an injection through the pic path
(``tile.inject``): both append behind the last alive slot of the container, into the slack of the
keys there or past the layout as tail; slack that the append dropped from the end of the
container is given back to its keys. An arrival in slack stays where it is and is no member
of the key yet: a product claims a slack slot only while its id is dead and otherwise goes to
the tail, the key's bounds taking the arrival in as its fill reaches it. An arrival inside a
key's fill reacts with the key at once (every step widens the bounds by the members of every
key before any
channel runs), one in the slack once the fill reaches it, a tail one from the next sort on.
Photons move at `c` along their momentum,
deposit no current and carry their weight through the exchange as their payload; the particle
snapshots carry the weight as a 13th field ``wgt`` and the spectra are weighted by it (for any
species: payload component 0 when the species has a payload, its ``w{i}`` otherwise). With
`n_t` tiles per edge the weight unit is that of the tile, so the depth of the box is `\sum w / n_t^2` and
``qed_ell_B`` is defined on the box edge with ``qed_H_lamC`` the tile edge.

Example
-------

``projects/qed-corona`` is the hard-state pair corona of Nättilä (2024): heated and injected pairs
in a magnetized slab make every photon themselves and photons escape; ``qed.py --conf corona.ini``
runs it and ``plot_corona.py`` / ``plot_weights.py`` draw the pair and energy budgets, the channel
rates, the escaping and internal spectra, the lepton distribution and the photon weights.
``qed-pic.py --conf corona-pic.ini`` is the same corona inside the PIC loop: a periodic box of
`2^3` tiles over any number of ranks, energized by the Langevin antenna of
``projects/pic-turbulence`` alone, the channels acting on the turbulent local field and photons
escaping against the global depth. ``qed_diagnostics.py`` reduces the depth, energies, counts and
histograms of all ranks to ``qed_<lap>.npz`` and ``qed_timeseries.npz`` next to the field,
particle and spectra snapshots; ``plot_qed_panel.py`` draws the diagnostic grid of a record with
the time series and the antenna-versus-escape energy balance, ``plot_qed_fields.py`` the field
and density slices, ``plot_qed_spectra.py`` the weighted spectra over time.

Tests
-----

``tests/cpp/test_qed*.c++`` pin the physics functions, check their samplers and the weight rules
of ``qed/weights.h`` (what an event moves and how often it fires); ``tests/py/test_qed_*``
are the analytic tests of the sampler: the closed-box equilibrium, Wien relaxation and Kompaneets
moments, inverse-Compton cooling, annihilation rates and lines, the photon–photon rate equation,
first-step exactness against the pair sum, self-absorption, synchrotron and Breit–Wheeler rates
and bookkeeping, pair balance with escape, weight invariance, a short corona, the layout across
arrivals and the channels following a lattice field; ``tests/py-multirank/test_qed_simulation.py``
runs the module over four ranks with photons crossing tiles.
