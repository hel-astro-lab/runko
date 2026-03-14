class MovingInjector:
    """Moving particle injector for shock simulations.

    Tracks the position of an injection front that moves in the +x direction.
    Every n_inj laps, injects particles into the stripe between the previous
    and current injector positions, accounting for upstream plasma drift.
    """

    def __init__(self, injloc, beta_inj, beta_flow, cfl, n_inj=1,
                 walloc=0.0, Lx_margin=10.0, Lx=None):
        self.injloc    = injloc
        self.beta_inj  = beta_inj
        self.beta_flow = beta_flow
        self.cfl       = cfl
        self.n_inj     = n_inj
        self.walloc    = walloc
        self.Lx_margin = Lx_margin
        self.Lx        = Lx
        self.moving    = True

    def inject(self, simulation, pgens, ppc):
        """Call every lap. Injects on every n_inj-th lap if still moving.

        Parameters
        ----------
        simulation : runko.Simulation
            The simulation object (provides lap counter and tile access).
        pgens : list of (int, callable)
            Pairs of (species_id, pgen_callable).
        ppc : int
            Particles per cell per species (number of injection passes).
        """
        if not self.moving:
            return
        if simulation.lap % self.n_inj != 0 or simulation.lap == 0:
            return

        stride  = self.n_inj * self.cfl
        x_left  = max(self.injloc - self.beta_flow * stride, self.walloc)
        x_right = self.injloc + self.beta_inj * stride

        if self.Lx is not None and x_right >= self.Lx - self.Lx_margin:
            self.moving = False
            return

        for tile in simulation.local_tiles():
            for _ in range(ppc):
                for species, pgen in pgens:
                    tile.batch_inject_in_x_stripe(species, pgen, x_left, x_right)

        self.injloc = x_right
