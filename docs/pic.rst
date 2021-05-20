Particle-in-cell 
================

Runko has a Particle-in-Cell (PIC) code implemented as a physics module.


Solving coupled Maxwell-Vlasov equations.
Relativistic collisionless plasma.

Instead of modeling the 6D phase space distribution the spatial and momentum dimensions are sampled with particles.
The computational particles interact with electromagnetic fields that are defined on a grid.

Position and velocities of the particles is updated according to the Lorentz force.
The electromagnetic fields are advanced in time accordiing to Ampere's and Faray's laws.



Initialization
--------------

e^- and e^+ are initialized on top of each other.
Initial charge density is then rho = 0.
No need to solve the nabla E = 4pi rho.


Cold mobile ions can be implemented in to the simulation by having rho \ne 0;
in practise this corresponds to having a ghost charges (that do not evolve) on the grid.


Time advancement
----------------

Advance B half time step forward according to Faray's law.

Compute Lorentz force on particles
Need to interpolate E and B fields from the grid to the location of the particles.
This is done with the interpolator.
Update particle velocity and location with a pusher.

Advance the second half of B field according to the Faray's law.

Advance a full time step of the electric field according to the Ampere's law.

Deposit particle current to the grid.
Done using the deposit solver.


Optionally smooth the resulting currents.
Done using the filter solver.

Add currents to the electric field.


