MPI Message Size (PIC)
======================


In Runko's PIC module, the functions `pack_all_particles()` and `pack_outgoing_particles()` (defined in `pic/particle.c++`) are responsible for packaging particles leaving an old Tile to be passed to a new Tile via an MPI Message. If there are lots of particles leaving the tile, then the performance can take a hit from constant resizing of the arrays. In that case, you can try and increase the `first_message_size` in `particle.h`. An optimal simulation would always send majority of the particles via the first message and whatever remains (fluctuations) via the extra message.


.. note::

   Note: the MPI_ERR_TRUNCATE bug described here should be outdated. It is left here in case it resurfaces.


If too many particles are leaving a Tile, then this can cause an `MPI_ERR_TRUNCATE` (error code 15). This usually occurs if the simulation conditions are very extreme with high particle numbers, and thus exceed the maxiumum that Runko can handle (equal to `first_message_size+extra_message_size` in `particle.h`, by default this is 4096+4096 particles). You can estimate the maximum number of particles leaving a cell with

`#outgoing_particles = #number_of_cells_in_tile_border * cfl * #species * parrticles_per_cell * factor`,

where `factor` is a fudge factor which quantifies the maximum expected overdensity in your simulation. If your estimate for `#outgoing_particles` exceeds the default value of `first_message_size` in `pic/particle.h`, then this may cause MPI errors. To fix this you must increase the value of `extra_message_size` in the same file in order to prevent any errors.

