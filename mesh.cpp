#include "mesh.hpp"


void Neighborhood_Cube::distribute(
    dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid, 
    uint64_t cell
) {

    if (grid.is_local(cell)) {

        // First distribute into the cell itself
        auto* const center_cell_data = grid[cell];
        if (center_cell_data == NULL) {
            cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
            abort();
        }

        center_cell_data->J += this->J(0,0,0);

#ifdef DEBUG
        cout << " DIS: 000 " << cell
            << " : " << center_cell_data->J(0) << " | " << fields[13](0)
            << " / " << center_cell_data->J(1) << " | " << fields[13](1)
            << " / " << center_cell_data->J(2) << " | " << fields[13](2)
            << " / " << center_cell_data->J(3) << " | " << fields[13](3)
            << endl;
#endif


        // next distribute to neighbors
        // TODO work is wasted going from -1 -> 0 indices; optimize

        const auto* const neighbors = grid.get_neighbors_of(cell);
        int ijk = 0;
        for (const auto neigh: *neighbors) {

            if(ijk == 13){ijk++;}; // skip (0,0,0)

            if (neigh == dccrg::error_cell){
                // TODO deal with boundary conditions
            } else if (grid.is_local(neigh)) {

                auto* const cell_data = grid[neigh];
                if (cell_data == NULL) {
                    cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
                    abort();
                }

                cell_data->J += fields[ijk];

#ifdef DEBUG
                cout << " DIS: " << cell
                    << " => loc " << neigh << " ijk: " << ijk 
                    << "   J= " << center_cell_data->J(0)
                    << " / "    << center_cell_data->J(1) 
                    << " / "    << center_cell_data->J(2) 
                    << " / "    << center_cell_data->J(3) 
                    << endl;
#endif

                // else the neighbor is not local and we have to send it
            } else {
#ifdef DEBUG
                cout << "DIS: " << cell <<  " => rem " << neigh << " : "
                    << fields[ijk](0) << " / " << fields[ijk](1) << " / "
                    << fields[ijk](2) << " / " << fields[ijk](2) << endl;
#endif

                center_cell_data->remote_neighbor_list[ijk] = neigh;

                // parse vec4 into std::array
                center_cell_data->incoming_currents[ijk*4+0] += fields[ijk](0);
                center_cell_data->incoming_currents[ijk*4+1] += fields[ijk](1);
                center_cell_data->incoming_currents[ijk*4+2] += fields[ijk](2);
                center_cell_data->incoming_currents[ijk*4+3] += fields[ijk](3);

            }
            ijk++;
        }

        // centering cell is not local but remote
    } else { 
        cerr << __FILE__ << ":" << __LINE__ << " Non local cell " << cell << endl;
    }


    return;
}



void Mesh::sort_particles_into_cells(
    dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    auto cells = grid.get_cells();
    const std::vector<uint64_t>& remote_neighbors
        = grid.get_remote_cells_on_process_boundary();
    cells.insert(cells.begin(), remote_neighbors.begin(), remote_neighbors.end());

    
    // move particles to the particle list of the cell the particles are currently inside of
    for (const auto& previous_cell: cells) {

        auto* const previous_data = grid[previous_cell];
        if (previous_data == NULL) {
            cerr << __FILE__ << ":" << __LINE__
                << " No data for cell " << previous_cell
                << endl;
            abort();
        }

        std::array<double, 3> grid_coords;
        std::array<double, 6> particle_coords_vel;
        vector<std::array<double, 6>>::size_type i = 0;

        for (int ipops=Population::ELECTRONS; ipops < Population::N_POPULATIONS; ipops++)
        {
            Population ptype = (Population) ipops; ///< type cast into Population type for later usage

            while (i < previous_data->particles(ptype).size()) {

                // handle grid wrap around
                // FIXME beautiful (not working) oneliner
                // previous_data->particles[i] = grid.geometry.get_real_coordinate(previous_data->electrons[i]);

                particle_coords_vel = previous_data->particles(ptype)[i];
                std::array<double, 3> particle_coords = {{particle_coords_vel[0],
                    particle_coords_vel[1],
                    particle_coords_vel[2]}};
                grid_coords = grid.geometry.get_real_coordinate(particle_coords);

                particle_coords_vel[0] = grid_coords[0];
                particle_coords_vel[1] = grid_coords[1];
                particle_coords_vel[2] = grid_coords[2];

                previous_data->particles(ptype)[i] = particle_coords_vel;

                // const uint64_t current_cell = grid.get_existing_cell(previous_data->particles[i]);
                const uint64_t current_cell = grid.get_existing_cell(particle_coords);

                // do nothing if particle hasn't changed cell
                if (current_cell == previous_cell) {
                    i++;
                    continue;
                }

                // change only local cell data because remote data will be overwritten anyway
                // add particle to the current cell's list
                if (grid.is_local(current_cell)) {
                    auto* const current_data = grid[current_cell];
                    // XXX Bug here? with wrap coordinates???
                    current_data->particles(ptype).push_back(previous_data->particles(ptype)[i]);
                    current_data->number_of(ptype) = current_data->particles(ptype).size();
                }

                // remove particle from its previous cell
                if (grid.is_local(previous_cell)) {
                    previous_data->particles(ptype).erase(previous_data->particles(ptype).begin() + i);
                    previous_data->number_of(ptype) = previous_data->particles(ptype).size();
                } else {
                    i++;
                }
            }
        }
    }


    return;
}


void Mesh::deposit_currents(
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    // get local and neighboring cells
    auto all_cells = grid.get_cells();
    for (const uint64_t cell: all_cells) {
#ifdef DEBUG
        cout << this->rank << ": mesh.deposit: my all cells are:"
            << cell << endl;
#endif

        auto* const cell_data = grid[cell];
        if (cell_data == NULL) {
            cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
            abort();
        }
        // wipe out previous charge and currents from mesh
        // this is done to ensure proper summation
        cell_data->J = vec4::Zero();

        // cout << this->rank << " mesh.deposit: J=" << cell_data->J << endl;
    }


    const std::vector<uint64_t>& remote_neighbors
        = grid.get_remote_cells_on_process_boundary();
#ifdef DEBUG
    for (const uint64_t& cell: remote_neighbors) {
        cout << this->rank << ": mesh.deposit: my remote cells on boundary are:"
            << cell << endl;

    }
#endif

    const std::vector<uint64_t>& local_cells_on_boundary
        = grid.get_local_cells_on_process_boundary();
    for (const uint64_t& cell: local_cells_on_boundary) {
#ifdef DEBUG
        cout << this->rank << ": mesh.deposit: my local cells on boundary are:"
            << cell << endl;
#endif

        auto* const cell_data = grid[cell];
        if (cell_data == NULL) {
            cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
            abort();
        }

        for (int ijk=0; ijk<=108; ijk++) {
            cell_data->incoming_currents[ijk] = 0.0;
        }

        // undo remote neighbors
        cell_data->remote_neighbor_list = {{0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0
        }};

    }

    const std::vector<uint64_t>& local_cells_not_boundary
        = grid.get_local_cells_not_on_process_boundary();
#ifdef DEBUG
    for (const uint64_t& cell: local_cells_not_boundary) {
        cout << this->rank << ": mesh.deposit: my local cells NOT on boundary are:"
            << cell << endl;
    }
#endif


    double q = 0.0;
    double m = 0.0;


    auto local_cells = grid.get_cells();
    for (const uint64_t& cell: local_cells) {

        auto* const cell_data = grid[cell];
        if (cell_data == NULL) {
            cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
            abort();
        }

        // neighborhood cube
        Neighborhood_Cube cube;
        cube.init_currents();


        // grid geometry
        std::array<double, 3> start_coords = grid.geometry.get_min(cell);
        std::array<double, 3> ds = grid.geometry.get_length(cell);

        // charge 
        double rhop = e/(ds[0]*ds[1]*ds[2]);

        for (int ipops=Population::ELECTRONS; ipops < Population::N_POPULATIONS; ipops++)
        {
            Population ptype = (Population) ipops; ///< type cast into Population type for later usage

            switch(ptype)
            {
                case Population::ELECTRONS:
                    q = 1.0; // note that minus is already in the formulas
                    m = me;
                    break;
                case Population::POSITRONS:
                    q = -1.0;
                    m = me;
                    break;
                default:
                    std::cerr << __FILE__ << ":" << __LINE__
                        << " Invalid population switch: " << ptype
                        << std::endl;
                    abort();
                    break;
            }

            for (auto& particle: cell_data->particles(ptype)) {

                double 
                    x = particle[0],
                    y = particle[1],
                    z = particle[2],
                    ux = particle[3],
                    uy = particle[4],
                    uz = particle[5];

                    // Lorentz transform
                    double gamma = sqrt(1.0 + ux*ux + uy*uy + uz*uz);
                    ux *= c/gamma;
                    uy *= c/gamma;
                    uz *= c/gamma;

                    double 
                        fp = (x-start_coords[0])/ds[0],
                        fq = (y-start_coords[1])/ds[1],
                        fr = (z-start_coords[2])/ds[2];

                        // basis for the current vector
                        vec4 dJ0, dJ;
                        dJ0 << 1.0, 0.5*ux, 0.5*uy, 0.5*uz;
                        dJ0 *= q*rhop;

#ifdef DEBUG
                        cout << rank << " MD: " << cell 
                            << " dJ0= " << dJ0(0) << " / " << dJ0(1)
                            << dJ0(2) << " / " << dJ0(3)
                            << " ||| " << rhop << " / " << fp 
                            << " / " << fq << " / " << fr << endl;
#endif

                        /* now add weight functions:
                           1st order cloud-in-the-cell model 

                           TODO think about AMR; is it done in Cube?
                           TODO add directly to Yee lattice
                           */
                        double wx, wy, wz;
                        for(int xdir=0; xdir<2; xdir++) {
                            wx = xdir == 0 ? 1.0-fp : fp;

                            for(int ydir=0; ydir<2; ydir++) {
                                wy = ydir == 0 ? 1.0-fq : fq;

                                for(int zdir=0; zdir<2; zdir++) {
                                    wz = zdir == 0 ? 1.0-fr : fr;

                                    dJ = dJ0*wx*wy*wz;
                                    cube.add_J(xdir, ydir, zdir, dJ);
                                }
                            }
                        }

            } // end of particle loop
        } // end of population loop


        // now deposit neighborhood cube into mesh
        // cout << this->rank << ": mesh.deposit: cube.distribute:" << endl;
        cube.distribute(grid, cell);
    }
}


void Mesh::yee_currents(
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    // get local cells
    auto cells = grid.get_cells();
    for (const uint64_t& cell: cells) {

        // read Cell data
        auto* const cell_data = grid[cell];

        // into this we build the current in the other edge of the cell
        // XXX vec4 Jp1 = vec4::Zero();
        Vector3d Jp1 = Vector3d::Zero();

        // this returns the face neighbors 100,010,001
        const auto* const neighbors = 
            grid.get_neighbors_of(cell, CP1_SHIFT);


        // XXX int dir = 1; // skip rho and start from Jx
        int dir = 0; 
        for (const auto neigh: *neighbors) {

            if (neigh == dccrg::error_cell){

                // TODO deal with boundary conditions
                continue;
            } else if(neigh == 0) {
                // TODO implement boundary condition switch
                // now we assume metals because Jp1 has zeros by default
                continue;
            } else {
                auto* const neigh_data = grid[neigh];
                Jp1(dir) = neigh_data->J(dir); 
            }
            dir++;
        }  // end of neighbors


        // now shift currents to cell faces
        // FIXME avoid this array->vector copy
        Vector3d J0(cell_data->J.tail(3).data());

        cell_data->JY = (J0 + Jp1)/2.0;

        // XXX copy nodal rho also to Yee vector
        // cell_data->JY(0) = cell_data->J(0);


    }

    return;
}


void Mesh::nodal_fields(
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    // get local cells
    auto cells = grid.get_cells();
    for (const uint64_t& cell: cells) {
        // XXX
        /*
           const auto* const neighbors = grid.get_neighbors_of(cell, Yee_shift);
           if (neigh == dccrg::error_cell){
        // TODO deal with boundary conditions 
        } else {
        auto* const cell_data = grid[neigh];
        }
        */
    }

    return;
}




