#include "field.hpp"

typedef Parameters P;


void Field_Cube::distribute( 
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


void Field_Solver::push_half_B(
        dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    // get local
    auto all_cells = grid.get_cells();
    for (const uint64_t cell: all_cells) {
        std::array<double, 3> dxyz = grid.geometry.get_length(cell);
        Vector3d ds; 
        ds << 1.0/dxyz[0], 1.0/dxyz[1], 1.0/dxyz[2];

        // vector containing shifted +1/+1/+1 elements
        Vector3d EYp1 = Vector3d::Zero();

        auto* const cell_data = grid[cell];

        // collect +1 E elements
        const auto* const neighbors = 
            grid.get_neighbors_of(cell, CP1_SHIFT);

        int dir = 0; 
        for (const auto neigh: *neighbors) {

            if (neigh == dccrg::error_cell){
                // TODO deal with boundary conditions
                continue;
            } else if(neigh == 0) {
                // TODO implement boundary condition switch
                // now we assume metals because EYp1 has zeros by default
                continue;
            } else {
                auto* const neigh_data = grid[neigh];
                EYp1(dir) = neigh_data->EY(dir); 
            }
            dir++;
        }  // end of neighbors

        /*
Theory:
E = En - E
d = c*dt/2*(dx, dy, dz)
B = BY_n-1

curl E = nabla x E
=+i(dy*Ez - dz*Ey)
-j(dx*Ez - dz*Ex)
+k(dx*Ey - dy*Ex)

B_n+1 = B_n1 -curl E
*/

        // vec3 dE = cell_data->EY - EYp1;
        cell_data->BY -= c*P::dt/(2.0) * ds.cross(cell_data->EY - EYp1); 

    }

    return; 
}


void Field_Solver::push_E(
        dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    // get local
    auto all_cells = grid.get_cells();
    for (const uint64_t cell: all_cells) {
        std::array<double, 3> dxyz = grid.geometry.get_length(cell);
        Vector3d ds;
        ds << 1.0/dxyz[0], 1.0/dxyz[1], 1.0/dxyz[2];

        // vector containing shifted -1/-1/-1 elements
        Vector3d BYm1 = Vector3d::Zero();

        auto* const cell_data = grid[cell];

        // collect -1 B elements
        const auto* const neighbors = 
            grid.get_neighbors_of(cell, CM1_SHIFT);

        int dir = 0; 
        for (const auto neigh: *neighbors) {

            if (neigh == dccrg::error_cell){
                // TODO deal with boundary conditions
                continue;
            } else if(neigh == 0) {
                // TODO implement boundary condition switch
                // now we assume metals because BYm1 has zeros by default
                continue;
            } else {
                auto* const neigh_data = grid[neigh];
                BYm1(dir) = neigh_data->BY(dir); 
            }
            dir++;
        }  // end of neighbors


        /*
           E_n+1 = E_n + dt*[ curl B - 4pi J ]
           */

        cell_data->EY += c*P::dt * ds.cross(cell_data->BY - BYm1)
            -4.0*pi*P::dt*cell_data->JY;

    }

    return; 
}


