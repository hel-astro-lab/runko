#ifndef COMM
#define COMM


#include "boost/mpi.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "cell.hpp"


using namespace std;
using namespace boost::mpi;
using namespace dccrg;


class Comm
{
public:

    // MPI parameters
    int 
        rank = 0,
        comm_size = 0;

    MPI_Comm 
        communicator;

    bool all_in_master = false;


    // create communicator with default parameters
    // this is the default constructor
    Comm(const int given_rank,
         const int given_comm_size,
         const MPI_Comm given_communicator
    ) {
            this->rank = given_rank;
            this->comm_size = given_comm_size;
            this->communicator = given_communicator;
    };


    // Move all cells in grid to rank 0
	template<
		class CellData,
		class Geometry
	> void move_all_to_master(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {
        
        this->all_in_master = true;

        auto cells = grid.get_cells();
        for (const auto& cell: cells) {
            grid.pin(cell, 0);
        }
        // TODO: set proper INIT mpi_datatype that sends nothing
        Cell::transfer_mode = Cell::PARTICLE_NUMBER;

        grid.balance_load(false); // send everything to rank 0; do not use Zoltan

      };

    // update particles with multi stage message transfer
	template<
		class CellData,
		class Geometry
	> void update_ghost_zone_particles(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {

        // reserve space for incoming particles in copies of remote neighbors
        Cell::transfer_mode = Cell::PARTICLE_NUMBER;
        grid.update_copies_of_remote_neighbors();

        const std::vector<uint64_t>& remote_neighbors
            = grid.get_remote_cells_on_process_boundary();

        for (const auto& remote_neighbor: remote_neighbors) {
            auto* const data = grid[remote_neighbor];
            data->resize();
        }

        // update particle data between neighboring cells on different processes
        Cell::transfer_mode = Cell::PARTICLES;
        grid.update_copies_of_remote_neighbors();


        return;
    };

	template<
		class CellData,
		class Geometry
	> void update_ghost_zone_currents(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {
        Cell::transfer_mode = Cell::REMOTE_NEIGHBOR_LIST;
        grid.update_copies_of_remote_neighbors();

        Cell::transfer_mode = Cell::INCOMING_CURRENTS;
        grid.update_copies_of_remote_neighbors();

        const std::vector<uint64_t> remote_neighbors
            = grid.get_remote_cells_on_process_boundary();

        for (const auto& remote_neighbor: remote_neighbors) {

            auto* const cell_data = grid[remote_neighbor];

            int ijk = 0; 
            for (uint64_t receive_neigh: cell_data->remote_neighbor_list) {
                if(receive_neigh != 0) {

                    #ifdef DEBUG
                    cout << rank << " UC: " << remote_neighbor 
                    << " => " << receive_neigh
                    << " : " << cell_data->incoming_currents[ijk*4] 
                    << " / " << cell_data->incoming_currents[ijk*4+1] 
                    << " / " << cell_data->incoming_currents[ijk*4+2] 
                    << " / " << cell_data->incoming_currents[ijk*4+3] 
                    << endl;
                    #endif


                    // disentangle from std::array to vec4
                    vec4 dJ;
                    dJ << cell_data->incoming_currents[ijk*4],
                          cell_data->incoming_currents[ijk*4+1],
                          cell_data->incoming_currents[ijk*4+2],
                          cell_data->incoming_currents[ijk*4+3];

                    auto* const receive_cell_data = grid[receive_neigh];
                    receive_cell_data->J += dJ;
                }

                ijk++;
            }
        }

        // finally update current vectors
        Cell::transfer_mode = Cell::CURRENT;
        grid.update_copies_of_remote_neighbors();


        return;
    };


    // Updating JY
	template<
		class CellData,
		class Geometry
	> void update_ghost_zone_yee_currents(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {
        Cell::transfer_mode = Cell::YEE_CURRENT;
        grid.update_copies_of_remote_neighbors();

        return;
    };

    // Updating BY
	template<
		class CellData,
		class Geometry
	> void update_ghost_zone_B(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {
        // FIXME update correct neighborhood
        Cell::transfer_mode = Cell::YEE_B;;
        grid.update_copies_of_remote_neighbors();

        return;
    };

    // Updating EY
	template<
		class CellData,
		class Geometry
	> void update_ghost_zone_E(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {

        // FIXME update correct neighborhood
        Cell::transfer_mode = Cell::YEE_E;
        grid.update_copies_of_remote_neighbors();

        return;
    };



    /* Multi stage load balancing:
        First we balance the load with Zoltan
        check which cells changed
        update these cells between ranks
    */
	template<
		class CellData,
		class Geometry
	> void load_balance(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {

        // if everything was forcefully pinned to master we unpin each cell
        if(this->all_in_master){
            grid.unpin_all_cells();
        }


        // start balance with Zoltan
        grid.initialize_balance_load(true);


        // cells that will be moved to this process
        const std::unordered_set<uint64_t>& added_cells
            = grid.get_cells_added_by_balance_load();
        // cells that will be moved from this process
        const std::unordered_set<uint64_t>& removed_cells
            = grid.get_cells_removed_by_balance_load();

        vector<uint64_t> all_transferred_cells;
        all_transferred_cells.insert(all_transferred_cells.end(), 
                added_cells.begin(),
                added_cells.end());
        all_transferred_cells.insert(all_transferred_cells.end(), 
                removed_cells.begin(),
                removed_cells.end());

        // transfer in parts
        Cell::transfer_mode = Cell::PARTICLE_NUMBER;
        grid.continue_balance_load();

        // resize according to incoming number of particles
        for (const uint64_t cell: all_transferred_cells) {
            // cout << rank << ": I got/lost cell " << cell << endl;

            auto* const cell_data = grid[cell];
            if (cell_data == NULL) {
                std::cerr << __FILE__ << ":" << __LINE__
                    << " No data for cell " << cell
                    << endl;
                abort();
            }
            cell_data->resize();
        }

        // Transfer actual particles
        Cell::transfer_mode = Cell::PARTICLES;
        grid.continue_balance_load();

        grid.finish_balance_load();

        // Finally update particles in the ghost zone
        // TODO: do we really need this or does load_balance do it?
        // TODO: update rest of the cell data too; general update function?
        this->update_ghost_zone_particles(grid);


        return;
    };


};



#endif
