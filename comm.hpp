#ifndef COMM_H
#define COMM_H


#include "boost/mpi.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wmissing-braces"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#pragma clang diagnostic pop

#include "cell.hpp"
#include "common.h"


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
    );


    // Move all cells in grid to rank 0
    void move_all_to_master(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);

    
    // update particles with multi stage message transfer
	void update_ghost_zone_particles(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


	void update_ghost_zone_currents(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


	void update_ghost_zone_yee_currents(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


    // Updating BY
	void update_ghost_zone_B(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);
    //
    // Updating EY
	void update_ghost_zone_E(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


    /* Multi stage load balancing:
        First we balance the load with Zoltan
        check which cells changed
        update these cells between ranks
    */
	void load_balance(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);

};


#endif
