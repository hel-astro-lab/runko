/*
Cell class to handle MPI communications
*/


#ifndef CELL_HPP
#define CELL_HPP

#include "vector"
#include "iostream"

#include <Eigen/Dense>
#include "mpi.h"

using namespace Eigen;

// FIXME multiple definitions of the same type also in pic.cpp
typedef Array<double, 4, 1> vec4;


class Cell
{
public:

	unsigned int number_of_particles = 0;

	// coordinates of particles in this cell
	std::vector<std::array<double, 6> > particles;

    // 4-current vector containing (rho, Jx, Jy, Jz)
    vec4 J;


    /* Yee lattice/mesh stuff
        0 - charge density rho
        1 2 3    - current vector on nodal points Jx, Jy, Jz
        4 5 6    - current vector on Yee lattice  JxY, JyY, JzY
        7 8 9    - E field on Yee lattice ExY, EyY, EzY
        10 11 12 - B field on Yee lattice BxY, ByY, BzY

        13 14 15 - E field on nodal lattice Ex, Ey, Ez
        16 17 18 - B field on nodal lattice Bx, By, Bz

    */
    std::array<double, 19> field = {{0.0, 
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0}};


    // incoming currents from other processes
    //  3x3x3 cube with four-vector elements = 108 values
    std::array<double, 108> incoming_currents;




    // auxiliary functions to help disentangle lattice values

    // note how we always give normal and const definition
    // this makes sure the variable pointed to by the returned 
    // pointer & won't be alterable and that the method does not 
    // alter the variable pointed to by the given pointer. (phew!)

    double& operator [](const std::size_t i){return this->field[i];}
    const double& operator [](const std::size_t i) const {return this->field[i];}

    // charge density
    double& rho(){ return this->field[0]; }
    const double& rho() const { return this->field[0]; }

    // nodal current 
    double& Jx(){ return this->field[1]; }
    const double& Jx() const { return this->field[1]; }

    double& Jy(){ return this->field[2]; }
    const double& Jy() const { return this->field[2]; }

    double& Jz(){ return this->field[3]; }
    const double& Jz() const { return this->field[3]; }

    // Yee currents
    double& JxY(){ return this->field[4]; }
    const double& JxY() const { return this->field[4]; }

    double& JyY(){ return this->field[5]; }
    const double& JyY() const { return this->field[5]; }

    double& JzY(){ return this->field[6]; }
    const double& JzY() const { return this->field[6]; }

    // Yee E
    double& ExY(){ return this->field[7]; }
    const double& ExY() const { return this->field[7]; }

    double& EyY(){ return this->field[8]; }
    const double& EyY() const { return this->field[8]; }

    double& EzY(){ return this->field[9]; }
    const double& EzY() const { return this->field[9]; }

    // Yee B
    double& BxY(){ return this->field[10]; }
    const double& BxY() const { return this->field[10]; }

    double& ByY(){ return this->field[11]; }
    const double& ByY() const { return this->field[11]; }

    double& BzY(){ return this->field[12]; }
    const double& BzY() const { return this->field[12]; }


    // nodal E field vector
    double& Ex(){ return this->field[13]; }
    const double& Ex() const { return this->field[13]; }

    double& Ey(){ return this->field[14]; }
    const double& Ey() const { return this->field[14]; }

    double& Ez(){ return this->field[15]; }
    const double& Ez() const { return this->field[15]; }


    // nodal B
    double& Bx(){ return this->field[16]; }
    const double& Bx() const { return this->field[16]; }

    double& By(){ return this->field[17]; }
    const double& By() const { return this->field[17]; }

    double& Bz(){ return this->field[18]; }
    const double& Bz() const { return this->field[18]; }

    // XXX: use me!
    double cell_type = 0;

    // list of remote neighbors
    std::array<uint64_t, 27> 
        remote_neighbor_list = {{0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0
                                }};

	/*
	The number of particles is transferred over MPI if
	this is false, otherwise the particle coordinates
	are transferred.
	*/
	// static bool transfer_particles;
	static int transfer_mode;
	static const int
		// data related to initialization
		INIT                             = 0,
		// number of particles to expect
		PARTICLES                        = 1,
		// number of particles to expect
		PARTICLE_NUMBER                  = 2,
		// data of fields
		FIELDS                           = 3,
		// neighborhood data
		REMOTE_NEIGHBOR_LIST             = 4,
		// data of incoming currents
		INCOMING_CURRENTS                = 5,
		// cell type data
		TYPE                             = 6;


	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		void* address = NULL;
		int count = -1;
		MPI_Datatype datatype = MPI_DATATYPE_NULL;

		switch(Cell::transfer_mode) {
		case Cell::PARTICLES:
			if (this->particles.size() > 0) {
				address = &(this->particles[0]);
			} else {
				// return a sane address just in case
				address = &(this->number_of_particles);
			}

			count = this->particles.size() * 6;
			datatype = MPI_DOUBLE;
			break;
		case Cell::PARTICLE_NUMBER:
			address = &(this->number_of_particles);
			count = 1;
			datatype = MPI_UNSIGNED;
			break;
		case Cell::FIELDS:
			address = &(this->field[0]);
			count = 19;
			datatype = MPI_DOUBLE;
			break;
		case Cell::REMOTE_NEIGHBOR_LIST:
			address = &(this->remote_neighbor_list[0]);
			count = 27;
			datatype = MPI_UINT64_T;
			break;
		case Cell::INCOMING_CURRENTS:
			// address = this->incoming_currents.data();
			address = &(this->incoming_currents[0]);
			count = 108;
			datatype = MPI_DOUBLE;
			break;
		case Cell::TYPE:
			address = &(this->cell_type);
			count = 1;
			datatype = MPI_DOUBLE;
			break;
		default:
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid transfer switch: " << Cell::transfer_mode
				<< std::endl;
			abort();
			break;
		}

		return std::make_tuple(address, count, datatype);
	}


	// reserves space for particle data coming over MPI.
	void resize()
	{
		this->particles.resize(this->number_of_particles);
	}



};

/* Define cell types for different boundary conditions
*/
#define NORMAL_CELL 0 
#define METAL_BOUNDARY_X_CELL 1
#define METAL_BOUNDARY_Y_CELL 2




#endif

