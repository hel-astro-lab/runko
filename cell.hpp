#ifndef CELL_HPP
#define CELL_HPP

#include "vector"
#include "array"
#include "iostream"


#include <Eigen/Dense>
#include "mpi.h"

using namespace Eigen;

/// 4x1 vector FIXME multiple definitions of the same type also in pic.cpp
typedef Array<double, 4, 1> vec4;
// typedef Array<double, 3, 1> vec3;

/// Full 4x4 GR vector
typedef Array<double, 4, 4> mat4;

// Define cell types for different boundary conditions
#define NORMAL_CELL 0               ///< Cell in a normal state
#define METAL_BOUNDARY_X_CELL 1     ///< Cell with a conducting X boundary
#define METAL_BOUNDARY_Y_CELL 2     ///< Cell with a conducting Y boundary



/// Modular Cell class for mpiGrid
/**
Cell class to handle MPI communications
*/

class Cell
{
public:
        
    /// current number of particles inside the Cell
	unsigned int number_of_particles = 0;

	/// coordinates of particles in this cell
	std::vector<std::array<double, 6> > particles;

    /// 4-current vector containing (rho, Jx, Jy, Jz)
    vec4 J;

    /// standard 3-vector for currents & field vectors at Yee lattice staggeration
    /**
        XXX Here we use 3+1 formalism which is more easily
        expressed via E and B, not F tensor
    */
    Vector3d JY, BY, EY;


    /* Yee staggered field tensor 
        note that E and B elements are defined in different locations
        due to the staggering
     mat4 FY = mat4::Zero();
    */


    // Yee lattice/mesh stuff array
    /*
        0 - charge density rho
        DONE 1 2 3    - current vector on nodal points Jx, Jy, Jz
        DONE 4 5 6    - current vector on Yee lattice  JxY, JyY, JzY
        7 8 9    - E field on Yee lattice ExY, EyY, EzY
        10 11 12 - B field on Yee lattice BxY, ByY, BzY

        13 14 15 - E field on nodal lattice Ex, Ey, Ez
        16 17 18 - B field on nodal lattice Bx, By, Bz
    std::array<double, 19> field = {{0.0, 
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0}};
    */


    /// incoming currents from other processes
    /**  3x3x3 cube with four-vector elements = 108 values */
    std::array<double, 108> incoming_currents;



    /** @name Auxiliary functions 
     *  auxiliary functions to help disentangle lattice values
     *
     *  note how we always give normal and const definition
     *  this makes sure the variable pointed to by the returned 
     *  pointer & won't be alterable and that the method does not 
     *  alter the variable pointed to by the given pointer. (phew!)
    */

    /*
    ///@{
    /// overload array operator for object
    double& operator [](const std::size_t i){return this->field[i];}
    const double& operator [](const std::size_t i) const {return this->field[i];}

    /// charge density
    double& rho(){ return this->field[0]; }
    const double& rho() const { return this->field[0]; }

    /// nodal current x
    double& Jx(){ return this->field[1]; }
    const double& Jx() const { return this->field[1]; }

    /// nodal current y
    double& Jy(){ return this->field[2]; }
    const double& Jy() const { return this->field[2]; }

    /// nodal current z
    double& Jz(){ return this->field[3]; }
    const double& Jz() const { return this->field[3]; }

    /// Yee currents x
    double& JxY(){ return this->field[4]; }
    const double& JxY() const { return this->field[4]; }

    /// Yee currents y
    double& JyY(){ return this->field[5]; }
    const double& JyY() const { return this->field[5]; }

    /// Yee currents z
    double& JzY(){ return this->field[6]; }
    const double& JzY() const { return this->field[6]; }

    /// Yee E x
    double& ExY(){ return this->field[7]; }
    const double& ExY() const { return this->field[7]; }

    /// Yee E y
    double& EyY(){ return this->field[8]; }
    const double& EyY() const { return this->field[8]; }

    /// Yee E z
    double& EzY(){ return this->field[9]; }
    const double& EzY() const { return this->field[9]; }

    /// Yee B x
    double& BxY(){ return this->field[10]; }
    const double& BxY() const { return this->field[10]; }

    /// Yee B y
    double& ByY(){ return this->field[11]; }
    const double& ByY() const { return this->field[11]; }

    /// Yee B z
    double& BzY(){ return this->field[12]; }
    const double& BzY() const { return this->field[12]; }

    /// nodal E field vector x
    double& Ex(){ return this->field[13]; }
    const double& Ex() const { return this->field[13]; }

    /// nodal E field vector y
    double& Ey(){ return this->field[14]; }
    const double& Ey() const { return this->field[14]; }

    /// nodal E field vector z
    double& Ez(){ return this->field[15]; }
    const double& Ez() const { return this->field[15]; }


    /// nodal B x
    double& Bx(){ return this->field[16]; }
    const double& Bx() const { return this->field[16]; }

    /// nodal B y
    double& By(){ return this->field[17]; }
    const double& By() const { return this->field[17]; }

    /// nodal B z
    double& Bz(){ return this->field[18]; }
    const double& Bz() const { return this->field[18]; }
    ///@}
    */


    /// Cell type XXX: use me!
    double cell_type = 0;

    /// list of remote neighbors
    std::array<uint64_t, 27> 
        remote_neighbor_list = {{0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0
                                }};

    /// transfer switch 
	/**
     Defines what is transferred over MPI
	*/
	static int transfer_mode;

    /// Define transfer modes
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
        // current four vector
        CURRENT                          = 6,
        // current four vector
        YEE_CURRENT                      = 7,
        // current four vector
        YEE_B                            = 8,
        // current four vector
        YEE_E                            = 9,
		// cell type data
		TYPE                             = 10;


    /// handle the MPI calls depending on cell state
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype();


	/// reserves space for particle data coming over MPI
	void resize()
	{
		this->particles.resize(this->number_of_particles);
	}

};





#endif

