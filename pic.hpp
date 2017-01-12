#ifndef PIC_H
#define PIC_H

// default libraries
#include "algorithm"
#include "iomanip"
#include "iostream"
#include "fstream"
#include "string"

// boost
#include "boost/array.hpp"
#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"

// load balancing
#include "zoltan.h"

// adaptive grid (loaded with pragmas to avoid error spam)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wmissing-braces"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#pragma clang diagnostic pop


// vector calculus
#include <Eigen/Dense>

// own header files
#include "cell.hpp" // Cell class for individual grid cells
#include "comm.hpp" // Domain decomposition communications
#include "inject.hpp" // Particle initializer & injector
#include "mesh.hpp" // Mesh and field related functions
#include "field.hpp" // actual field propagators
#include "parameters.hpp" // simulation variables
#include "particles.hpp" // particle spatial & momentum pushers
#include "io.hpp" // Simulation saver


// name spaces
using namespace std;
using namespace boost;
using namespace boost::mpi;
using namespace dccrg;
using namespace Eigen;


// rank-specific switch for MPI data transfers
int Cell::transfer_mode = Cell::INIT;


// Vector and Matrix calculus type definitions
typedef Array<double, 4, 1> Vectord4;
typedef Array<double, 4, 4> Matrixd44;
typedef Array<double, 3, 2> Matrixd32; // old compatibility type



/* Main simulation loop

    Here we initialize
        MPI
        Zoltan
        grid (based on dccrg)
        
    inject particles
    construct initial fields

    And then roll with simulation loop

*/


int main(int argc, char* argv[]);


#endif

