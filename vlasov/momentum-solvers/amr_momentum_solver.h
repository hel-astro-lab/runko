#pragma once


#include "../amr/mesh.h"
#include "../tile.h"

//#include <cmath> 
#include <Eigen/Dense>
//#include "tile.h"
//#include "../em-fields/tile.h"
//#include "amr/mesh.h"
//#include "amr/numerics.h"
//#include "amr/refiner.h"
//#include "../tools/cppitertools/zip.hpp"
//#include "../tools/cppitertools/enumerate.hpp"
//#include "../tools/signum.h"
//#include "../units.h"
//#include "amr_analyzator.h"


//using namespace Eigen;
//using std::floor;
//using iter::zip;
//using toolbox::sign;
//using units::pi;

namespace vlv {

using namespace Eigen;


/*! General interface for AMR Momentum space solvers
 *
 * Solve is the general function that feeds stuff to the actual solveMesh()
 * that does all the dirty work. Hence, new solvers should only introduce
 * solveMesh() functions and let solve do what ever it does to loop the meshes
 * correctly from the container.
 *
 */
template< typename T, int D, int V >
class MomentumSolver {

  public:
    typedef std::array<T, 3> vec;

    MomentumSolver() {};

    virtual ~MomentumSolver() = default;


    /// Get snapshot current J_i^n+1 from momentum distribution
    void updateFutureCurrent( vlv::Tile<D>& tile, T cfl);

    void solve( vlv::Tile<D>& tile, T step_size = T(1));

    virtual void solveMesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& E,
        vec& B,
        T qm,
        T cfl) = 0;

};

} // end of namespace vlv

