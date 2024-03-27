#pragma once

#include "vlv/amr/mesh.h"
#include "vlv/tile.h"
#include "tools/vector.h"

//#include <Eigen/Dense>


namespace vlv {
//using namespace Eigen;


/// small container for all the solver parameters
namespace tools {
  template<typename T>
  struct Params 
  {
    T qm;
    T cfl;
    T xloc;
  };
}


/*! General interface for AMR Momentum space solvers
 *
 * Solve is the general function that feeds stuff to the actual solve_mesh()
 * that does all the dirty work. Hence, new solvers should only introduce
 * solve_mesh() functions and let solve do what ever it does to loop the meshes
 * correctly from the container.
 *
 */
template< typename T, int D, int V >
class MomentumSolver {

  public:
    using vec = std::array<T, 3>;
    //using Vec3E = Vector3d;
    using Vec3E = toolbox::Vec3<T>; // replace with own vector type

    MomentumSolver() = default;

    virtual ~MomentumSolver() = default;


    /// Get snapshot current J_i^n+1 from momentum distribution
    void update_future_current( vlv::Tile<D>& tile, T cfl);

    void solve( vlv::Tile<D>& tile, T step_size = T(1));

    virtual void solve_mesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& E,
        vec& B,
        vlv::tools::Params<T>& params
        ) = 0;
};

} // end of namespace vlv

