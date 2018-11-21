#pragma once

#include "amr_momentum_solver.h"
#include <Eigen/Dense>

namespace vlv {

using namespace Eigen;

/// \brief back-substituting semi-Lagrangian adaptive advection solver
template<typename T, int D, int V>
class AmrMomentumLagrangianSolver : 
  virtual public MomentumSolver<T, D, V> 
{

  public:

    typedef std::array<T, 3> vec;

    AmrMomentumLagrangianSolver() {};

    virtual ~AmrMomentumLagrangianSolver() = default;

    void solve_mesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& Einc,
        vec& Binc,
        tools::Params<T>& params
        ) override;

    T backward_advect(
        std::array<uint64_t, 3>& index,
        int rfl,
        const toolbox::AdaptiveMesh<T, 3>& mesh0,
        Vector3f& E,
        Vector3f& B,
        tools::Params<T>& params
        );

    virtual Vector3f lorentz_force(
        Vector3f& /*uvel*/,
        Vector3f& E,
        Vector3f& /*B*/,
        T qm,
        T cfl);

    virtual Vector3f other_forces(
        Vector3f& uvel,
        tools::Params<T>& params
        );
};


} // end of namespace vlv
