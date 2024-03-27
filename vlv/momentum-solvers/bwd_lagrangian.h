#pragma once

#include "vlv/momentum-solvers/amr_momentum_solver.h"
//#include <Eigen/Dense>

namespace vlv {

//using namespace Eigen;

/// \brief back-substituting semi-Lagrangian adaptive advection solver
template<typename T, int D, int V>
class AmrMomentumLagrangianSolver : 
  virtual public MomentumSolver<T, D, V> 
{

  public:
    using typename MomentumSolver<T,D,V>::vec;
    using typename MomentumSolver<T,D,V>::Vec3E;


    AmrMomentumLagrangianSolver() = default;

    ~AmrMomentumLagrangianSolver() override = default;

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
        Vec3E& E,
        Vec3E& B,
        tools::Params<T>& params
        );

    virtual Vec3E lorentz_force(
        Vec3E& /*uvel*/,
        Vec3E& E,
        Vec3E& /*B*/,
        T qm,
        T cfl);

    virtual Vec3E other_forces(
        Vec3E& uvel,
        tools::Params<T>& params
        );
};


} // end of namespace vlv
