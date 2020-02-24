#pragma once

#include "bwd_lagrangian.h"

namespace vlv {


using namespace Eigen;


/// \brief back-substituting semi-Lagrangian adaptive advection solver with gravity
template<typename T, int D, int V>
class GravityAmrMomentumLagrangianSolver : 
  virtual public AmrMomentumLagrangianSolver<T,D,V> 
{

  public:
    using typename MomentumSolver<T,D,V>::vec;
    using typename MomentumSolver<T,D,V>::Vec3E;

    T g0; // strength of gravity
    T Lx; // box size

    GravityAmrMomentumLagrangianSolver(
        T g0,
        T Lx
        ) :
      g0(g0), Lx(Lx)
    {};

    ~GravityAmrMomentumLagrangianSolver() override = default;

    /// Gravity
    inline Vec3E other_forces(
        Vec3E& uvel,
        vlv::tools::Params<T>& params) 
    override
    {
      T gam  = std::sqrt(1.0 + uvel.dot(uvel));
      //std::cout << "using special other force at xloc=" << params.xloc << " gam=" << gam;
      //std::cout << " with g0=" << g0 << " and Lx=" << Lx << "\n";

      //Vec3E ret( -g0*gam*(2.0*params.xloc/Lx - 1.0), 0, 0); // flux tube
      Vec3E ret( -g0*gam*(params.xloc/Lx)/params.cfl, 0, 0); // atmosphere
      return ret;
    }

    /// Relativistic Lorentz force
    inline Vec3E lorentz_force(
        Vec3E& /*uvel*/,
        Vec3E& E,
        Vec3E& /*B*/,
        T qm,
        T cfl)
    override {
      //std::cout << "using special Lorentz force, otherwise same\n";
        
      // electrostatic push
      //
      // Boris scheme for b=0 translates to
      // u = (cfl*u_0 + e + e)/cfl = u_0 + E/cfl
      //
      // with halving taken into account in definition of Ex
      return -qm*E/cfl;
    }

};


} // end of namespace vlv
