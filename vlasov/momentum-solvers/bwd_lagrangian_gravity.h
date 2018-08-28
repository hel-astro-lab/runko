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

    GravityAmrMomentumLagrangianSolver() {};

    virtual ~GravityAmrMomentumLagrangianSolver() = default;

    typedef std::array<T, 3> vec;

    /// Relativistic Lorentz force
    inline Vector3f lorentz_force(
        Vector3f& /*uvel*/,
        Vector3f& E,
        Vector3f& /*B*/,
        T qm,
        T cfl)
    override {

      std::cout << "using special Lorentz force, otherwise same\n";
        
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
