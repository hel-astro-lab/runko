#include "bwd_lagrangian.h"


#include <cmath> 
#include <Eigen/Dense>

#include "../amr/mesh.h"
#include "../amr/numerics.h"
#include "../amr/refiner.h"



using namespace Eigen;


namespace vlv {
  namespace tools {
template<typename T>
T modf(T arg, T* iptr) = delete;

template<>
float modf(float arg, float* iptr) { return modff(arg, iptr); }

template<>
double modf(double arg, double* iptr) { return modf(arg, iptr); }

}
}


template<typename T, int D, int V>
void vlv::AmrMomentumLagrangianSolver<T,D,V>::solveMesh( 
        toolbox::AdaptiveMesh<T, 3>& mesh0,
        toolbox::AdaptiveMesh<T, 3>& mesh1,
        vec& Einc,
        vec& Binc,
        tools::Params<T>& params)
{

  toolbox::Adapter<T,3> adapter;
  adapter.cells_to_refine.clear();
  adapter.cells_to_unrefine.clear();


  // empty the target mesh
  // TODO: is this efficient; should we recycle instead?
  mesh1.clear();


  std::array<uint64_t, 3> index;
  // std::array<T, 3> grad;

  Vector3f B(Binc.data());  
  Vector3f E(Einc.data());  


  // level zero fill
  auto val     = T(0);
  T max_val = mesh0.max_value();

  T refine_indicator, unrefine_indicator;
  auto len = mesh1.get_size(0);

  // XXX remove higher-dimension updates (compresses inner loops away)
  // TODO does not work if the fill is not correctly done
  for(int i = 2; i>V-1; i--) len[i] = 1;

  // normalize
  //T norm0 = integrate_moment( mesh0,
  //          [](std::array<T,3>& uvel) -> T { return T(1);}
  //          );


  // cfl bound guards
  auto cids = mesh0.get_cells(true);
    
  // if mesh is empty, bail out
  if (cids.empty()) return;
  auto min_ind = mesh0.get_indices( cids.front() );
  auto max_ind = mesh0.get_indices( cids.back()  );

  int cfl_halo = 10; // how many CFL steps are allowed in backward substitution

  for(uint64_t t=0; t<len[2]; t++) {
    index[2] = t;
    for(uint64_t s=0; s<len[1]; s++) {
      index[1] = s;
      for(uint64_t r=0; r<len[0]; r++) {
        index[0] = r;

        // cfl bound guards
        if( r < (min_ind[0] - cfl_halo) ) continue;
        if( r > (max_ind[0] + cfl_halo) ) continue;

        uint64_t cid = mesh1.get_cell_from_indices(index, 0);
        val = backward_advect(index, 0, mesh0, E, B, params);

        // refinement
        // TODO specialize to value & gradient instead of just value
        refine_indicator   = val/max_val;
        unrefine_indicator = val/max_val;
        adapter.checkCell(mesh1, cid, refine_indicator, unrefine_indicator);

        mesh1.set(cid, val);
      }
    }
  }



  // create new leafs
  // TODO fixme to use max ref. lvl of mesh
  // for(size_t sweep=1; sweep<=mesh1.maximum_refinement_level; sweep++){
  for(int sweep=1; sweep<=mesh1.top_refinement_level; sweep++){

    adapter.refine(mesh1);

    // fmt::print(" sol: cells created {}\n", adapter.cells_to_refine.size());
    // fmt::print(" sol: cells removed {}\n", adapter.cells_removed.size());

    adapter.cells_to_refine.clear();
    adapter.cells_to_unrefine.clear();

    // next we refine
    for(auto&& cid : adapter.cells_created) {
      int rfl = mesh1.get_refinement_level(cid);
      auto index2 = mesh1.get_indices(cid);

      // fmt::print("creating {} at {}\n", cid, rfl);

      val = backward_advect(index2, rfl, mesh0, E, B, params);
      mesh1.set(cid, val);

      refine_indicator   = val/max_val;
      unrefine_indicator = val/max_val;
      adapter.checkCell(mesh1, cid, refine_indicator, unrefine_indicator);
    }

    // now unrefine 
    adapter.unrefine(mesh1);

  }


  // normalize back to original weight
  //T norm1 = integrate_moment( mesh1,
  //          [](std::array<T,3>& uvel) -> T { return T(1);}
  //          );
  //mesh1 *= norm0/norm1;


  return;
}

// backward advected Lorentz force
template<typename T, int D, int V>
T vlv::AmrMomentumLagrangianSolver<T,D,V>::backward_advect(
    std::array<uint64_t, 3>& index,
    int rfl,
    const toolbox::AdaptiveMesh<T, 3>& mesh0,
    //toolbox::AdaptiveMesh<T, 3>& mesh1,
    Vector3f& E,
    Vector3f& B,
    tools::Params<T>& params)
{
  T val; // return value
  vec u    = mesh0.get_center(index, rfl);  // velocity
  vec du   = mesh0.get_length(rfl);         // box size, i.e., \Delta u
  auto len = mesh0.get_size(rfl);


  // get shift of the characteristic solution from Lorentz force
  Vector3f uvel( u.data() );
  Vector3f F = lorentz_force(uvel, E, B, params.qm, params.cfl);

  // add other forces; default to zero 
  Vector3f Fi = other_forces(uvel, params);
  F += Fi;


  // advection in units of cells 
  // NOTE: We normalize with CFL because velocities are in units 
  // of c and we need them in units of grid speed.
  // XXX: CHECK
  std::array<T,3> shift, cell_shift;
  for(int i=0; i<V; i++) shift[i] = F(i) / du[i];

  // advected tiles
  std::array<int,3> index_shift;
  for(int i=0; i<V; i++) index_shift[i] = static_cast<int>( trunc(shift[i]) );

  // advection inside cell
  T tmp;
  for(int i=0; i<V; i++) cell_shift[i] = tools::modf<T>(shift[i], &tmp);

  // new grid indices
  std::array<uint64_t, 3> index_new;
  for(int i=0;   i<V;   i++) index_new[i] = index[i] + index_shift[i];
  for(int i = 2; i>V-1; i--) index_new[i] = index[i];


  // set boundary conditions (zero outside the grid limits)
  for(int i=0; i<V; i++){
    if( (index_new[i] <2) || (index_new[i] >= len[i]-2) ) return T(0);
  }

  // interpolation branch
  val = toolbox::interp_cubic<T,V>(mesh0, index_new, cell_shift, rfl);
  //val = toolbox::interp_linear<T,V>(mesh0, index_new, cell_shift, rfl);

  return val;
}

/// Relativistic Lorentz force
template<typename T, int D, int V>
inline Vector3f vlv::AmrMomentumLagrangianSolver<T,D,V>::lorentz_force(
  Vector3f& /*uvel*/,
  Vector3f& E,
  Vector3f& /*B*/,
  T qm,
  T cfl)
{
  // electromagnetic combined push
  /*
  Vector3f Ehalf = dt*E/2; // half step in E
  Vector3f us = uvel - Ehalf; // P^*, i.e., velocity in the middle of the step

  // B-field rotation matrix
  Vector3f b = B.normalized();
  T gamma = sqrt(1.0 + us.transpose()*us );
  T wt = dt*(qm*B.norm() / gamma); // relativistic cyclotron frequency

  Matrix3f Rot;
  Rot <<
    b(0)*b(0)+(1-b(0)*b(0))*cos(wt),    b(0)*b(1)*(1-cos(wt))-b(2)*sin(wt), b(0)*b(2)*(1-cos(wt))+b(1)*sin(wt),
    b(0)*b(1)*(1-cos(wt))+b(2)*sin(wt), b(1)*b(1)+(1-b(1)*b(1))*cos(wt),    b(1)*b(2)*(1-cos(wt))-b(0)*sin(wt),
    b(0)*b(2)*(1-cos(wt))-b(1)*sin(wt), b(1)*b(2)*(1-cos(wt))+b(0)*sin(wt), b(2)*b(2)+(1-b(2)*b(2))*cos(wt);

  return qm*( -Ehalf + Rot*us - uvel )*dt;
  */

  // electrostatic push
  //
  // Boris scheme for b=0 translates to
  // u = (cfl*u_0 + e + e)/cfl = u_0 + E/cfl
  //
  // with halving taken into account in definition of Ex
  return -qm*E/cfl;
}


/// default zero force for to be overloaded by more complicated solvers
template<typename T, int D, int V>
inline Vector3f vlv::AmrMomentumLagrangianSolver<T,D,V>::other_forces(
    Vector3f& /*uvel*/,
    tools::Params<T>& /*params*/
    )
{
  Vector3f ret = Vector3f::Zero();
  return ret;
}



//--------------------------------------------------
// explicit template instantiation
template class vlv::AmrMomentumLagrangianSolver<Realf, 1, 1>; //1D1V

