#pragma once

#include "../../em-fields/boundaries/conductor.h"
#include "../tile.h"
#include <random>


namespace pic {

/// spherical stellar surface that injects particles close to the surface
template<size_t D>
class Star :
  public fields::Conductor<D>
{

private:

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float_p> uni_dis;

  // using raw pointer instead of smart ptrs; it does not take ownership of the object
  // so the container is not deleted when the temporary storage goes out of scope.
  using ConPtr = pic::ParticleContainer<D>* ;

public:

  // members
  using fields::Conductor<D>::radius;
  using fields::Conductor<D>::period;
  using fields::Conductor<D>::B0;
  using fields::Conductor<D>::chi_mu;
  using fields::Conductor<D>::chi_om;
  using fields::Conductor<D>::phase_mu;
  using fields::Conductor<D>::phase_om;
  using fields::Conductor<D>::cenx;
  using fields::Conductor<D>::ceny;
  using fields::Conductor<D>::cenz;
  using fields::Conductor<D>::delta;
  using fields::Conductor<D>::angular_velocity;
  using fields::Conductor<D>::radius_pc;
  using fields::Conductor<D>::delta_pc;
  using fields::Conductor<D>::Nx;
  using fields::Conductor<D>::Ny;
  using fields::Conductor<D>::Nz;

  // methods
  using fields::Conductor<D>::dipole;
  using fields::Conductor<D>::insert_em;
  using fields::Conductor<D>::update_b;
  using fields::Conductor<D>::update_e;

  double temp_pairs = 0.2; // pair injection temperature
  double temp_phots = 0.001; // photon injection temperature

  double ninj_pairs = 0.05; // pair injection rate (factor in front of E/q
  double ninj_phots = 0.0;  // photon injection rate per cell per step

  double ninj_min_pairs = 0.01; // minimum pairs per cell per step to inject
  double ninj_min_phots = 0.0;  // minimum photons per cell per step to inject

  Star() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)
  { }

  // random numbers between [0, 1[
  float_p rand() { return uni_dis(gen); };

  void solve(pic::Tile<D>&  tile);

};


} // end of ns pic
