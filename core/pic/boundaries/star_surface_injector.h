#pragma once

#include <random>

#include "core/emf/boundaries/conductor.h"
#include "core/pic/tile.h"

namespace pic {

/// spherical stellar surface that injects particles close to the surface
template<size_t D>
class Star :
  public emf::Conductor<D>
{

private:

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float> uni_dis;

  // using raw pointer instead of smart ptrs; it does not take ownership of the object
  // so the container is not deleted when the temporary storage goes out of scope.
  using ConPtr = pic::ParticleContainer<D>* ;

public:

  // members
  using emf::Conductor<D>::radius;
  using emf::Conductor<D>::period;
  using emf::Conductor<D>::B0;
  using emf::Conductor<D>::chi_mu;
  using emf::Conductor<D>::chi_om;
  using emf::Conductor<D>::phase_mu;
  using emf::Conductor<D>::phase_om;
  using emf::Conductor<D>::cenx;
  using emf::Conductor<D>::ceny;
  using emf::Conductor<D>::cenz;
  using emf::Conductor<D>::delta;
  using emf::Conductor<D>::angular_velocity;
  using emf::Conductor<D>::radius_pc;
  using emf::Conductor<D>::delta_pc;
  using emf::Conductor<D>::Nx;
  using emf::Conductor<D>::Ny;
  using emf::Conductor<D>::Nz;

  // methods
  using emf::Conductor<D>::dipole;
  using emf::Conductor<D>::insert_em;
  using emf::Conductor<D>::update_b;
  using emf::Conductor<D>::update_e;

  double temp_pairs = 0.2; // pair injection temperature
  double temp_phots = 0.001; // photon injection temperature

  double ninj_pairs = 0.05; // pair injection rate (factor in front of E/q
  double ninj_phots = 0.0;  // photon injection rate per cell per step

  double ninj_min_pairs = 0.01; // minimum pairs per cell per step to inject
  double ninj_min_phots = 0.0;  // minimum photons per cell per step to inject

  int height_atms = 1; // height of the atmosphere in cells
  float wep = 1.0f;    // weight of added electrons and positrons
  float wph = 1.0f;    // weight of added photons

  Star() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)
  { }

  // random numbers between [0, 1[
  float rand() { return uni_dis(gen); };

  void solve(pic::Tile<D>&  tile);

};


} // end of ns pic
