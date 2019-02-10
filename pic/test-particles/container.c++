#include "container.h"
#include <mpi.h>
#include <cmath>




pic::TestParticleContainer::TestParticleContainer() :
    ParticleContainer() 
{
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);

  // NOTE: we assume this in id key generation
  assert(mpi_world_size < 1e6);

  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

}

double pic::TestParticleContainer::keygen() 
{

  // base
  double unique_key = static_cast<double>(key);
  key++;

  // rank separator
  unique_key += static_cast<double>(rank)/1.0e6;

  return unique_key;
}

void pic::TestParticleContainer::add_test_particle (
    std::vector<double> prtcl_loc,
    std::vector<double> prtcl_vel)
{
  // get unique running key
  double unique_key = keygen();

  // add normally
  ParticleContainer::add_particle(prtcl_loc, prtcl_vel, unique_key);
}













