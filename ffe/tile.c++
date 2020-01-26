#include <iostream>

#include "tile.h"


namespace ffe {
  using namespace mpi4cpp;


/// Get current integration time snapshot of Yee lattice
// FIXME: does not need to be overloaded
template<std::size_t D>
fields::YeeLattice& Tile<D>::get_yee(int /*i*/) 
{
  return this->yee.at(0);
}


template<std::size_t D>
ffe::SkinnyYeeLattice& Tile<D>::get_step(int n) 
{
  if     (n == 0) return this->step0; 
  else if(n == 1) return this->step1; 
  else if(n == 2) return this->step2; 
  else if(n == 3) return this->step3; 

  // silently return step0 else; 
  // in reality bindings captures this and throw
  return this->step0;
}


//--------------------------------------------------
// explicit template instantiation
template class Tile<2>;
//template class Tile<3>;

} // end of ns ffe
