#pragma once

#include "external/corgi/corgi.h"
#include "core/vlv/tile.h"


namespace vlv {

namespace piston {

template<size_t D>
class Tile :
  virtual public vlv::Tile<D>
{

  public:

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
    emf::Tile<D>(nx,ny,nz),
    vlv::Tile<D>(nx,ny,nz)
  { }


  void step_location(corgi::Grid<D>& /*grid*/) override {
      //std::cout<<"BC spatial step\n";

      //vlv::AmrSpatialLagrangianSolver<float> ssol;
      //ssol.solve(*this, grid);
    }
    
  void reflect(corgi::Grid<D>& grid);

};

} // end of namespace piston
} // end of namespace vlv
