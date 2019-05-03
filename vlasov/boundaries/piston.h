#pragma once

#include "../../corgi/corgi.h"
#include "../../vlasov/tile.h"


namespace vlv {

namespace piston {

template<size_t D>
class Tile :
  virtual public vlv::Tile<D>
{

  public:

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
    fields::Tile<D>(nx,ny,nz),
    vlv::Tile<D>(nx,ny,nz)
  { }


  void step_location(corgi::Node<D>& /*grid*/) override {
      //std::cout<<"BC spatial step\n";

      //vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
      //ssol.solve(*this, grid);
    }
    
  void reflect(corgi::Node<D>& grid);

};

} // end of namespace piston
} // end of namespace vlv
