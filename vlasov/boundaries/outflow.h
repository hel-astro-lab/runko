#pragma once

#include "../../corgi/corgi.h"
#include "../../em-fields/damping_tile.h"
#include "../../vlasov/tile.h"


namespace vlv {

namespace outflow {

template<size_t D, int S>
class Tile :
  virtual public fields::damping::Tile<D, S>,
  virtual public vlv::Tile<D>,
  virtual public fields::Tile<D>,
  virtual public corgi::Tile<D>
{

  public:

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
    corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz),
    fields::damping::Tile<D,S>(nx,ny,nz),
    vlv::Tile<D>(nx,ny,nz)
  { }

  /// destructor
  ~Tile() override = default;


  void stepLocation(corgi::Node<D>& /*grid*/) override {
      //std::cout<<"BC spatial step\n";

      //vlasov::AmrSpatialLagrangianSolver<Realf> ssol;
      //ssol.solve(*this, grid);
    }
    
  /// wall location
  using fields::damping::Tile<D, S>::fld1;
  using fields::damping::Tile<D, S>::fld2;

};

} // end of namespace outflow
} // end of namespace vlasov
