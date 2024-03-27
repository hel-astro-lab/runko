#pragma once

#include "external/corgi/corgi.h"
#include "emf/boundaries/damping_tile.h"
#include "vlv/tile.h"
#include "vlv/spatial-solvers/amr_spatial_solver.h"

namespace vlv {

namespace outflow {

template<size_t D, int S>
class Tile :
  virtual public emf::damping::Tile<D, S>,
  virtual public vlv::Tile<D>
{

  public:

  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
    emf::Tile<D>(nx,ny,nz),
    emf::damping::Tile<D,S>(nx,ny,nz),
    vlv::Tile<D>(nx,ny,nz)
  { }

  /// switch flag to determine if we move plasma or not
  bool advance = false;

  void step_location(corgi::Grid<D>& grid) override {
      //std::cout<<"BC spatial step\n";

    if(advance) {
      vlv::AmrSpatialLagrangianSolver<float> ssol;
      ssol.solve(*this, grid);
    }
  }
    
  /// wall location
  using emf::damping::Tile<D, S>::fld1;
  using emf::damping::Tile<D, S>::fld2;

};

} // end of namespace outflow
} // end of namespace vlv
