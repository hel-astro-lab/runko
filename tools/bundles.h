#pragma once

#include <array>
#include <vector>

#include "definitions.h"


namespace toolbox {
  /* \brief Bundle of pencils
   *
   *
   */
  class Bundle {

    /// values along the pencil
    std::vector<Realf> pencil;

    /// guiding grid for the pencil
    std::vector<Realf> grid;

    public:

    void resize( size_t N);

    size_t size();

    void loadZeroBlock(size_t q);

    void loadBlock(size_t q, vblock_t block);

    void loadGrid(size_t q, Realf val);

    std::vector<Realf> getGrid();

    std::vector<Realf> getPencil();

    bool isNonZero(size_t q);

    vblock_t getSlice(size_t q);

    Realf getDx(size_t q);

  }; // end of Bundle class



  /// Abstract base class for bundle interpolator
  class BundleInterpolator {
    public:
      /// internal bundle that we interpolate
      Bundle bundle;

      /// force acting on the fluid
      Bundle delta;

      /// time step
      Realf dt = 0.0;

      /// numerical zero
      Realf nzero = 1.0e-4;


      virtual ~BundleInterpolator() { };

      void setBundle(Bundle _bundle);

      Bundle getBundle( );

      void setDelta( Bundle _delta );

      vblock_t getDeltaSlice(size_t i);

      virtual Bundle interpolate() = 0;
  };



  /// Second order Lagrangian interpolator
  class BundleInterpolator2nd : public BundleInterpolator {
    public:
      Bundle interpolate( );
  };


  /// Fourth order Lagrangian interpolator
  class BundleInterpolator4th : public BundleInterpolator {
    public:
      Bundle interpolate( );
  };

  /// Fourth order Lagrangian positive non-oscillatory scheme
  class BundleInterpolator4PIC : public BundleInterpolator {
    public:
      Bundle interpolate( );
  };


} // end of toolbox namespace
