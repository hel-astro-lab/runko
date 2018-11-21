#pragma once

#include <array>
#include <vector>

#include "../definitions.h"


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

    void load_zero_block(size_t q);

    void load_block(size_t q, vblock_t block);

    void load_grid(size_t q, Realf val);

    std::vector<Realf> get_grid();

    std::vector<Realf> get_pencil();

    bool is_non_zero(size_t q);

    vblock_t get_slice(size_t q);

    Realf get_dx(size_t q);

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

      void set_bundle(Bundle _bundle);

      Bundle get_bundle( );

      void set_delta( Bundle _delta );

      vblock_t get_delta_slice(size_t i);

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
