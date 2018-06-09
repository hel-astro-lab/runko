#pragma once

#include <iostream>
#include <vector>
#include <tuple>

#include "cell.h"
#include "communicate.h"

#include "fftw3.h"


namespace pic {

//! Spatial current filter using fftw
//
// Because of the fftw3 usage this class is internally heavily relying on C code
// instead of C++.
class Filter {

  /// internal indexing, 
  // NOTE: this is different from the rest of the code but done like this
  // in order to appear consistent with the fttw examples.
  inline int index(int i, int j) {
    return i*width + j;
  }


  /// neighbor-padded arrays for convolution
  fftw_complex *jx, *jy, *jz;

  /// actual (zero-padded) convolution kernel
  fftw_complex *kernel;

  /// fftw3 transform plans
  fftw_plan p_kernel, p_forw, p_back;



  public:

  /// image width
  int width;

  /// image height
  int height;

  // image depth
  int depth=0;

  /// Build filter assuming information from 3x3x1 tiles
  Filter(int NxMesh, int NyMesh) :
      width( 3*NxMesh),
      height(3*NyMesh),
      depth( 1) 
  {

    // allocate
    jx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * height * width);
    jy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * height * width);
    jz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * height * width);

    kernel = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * height * width);


    // in-place complex to complex transform
    // TODO: change to real2complex and complex2real plans
    // TODO: use fftw guru interface to apply same plan for all arrays (in same direction)
    p_kernel = fftw_plan_dft_2d(height, width, kernel, kernel, FFTW_FORWARD,  FFTW_MEASURE);

    p_forw   = fftw_plan_dft_2d(height, width, jx, jx, FFTW_FORWARD,  FFTW_MEASURE);
    p_back   = fftw_plan_dft_2d(height, width, jx, jx, FFTW_BACKWARD, FFTW_MEASURE);

  }


  /// explicit destructor for fftw arrays and plans
  virtual ~Filter() {
    fftw_free(jx);
    fftw_free(jy);
    fftw_free(jz);
    fftw_free(kernel);

    fftw_destroy_plan(p_forw);
    fftw_destroy_plan(p_back);
  }


  /// zero-padded centralized index
  //
  //NOTE:
  // The FFT is relative to the middle and in its scale there are negative points.
  // In the memory the points are 0...n-1, but the FFT treats them as 
  // -ceil(n/2)...floor(n/2), 
  // where 0 is -ceil(n/2) and n-1 is floor(n/2)
  //
  inline virtual std::tuple<int,int> zero_padded_index(int i, int j) 
  {
    int iff = i - ceil(height/2);
    int jff = j - ceil(width /2);
    
    return std::make_tuple(iff, jff);
  }

  /// initialize kernel given the rank
  virtual void init_kernel() 
  {

    double val;
    auto zcenter = std::make_tuple(0,0); // zero-padded center
    for(int i  = 0 ; i < height ; ++i) {
      for(int j = 0 ; j < width ; ++j) {
        auto zindx = zero_padded_index(i,j);
        val = zindx == zcenter ? 1.0 : 0.0;
        
        kernel[ index(i,j) ][0] = val; // real part
        kernel[ index(i,j) ][1] = 0.0; // complex part
      }
    }
  }

  // normalize fft transformation
  void normalize() 
  {
    for(int i  = 0 ; i < height ; ++i)  
      for(int j = 0 ; j < width ; ++j)  
        jx[ index(i,j) ][0] /= width*height;

  }


  /// transformation
  virtual void fft_kernel()         { fftw_execute(p_kernel); }

  virtual void fft_image_forward()  { fftw_execute(p_forw);   }

  virtual void fft_image_backward() { 
    fftw_execute(p_back);
    normalize();
  }


  /// Copy currents from neighbors into fftw array
  void get_padded_current(
      pic::PicCell& cell, 
      corgi::Node& grid
      )
  {


    int k = 0;
    for (int i=-1; i<=1; i++)
    for (int j=-1; j<=1; j++) {
    //for (int k=-1; k<=1; k++) { // TODO: hack to get 2d tiles working
      std::cout << "from: (" << i << "," << j << "," << k << ")" << '\n';
      
      if( (i == 0) & (j == 0) ){
        pic::ParticleBlock& neigh = cell.container; // local cell
      } else {
        pic::ParticleBlock& neigh  = get_external_data(i, j, cell, grid); //external neighbor
      }


    }
  
  }


  // --------------------------------------------------
  // auxiliary/utility functions for debugging

  void set_kernel(std::vector<double>& in)
  {
    if(in.size() != (size_t)height*width) std::cout << "error in size!\n";

    int q = 0;
    for(int i = 0 ; i < height ; ++i)  
    for(int j = 0 ; j < width ; ++j, q++)  
      kernel[ index(i,j) ][0] = in[q];
  }


  void set_image(std::vector<double>& in)
  {
    if(in.size() != (size_t)height*width) std::cout << "error in size!\n";

    int q = 0;
    for(int i = 0 ; i < height ; ++i)  
    for(int j = 0 ; j < width ; ++j, q++)  
      jx[ index(i,j) ][0] = in[q];
  }


  std::vector<double> get_kernel()
  {
    std::vector<double> ret;

    for(int i = 0 ; i < height ; ++i)  
    for(int j = 0 ; j < width ; ++j)  
      ret.push_back( kernel[ index(i,j) ][0] );

    return ret;
  }


  std::vector<double> get_image()
  {
    std::vector<double> ret;

    for(int i = 0 ; i < height ; ++i)  
    for(int j = 0 ; j < width ; ++j)  
      ret.push_back( jx[ index(i,j) ][0] );

    return ret;
  }




};

} // end of namespace pic
