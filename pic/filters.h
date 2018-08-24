#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>

#include "tile.h"
#include "communicate.h"

#include "../units.h"
#include "fftw3.h"



namespace pic {


inline fields::YeeLattice& get_neighbor_yee(
    int i, int j,
    pic::Tile<2>& tile, 
    corgi::Node<2>& node)
{
  auto cneigh = std::dynamic_pointer_cast<fields::Tile<2>>(
        node.getTilePtr( tile.neighs(i, j) ));
  return cneigh->getYee();
}



using units::pi;

//! Spatial current filter using fftw
//
// Because of the fftw3 usage this class is internally heavily relying on C code
// instead of C++.
class Filter {

  private:

  /// image size
  int Nx, Ny, Nz;

    
  /// circular/wrap indexing
  inline int circular(int x, int M)
  {
    if(x<0)    return x+M;
    if(x >= M) return x-M;
    return x;
  }

  /// internal circular indexing, 
  // NOTE: this is different from the rest of the code but done like this
  // in order to appear consistent with the fttw examples.
  inline int index(int i, int j) {
    i = circular(i, Nx);
    j = circular(j, Ny);

    assert(i >= 0 && i < Nx);
    assert(j >= 0 && j < Ny);

    return i + j*Nx;
  }


  /// neighbor-padded arrays for convolution
  fftwf_complex *jx, *jy, *jz;

  /// actual (zero-padded) convolution kernel
  fftwf_complex *kernel;

  /// fftw3 transform plans
  fftwf_plan p_kernel, 
            p_forw_jx, p_forw_jy, p_forw_jz, 
            p_back_jx, p_back_jy, p_back_jz;

  public:

  /// Build filter assuming information from 3x3x1 tiles
  Filter(int NxMesh, int NyMesh) :
      Nx( 3*NxMesh ),
      Ny( 3*NyMesh ),
      Nz( 1) 
  {

    // allocate
    //jx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny*Nz);
    //jy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny*Nz);
    //jz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny*Nz);
    jx = (fftwf_complex*) fftwf_alloc_complex(Nx*Ny*Nz);
    jy = (fftwf_complex*) fftwf_alloc_complex(Nx*Ny*Nz);
    jz = (fftwf_complex*) fftwf_alloc_complex(Nx*Ny*Nz);

    for(int j = 0 ; j < Ny ; ++j) {
      for(int i = 0 ; i < Nx ; ++i) {
        jx[ index(i,j) ][0] = 0.0;
        jx[ index(i,j) ][1] = 0.0;

        jy[ index(i,j) ][0] = 0.0;
        jy[ index(i,j) ][1] = 0.0;

        jz[ index(i,j) ][0] = 0.0;
        jz[ index(i,j) ][1] = 0.0;
      }
    }


    //kernel = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny*Nz);
    kernel = (fftwf_complex*) fftwf_alloc_complex(Nx*Ny*Nz);


    // in-place complex to complex transform
    // TODO: change to real2complex and complex2real plans
    // TODO: use fftw guru interface to apply same plan for all arrays (in same direction)
    p_kernel    = fftwf_plan_dft_2d( Nx,Ny, kernel, kernel, FFTW_FORWARD,  FFTW_MEASURE);

    p_forw_jx   = fftwf_plan_dft_2d( Nx,Ny, jx, jx, FFTW_FORWARD,  FFTW_MEASURE );
    p_forw_jy   = fftwf_plan_dft_2d( Nx,Ny, jy, jy, FFTW_FORWARD,  FFTW_MEASURE );
    p_forw_jz   = fftwf_plan_dft_2d( Nx,Ny, jz, jz, FFTW_FORWARD,  FFTW_MEASURE );

    p_back_jx   = fftwf_plan_dft_2d( Nx,Ny, jx, jx, FFTW_BACKWARD, FFTW_MEASURE );
    p_back_jy   = fftwf_plan_dft_2d( Nx,Ny, jy, jy, FFTW_BACKWARD, FFTW_MEASURE );
    p_back_jz   = fftwf_plan_dft_2d( Nx,Ny, jz, jz, FFTW_BACKWARD, FFTW_MEASURE );


    /*
    fftw_print_plan(p_forw_jx);
    fftw_print_plan(p_forw_jy);
    fftw_print_plan(p_forw_jy);

    fftw_print_plan(p_back_jx);
    fftw_print_plan(p_back_jy);
    fftw_print_plan(p_back_jz);
    */
  }


  /// explicit destructor for fftw arrays and plans
  virtual ~Filter() {
    fftwf_free(jx);
    fftwf_free(jy);
    fftwf_free(jz);
    fftwf_free(kernel);

    fftwf_destroy_plan(p_forw_jx);
    fftwf_destroy_plan(p_forw_jy);
    fftwf_destroy_plan(p_forw_jz);

    fftwf_destroy_plan(p_back_jx);
    fftwf_destroy_plan(p_back_jy);
    fftwf_destroy_plan(p_back_jz);
  }


  /// zero-padded centralized index
  //
  //NOTE:
  // The FFT is relative to the middle and in its scale there are negative points.
  // In the memory the points are 0...n-1, but the FFT treats them as 
  // -ceil(n/2)...floor(n/2), 
  // where 0 is -ceil(n/2) and n-1 is floor(n/2)
  //
  //NOTE: this is not how fftw3 does its indexing anymore!
  //
  //inline virtual std::tuple<int,int> zero_padded_index(int i, int j) 
  //{
  //  int iff = ceil(height/2);
  //  int jff = ceil(width /2);
  //  
  //  return std::make_tuple(iff, jff);
  //}

  /// zero-wrapped index convention used by fftw3
  inline virtual std::tuple<int,int> zero_wrapped_index(int i, int j) 
  {
    int iff = i >= 0 ? i : Nx + i;
    int jff = j >= 0 ? j : Ny + j;

    return std::make_tuple(iff, jff);
  }


  /// initialize kernel 
  //
  // In practice this is just a complex way to set kernel[0,0] = 1,
  // but we do it like this for clarity as the syntax is used later on.
  virtual void init_kernel() 
  {
    Realf val;

    // kernel size (even number because of initialization)
    //int knx = height/3;
    //int kny = width /3;

    // halo region size
    int h1 = Nx/2; // division is floor(x/y) automatically
    int w1 = Ny/2;
    //int h2 = (height + 2 - 1)/2; // ceil(x/y)
    //int w2 = (width  + 2 - 1)/2;

    int wi,wj;

    for(int j=-w1; j<=w1 ; ++j) {
      for(int i =-h1; i<=h1 ; ++i) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);
          
        val = 0.0;
        if ((i ==  0) && (j ==  0)) val = 1.0;

        //if ((i ==  1) && (j ==  0)) val = 1.0;
        //if ((i ==  0) && (j ==  1)) val = 1.0;
        //if ((i ==  1) && (j ==  1)) val = 1.0;
        //if ((i == -1) && (j ==  0)) val = 1.0;
        //if ((i ==  0) && (j == -1)) val = 1.0;
        //if ((i == -1) && (j == -1)) val = 1.0;
        //if ((i == 1 ) && (j == -1)) val = 1.0;
        //if ((i ==-1 ) && (j ==  1)) val = 1.0;
        
        kernel[ index(wi,wj) ][0] = val; // real part
        kernel[ index(wi,wj) ][1] = 0.0; // complex part
      }
    }
  }


  /// Apply 3-point digital filter directly
  virtual void direct_convolve_3point()
  {

    // 3-point digital filter
    std::vector<Realf> coeffs = {{ 1., 2., 1.,
                                    2., 4., 2.,
                                    1., 2., 1. }};
    // normalize
    Realf norm = 0.0;
    for(auto& coeff : coeffs) norm += coeff;
    for(auto& coeff : coeffs) coeff /= norm;

    std::vector<Realf> image1, image2, image3;
    image1.resize(Nx*Ny*Nz); 
    image2.resize(Nx*Ny*Nz); 
    image3.resize(Nx*Ny*Nz); 

    for (int j=0; j < Ny;  ++j) {
      for (int i=0; i < Nx; ++i) {
        image1[ index(i,j) ] = jx[ index(i,j) ][0];
        image2[ index(i,j) ] = jy[ index(i,j) ][0];
        image3[ index(i,j) ] = jz[ index(i,j) ][0];
      }
    }

    direct_convolve(image1.data(), coeffs.data(), sqrt( coeffs.size() ) );
    direct_convolve(image2.data(), coeffs.data(), sqrt( coeffs.size() ) );
    direct_convolve(image3.data(), coeffs.data(), sqrt( coeffs.size() ) );

    // normalize and copy back
    norm = 1.0;
    for (int j=0; j < Ny;  ++j) {
      for (int i=0; i < Nx; ++i) {
        jx[ index(i,j) ][0] = image1[ index(i,j) ]/norm;
        jy[ index(i,j) ][0] = image2[ index(i,j) ]/norm;
        jz[ index(i,j) ][0] = image3[ index(i,j) ]/norm;
      }
    }


  }


  /// initialize 3-point digital filter into kernel array for FFT transformations
  virtual void init_3point_kernel(int times)
  {

    // kernel
    //int K = 3; // three point kernel

    // 3-point digital filter
    std::vector<Realf> coeffs = {{ 1., 2., 1.,
                                    2., 4., 2.,
                                    1., 2., 1. }};

    Realf norm = 0.0;
    for(auto& coeff : coeffs) norm += coeff;
    for(auto& coeff : coeffs) coeff /= norm;


    // create temporary Real number image array
    // NOTE: can not easily copy kernel into pure real part due to interleaved nature
    std::vector<Realf> image;
    image.resize(Nx*Ny*Nz); 

    for (int j=0; j < Ny;  ++j) 
    for (int i=0; i < Nx; ++i)
        image[ index(i,j) ] = kernel[ index(i,j) ][0];

    // perform convolution N times
    for(int N=0; N < times; N++) direct_convolve(image.data(), coeffs.data(), sqrt( coeffs.size() ) );
    
    // normalize
    // norm = 0.0;
    //for (int i=0; i < height; ++i)
    //  for (int j=0; j < width;  ++j) 
    //    norm += image[ index(i,j) ];
    norm = 1.0;

    // normalize and copy back
    for (int j=0; j < Ny;  ++j) 
      for (int i=0; i < Nx; ++i)
        kernel[ index(i,j) ][0] = image[ index(i,j) ]/norm;

  }


  /// Direct circular convolve two arrays
  // NOTE: assumes height x width for image size
  //
  // in, out are m x n images (integer data)
  // K is the kernel size (KxK) - currently needs to be an odd number, e.g. 3
  // coeffs[K][K] is a 2D array of integer coefficients
  // scale is a scaling factor to normalise the filter gain
  //
  void direct_convolve(
      Realf* image,
      Realf* kernel,
      int K)
  {

    // out array
    Realf data;
    std::vector<Realf> out;
    out.resize(Nx*Ny*Nz);

    //for (int j = K/2; j < width -K/2; ++j) // iterate through image
    for (int j=0; j < Ny; ++j) // iterate through circular image
    {
      //for (int i = K/2; i < height - K/2; ++i) // iterate through image
      for (int i=0; i < Nx; ++i) // iterate through circular image
      {
        Realf sum = 0.0; // sum will be the sum of input data * coeff terms

        // convolution of single point
        for (int jj = -K/2; jj <= K/2; ++jj)
        {
          for (int ii = -K/2; ii <= K/2; ++ii) // iterate over kernel
          {
            data = image[ index(i+ii, j+jj) ]; 
            Realf coeff = kernel[ (ii + K/2)*K + (jj + K/2) ];

            sum += data * coeff;
          }
        }
        out[ index(i,j) ] = sum; // sum of convolution products and store in output
      }
    } // end of conv

    // copy (real part) back
    for (int j=0; j < Ny;  ++j) 
      for (int i=0; i < Nx; ++i)
        image[ index(i,j) ] = out[ index(i,j) ];

  } 




  /// initialize Gaussian kernel given the sigma
  virtual void init_gaussian_kernel(Realf sigx, Realf sigy) 
  {
    // kernel size (even number because of initialization)
    //int knx = height/3;
    //int kny = width /3;

    // halo region size
    int h1 = Nx/2; // division is floor(x/y) automatically
    int w1 = Ny /2;
    int h2 = (Nx + 2 - 1)/2; // ceil(x/y)
    int w2 = (Ny + 2 - 1)/2;

    int wi,wj;
    Realf val;

    Realf sum = 0.0;
    for(int j=-w1; j<w2; ++j) {
      for(int i =-h1; i<h2; ++i) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);

        assert(wi >= 0 && wi < Nx);
        assert(wj >= 0 && wj < Ny);

        val = 1.0;
        val *= exp( -0.5*((Realf)(i*i))/sigx/sigx);
        val *= exp( -0.5*((Realf)(j*j))/sigy/sigy);

        kernel[ index(wi,wj) ][0] = val; // real part
        kernel[ index(wi,wj) ][1] = 0.0; // complex part
        sum += val;
      }
    }

    // normalize 
    //for(int i=0; i<height; ++i)
    //for(int j=0; j<width;  ++j)
    //  kernel[ index(i,j) ][0] /= sum;

  }


  /// initialize 2d box sinc filter (in frequency space)
  /*
  virtual void init_sinc_kernel(double X, double Y) 
  {
    double val, u,v;
    int h2 = floor(height/2);
    int w2 = floor(width /2);
    int wi, wj;

    double sum = 0.0;
    for(int i =-h2; i<=h2 ; ++i) {
      for(int j=-w2; j<=w2 ; ++j) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);

        u = (double)i;
        v = (double)j;

        val = 0.0;
        if ((abs(i) < height/3) || (abs(j) < width/3)){
          val = X*sin(pi*X*u)/(pi*X*u) * 
                Y*sin(pi*Y*v)/(pi*Y*v);
        }
        //val = (X*sin(pi*X*u)/(pi*X*u)) * (Y*sin(pi*Y*v)/(pi*Y*v));
        if ((i == 0) && (j == 0)) val = 4.0*pi*pi*X*Y; // center to avoid 0/0
        else if (i == 0) val = 2.0*pi*X*Y*sin(pi*Y*v)/(pi*Y*v);
        else if (j == 0) val = 2.0*pi*X*Y*sin(pi*X*u)/(pi*X*u);

        // windowing
        val *= 0.42 - 0.5*cos(2*pi*3.0/height) + 0.08*cos(4.0*pi*3.0/height);
        val *= 0.42 - 0.5*cos(2*pi*3.0/width ) + 0.08*cos(4.0*pi*3.0/width);

        kernel[ index(wi,wj) ][0] = val; // real part
        kernel[ index(wi,wj) ][1] = 0.0; // complex part
        sum += val;
      }
    }

    // normalize 
    for(int i =0; i<height ; ++i)
    for(int j=0; j<width ; ++j)
      kernel[ index(i,j) ][0] /= sum;

  }
  */


  /// Low-pass filter in frequency space
  // uses cutoff to describe how many array elements are filtered
  virtual void init_lowpass_fft_kernel(int cutoff) 
  {

    // kernel size (even number because of initialization)
    //int knx = height/3;
    //int kny = width /3;

    // halo region size
    int h1 = Nx/2; // division is floor(x/y) automatically
    int w1 = Ny/2;
    int h2 = (Nx + 2 - 1)/2; // ceil(x/y)
    int w2 = (Ny + 2 - 1)/2;

    int wi,wj;
    Realf val;

    Realf sum = 0.0;
    for(int i =-h1; i<h2; ++i) {
      for(int j=-w1; j<w2; ++j) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);

        assert(wi >= 0 && wi < Nx);
        assert(wj >= 0 && wj < Ny);

        val = 0.0;
        if (sqrt(i*i + j*j) < cutoff) val = 1.0; // circle (Bessel)
        //if ((sqrt(i*i) < cutoff) && (sqrt(j*j) < cutoff)) val = 1.0; // box (2D Sinc)

        kernel[ index(wi,wj) ][0] = val; // real part
        kernel[ index(wi,wj) ][1] = 0.0; // complex part
        sum += val;
      }
    }
  }


  // normalize fft transformation
  void normalize() 
  {
    for(int k = 0 ; k < Nz ; ++k) {
      for(int j = 0 ; j < Ny ; ++j) {
        for(int i  = 0 ; i < Nx ; ++i) {
          jx[ index(i,j) ][0] /= Nx*Ny*Nz;
          jy[ index(i,j) ][0] /= Nx*Ny*Nz;
          jz[ index(i,j) ][0] /= Nx*Ny*Nz;
        }
      }
    }
  }




  /// FFT kernel (once is enough)
  virtual void fft_kernel()         { fftwf_execute(p_kernel); }

  /// FFT image forward
  virtual void fft_image_forward()  
  { 
    fftwf_execute(p_forw_jx);   
    fftwf_execute(p_forw_jy);   
    fftwf_execute(p_forw_jz);   
  }

  /// FFT image backwards and normalize 
  virtual void fft_image_backward() { 
    fftwf_execute(p_back_jx);
    fftwf_execute(p_back_jy);
    fftwf_execute(p_back_jz);
    normalize();
  }


  /// Multiply kernel and image
  void apply_kernel()
  {
    Realf x1, y1, x2, y2, x3, y3;
    Realf u, v;
    for(int j = 0 ; j < Ny ; ++j) {
      for(int i  = 0 ; i < Nx ; ++i) {

        x1 = jx[ index(i,j) ][0];
        y1 = jx[ index(i,j) ][1];

        x2 = jy[ index(i,j) ][0];
        y2 = jy[ index(i,j) ][1];

        x3 = jz[ index(i,j) ][0];
        y3 = jz[ index(i,j) ][1];

        u = kernel[ index(i,j) ][0];
        v = kernel[ index(i,j) ][1];

        jx[ index(i,j) ][0] = x1*u - y1*v;
        jx[ index(i,j) ][1] = x1*v + y1*u;

        jy[ index(i,j) ][0] = x2*u - y2*v;
        jy[ index(i,j) ][1] = x2*v + y2*u;

        jz[ index(i,j) ][0] = x3*u - y3*v;
        jz[ index(i,j) ][1] = x3*v + y3*u;
      }
    }


  }


  /// Copy currents from neighbors into fftw array
  void get_padded_current(
      pic::Tile<2>& tile, 
      corgi::Node<2>& node
      )
  {

    //int NxMesh = tile.Nx;
    //int NyMesh = tile.Ny;
    //int NzMesh = tile.Nz;
    int indx;

    //int k = 0;
    for (int i=-1; i<=1; i++)
    for (int j=-1; j<=1; j++) {
    //for (int k=-1; k<=1; k++) { // TODO: hack to get 2d tiles working
      //std::cout << "from: (" << i << "," << j << "," << k << ")" << '\n';
      
      fields::YeeLattice& mesh = ((i==0)&&(j==0)) ? tile.getYee() : get_neighbor_yee(i,j,tile,node);

      /*
      std::cout << "--------------------------------------------------\n";
      std::cout << " Nx: " << mesh.Nx;
      std::cout << " Ny: " << mesh.Ny;
      std::cout << " Nz: " << mesh.Nz;
      std::cout << " i: " << i;
      std::cout << " j: " << j;
      std::cout << "\n";
      std::cout << "==============================\n";
      */

      int s = 0;
      for(int r=0; r<(int)mesh.Ny; r++) {
        for(int q=0; q<(int)mesh.Nx; q++) {
          //for(int s=0; r<Nz; s++) {
          
          indx = index((i+1)*mesh.Nx + q, (j+1)*mesh.Ny + r);

          /*
      std::cout << " q: " << q;
      std::cout << " r: " << r;
      std::cout << " s: " << s;
      std::cout << " indx: " << indx;
      std::cout << " ix: " << (i+1)*mesh.Nx + q;
      std::cout << " jy: " << (j+1)*mesh.Ny + r;
      */

          jx[ indx ][0] = mesh.jx(q,r,s);
          jy[ indx ][0] = mesh.jy(q,r,s);
          jz[ indx ][0] = mesh.jz(q,r,s);

          /*
      std::cout << " cur jx: " << jx[indx][0];
      std::cout << " cur jy: " << jy[indx][0];
      std::cout << " cur jz: " << jz[indx][0];
      std::cout << "\n";
      */

        }
      }


    }
  }

  /// set (filtered) current back to Yee lattice
  // NOTE: we export the current to jx1/jy1/jz1 instead of jx/jy/jz to 
  //       ensure thread safety
  void set_current( pic::Tile<2>& tile)
  {
    fields::YeeLattice& mesh = tile.getYee();

    int indx;

    // set only "central" tile values
    int i = 0;
    int j = 0;

    int halo = 3;
    int s = 0; // TODO: third index for 3D case

    //for(int s=0; r<Nz; s++) {
    for(int r=-halo; r<(int)mesh.Ny+halo; r++) {
      for(int q=-halo; q<(int)mesh.Nx+halo; q++) {

        indx = index((i+1)*mesh.Nx + q, (j+1)*mesh.Ny + r);

        /*
        std::cout << "setting " << q << " " << r << " " << s;
        std::cout << " to indx " << indx;
        std::cout << " jx " <<  jx[ indx ][0];
        std::cout << " jy " <<  jy[ indx ][0];
        std::cout << " jz " <<  jz[ indx ][0];
        std::cout << " \n";
        */

        mesh.jx1(q,r,s) = jx[ indx ][0];
        mesh.jy1(q,r,s) = jy[ indx ][0];
        mesh.jz1(q,r,s) = jz[ indx ][0];
      }
    }

  }


  // --------------------------------------------------
  // auxiliary/utility functions for debugging

  void set_kernel(std::vector<Realf>& in)
  {
    if(in.size() != (size_t)Nx*Ny*Nz) std::cout << "error in size!\n";

    int q = 0;
    for(int j = 0 ; j < Ny ; ++j)  
    for(int i = 0 ; i < Nx ; ++i, q++)  
      kernel[ index(i,j) ][0] = in[q];
  }


  void set_image(std::vector<Realf>& in)
  {
    if(in.size() != (size_t)Nx*Ny*Nz) std::cout << "error in size!\n";

    int q = 0;
    for(int j = 0 ; j < Ny ; ++j)  
    for(int i = 0 ; i < Nx ; ++i, q++)  
      jx[ index(i,j) ][0] = in[q];
  }


  std::vector<Realf> get_kernel()
  {
    std::vector<Realf> ret;

    for(int j = 0 ; j < Ny ; ++j)  
    for(int i = 0 ; i < Nx ; ++i)  
      ret.push_back( kernel[ index(i,j) ][0] );

    return ret;
  }


  std::vector<Realf> get_image()
  {
    std::vector<Realf> ret;

    for(int j = 0 ; j < Ny ; ++j)  
    for(int i = 0 ; i < Nx ; ++i)  
      ret.push_back( jx[ index(i,j) ][0] );

    return ret;
  }




};

} // end of namespace pic
