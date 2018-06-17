#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <assert.h>

#include "cell.h"
#include "communicate.h"

#include "../units.h"
#include "fftw3.h"



namespace pic {

using units::pi;

//! Spatial current filter using fftw
//
// Because of the fftw3 usage this class is internally heavily relying on C code
// instead of C++.
class Filter {
    
  /// circular/wrap indexing
  inline int circular(int x, int M)
  {
    if (x<0)
      return x+M;
    if(x >= M)
      return x-M;
    return x;
  }

  /// internal circular indexing, 
  // NOTE: this is different from the rest of the code but done like this
  // in order to appear consistent with the fttw examples.
  inline int index(int i, int j) {
    i = circular(i, height);
    j = circular(j, width);

    assert(i >= 0 && i < height);
    assert(j >= 0 && j < width);

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
    //int iff = (i - ceil(height/2) >= 0) ? i : height + i;
    //int jff = (j - ceil(width /2) >= 0) ? i : width  + j;
    int iff = i >= 0 ? i : height + i;
    int jff = j >= 0 ? j : width  + j;

    return std::make_tuple(iff, jff);
  }


  /// initialize kernel 
  //
  // In practice this is just a complex way to set kernel[0,0] = 1,
  // but we do it like this for clarity as the syntax is used later on.
  virtual void init_kernel() 
  {
    double val;

    // kernel size (even number because of initialization)
    int knx = height/3;
    int kny = width /3;

    // halo region size
    int h1 = height/2; // division is floor(x/y) automatically
    int w1 = width /2;
    int h2 = (height + 2 - 1)/2; // ceil(x/y)
    int w2 = (width  + 2 - 1)/2;

    int wi,wj;

    for(int i =-h1; i<=h1 ; ++i) {
      for(int j=-w1; j<=w1 ; ++j) {
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



  // in, out are m x n images (integer data)
  // K is the kernel size (KxK) - currently needs to be an odd number, e.g. 3
  // coeffs[K][K] is a 2D array of integer coefficients
  // scale is a scaling factor to normalise the filter gain
  virtual void init_3point_kernel(int times)
  {

    // kernel
    int K = 3; // three point kernel

    // 3-point digital filter
    std::vector<double> coeffs = {{ 1., 2., 1.,
                                    2., 4., 2.,
                                    1., 2., 1. }};

    // normalize
    // double scale = 0.0;
    // for (double c : coeffs) scale += c;

    // out array
    //double data;
    //std::vector<double> out;
    //out.resize(width*height);
      

    // create temporary Real number image array
    // NOTE: can not easily copy kernel into pure real part due to interleaved nature
    std::vector<double> image;
    image.resize(height*width); 

    for (int i=0; i < height; ++i)
      for (int j=0; j < width;  ++j) 
        image[ index(i,j) ] = kernel[ index(i,j) ][0];



    // perform convolution N times
    for(int N=0; N < times; N++) {

      /*
      //for (int i = K/2; i < height - K/2; ++i) // iterate through image
      for (int i=0; i < height; ++i) // iterate through circular image
      {
        //for (int j = K/2; j < width -K/2; ++j) // iterate through image
        for (int j=0; j < width; ++j) // iterate through circular image
        {
          double sum = 0.0; // sum will be the sum of input data * coeff terms
      
          // convolution of single point
          for (int ii = -K/2; ii <= K/2; ++ii) // iterate over kernel
          {
            for (int jj = -K/2; jj <= K/2; ++jj)
            {
              //int data = in[i + ii][j +jj];
              data = kernel[ index(i+ii, j+jj) ][0]; // real part
              //double coeff = coeffs[ii + K / 2][jj + K / 2];
              double coeff = coeffs[ (ii + K/2)*K + (jj + K/2) ];
      
              sum += data * coeff;
            }
          }

          //out[i][j] = sum / scale; // scale sum of convolution products and store in output
          out[ index(i,j) ] = sum/scale; // scale sum of convolution products and store in output
        }
      } // end of conv

      // copy (real part) back
      for (int i=0; i < height; ++i)
        for (int j=0; j < width;  ++j) 
          kernel[ index(i,j) ][0] = out[ index(i,j) ];
      */

      direct_convolve(image.data(), coeffs.data(), K);

    } // end of ntimes loop

    
    // normalize
    double norm = 0.0;
    for (int i=0; i < height; ++i)
      for (int j=0; j < width;  ++j) 
        norm += image[ index(i,j) ];

    // normalize and copy back
    for (int i=0; i < height; ++i)
      for (int j=0; j < width;  ++j) 
        kernel[ index(i,j) ][0] = image[ index(i,j) ]/norm;

  }
  


  /// Direct circular convolve two arrays
  // NOTE: assumes height x width for image size
  void direct_convolve(
      double* image,
      double* kernel,
      int K)
  {

    // out array
    double data;
    std::vector<double> out;
    out.resize(width*height);

    //for (int i = K/2; i < height - K/2; ++i) // iterate through image
    for (int i=0; i < height; ++i) // iterate through circular image
    {
      //for (int j = K/2; j < width -K/2; ++j) // iterate through image
      for (int j=0; j < width; ++j) // iterate through circular image
      {
        double sum = 0.0; // sum will be the sum of input data * coeff terms

        // convolution of single point
        for (int ii = -K/2; ii <= K/2; ++ii) // iterate over kernel
        {
          for (int jj = -K/2; jj <= K/2; ++jj)
          {
            //int data = in[i + ii][j +jj];
            data = image[ index(i+ii, j+jj) ]; 
            double coeff = kernel[ (ii + K/2)*K + (jj + K/2) ];

            sum += data * coeff;
          }
        }
        out[ index(i,j) ] = sum; // sum of convolution products and store in output
      }
    } // end of conv

    // copy (real part) back
    for (int i=0; i < height; ++i)
      for (int j=0; j < width;  ++j) 
        image[ index(i,j) ] = out[ index(i,j) ];

  } 




  /// initialize Gaussian kernel given the sigma
  /*
  virtual void init_gaussian_kernel(double sigma) 
  {
    double val;
    int h2 = floor(height/2);
    int w2 = floor(width /2);
    int wi, wj;

    double sum = 0.0;
    for(int i =-h2; i<=h2 ; ++i) {
      for(int j=-w2; j<=w2 ; ++j) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);

        val = 0.0;
        if (i*i + j*j < 200.0) {
          val = exp(- 0.5*((double)(i*i))/sigma/sigma)
              * exp(- 0.5*((double)(j*j))/sigma/sigma);
        }

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
  /*
  virtual void init_lowpass_fft_kernel(int cutoff) 
  {
    int h2min = (int)floor((double)height/2);
    int w2min = (int)floor((double)width /2);
    int h2max = (int)ceil( (double)height/2);
    int w2max = (int)ceil( (double)width /2);

    int h3 = height/3;
    int w3 = width /3;
    int wi, wj;

    std::cout << "h2: " << h2min << " " << h2max << " w2: " << w2min << " " << w2max << '\n';
    std::cout << "h3: " << h3 << " w3: " << w3 << '\n';

    // initialize to zero
    for(int i =-h2min; i<h2max; ++i) {
      for(int j=-w2min; j<w2max; ++j) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);

        kernel[ index(wi,wj) ][0] = 0.0; // real part
        kernel[ index(wi,wj) ][1] = 0.0; // complex part
      }
    }

    // make central bins 1.0
    double sum = 0.0;
    for(int i =-h2min+cutoff; i<h2max-cutoff; ++i) {
      for(int j=-w2min+cutoff; j<w2max-cutoff; ++j) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);

        kernel[ index(wi,wj) ][0] = 1.0; // real part
        //kernel[ index(wi,wj) ][1] = 0.0; // complex part
        sum += 1.0;
      }
    }

    // filter out lowest waves to avoid cyclic boundaries
    for(int i =-h3; i<h3; ++i) {
      for(int j=-w3; j<w3; ++j) {
        auto zindx = zero_wrapped_index(i,j);
        wi = std::get<0>(zindx);
        wj = std::get<1>(zindx);

        kernel[ index(wi,wj) ][0] = 0.0; // real part
        kernel[ index(wi,wj) ][1] = 0.0; // complex part
      }
    }

    // normalize 
    //for(int i =0; i<height ; ++i)
    //for(int j=0; j<width ; ++j)
    //  kernel[ index(i,j) ][0] /= sum;

  }
  */


  // normalize fft transformation
  void normalize() 
  {
    for(int i  = 0 ; i < height ; ++i)  
      for(int j = 0 ; j < width ; ++j)  
        jx[ index(i,j) ][0] /= width*height;

  }


  /// FFT kernel (once is enough)
  virtual void fft_kernel()         { fftw_execute(p_kernel); }

  /// FFT image forward
  virtual void fft_image_forward()  { fftw_execute(p_forw);   }

  /// FFT image backwards and normalize 
  virtual void fft_image_backward() { 
    fftw_execute(p_back);
    normalize();
  }


  /// Multiply kernel and image
  void apply_kernel()
  {
    double x, y, u, v;
    for(int i  = 0 ; i < height ; ++i) {
      for(int j = 0 ; j < width ; ++j) {

        x =     jx[ index(i,j) ][0];
        y =     jx[ index(i,j) ][1];
        u = kernel[ index(i,j) ][0];
        v = kernel[ index(i,j) ][1];

        jx[ index(i,j) ][0] = x*u - y*v;
        jx[ index(i,j) ][1] = x*v - y*u;
      }
    }


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
