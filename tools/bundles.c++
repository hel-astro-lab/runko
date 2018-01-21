
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

#include <iostream>

#include "bundles.h"

using std::min;
using std::max;


using namespace toolbox;


/// Signum of value
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}



/// resize the container 
void Bundle::resize( size_t N) { 
    pencil.resize(N);
    grid.resize(N);
};

size_t Bundle::size() {
    return grid.size();
};

/// load bundle full of zeros
void Bundle::loadZeroBlock(size_t q) {
    pencil[q] = 0.0;
};

/// load block in x/y/z order 
void Bundle::loadBlock(size_t q, vblock_t block) {
    pencil[q] = block[0]; 
};

/// load values to the grid and transform the incoming cube according to dim
void Bundle::loadGrid(size_t q, Realf val) {
    grid[q] = val;
};

/// return the guiding grid
std::vector<Realf> Bundle::getGrid() {
    return grid;
};

/// return the pencil values
std::vector<Realf> Bundle::getPencil() {
    return pencil;
};

/// If current bundle slice is non-zero
bool Bundle::isNonZero(size_t q) {
    if ( pencil[q] == 0.0 ) { return false; };
    
    return true;
};

/// return q:th slice of the bundle
vblock_t Bundle::getSlice(size_t q) {
    vblock_t ret;
    ret[0] = pencil[q];

    return ret;
};

/// get grid size
Realf Bundle::getDx(size_t q) {
    return std::abs( grid[q+1] - grid[q] );
};




// --------------------------------------------------
void BundleInterpolator::setBundle(Bundle _bundle) {
  bundle = _bundle;
};

Bundle BundleInterpolator::getBundle( ) {
  return bundle;
};

void BundleInterpolator::setDelta( Bundle _delta ) {
  delta = _delta;
};

vblock_t BundleInterpolator::getDeltaSlice(size_t i) {
  return delta.getSlice(i);
};



Bundle BundleInterpolator2nd::interpolate( ) {

  // prepare target bundle
  Bundle ret;
  ret.resize( bundle.size() );

  // compute flux (inner region)
  vblock_t block, fp1, f0, Delta;

  ret.loadZeroBlock(0);
  for(size_t i=1; i<bundle.size()-1; i++) {
    fp1     = bundle.getSlice(i+1);
    f0      = bundle.getSlice(i  );

    // get shift 
    Delta     = getDeltaSlice(i);
    Delta[0] *= dt / bundle.getDx(i);

    // 2nd order conservative Lagrangian interpolation
    block[0] = Delta[0]          * ( fp1[0] + f0[0] )*0.5 
             - Delta[0]*Delta[0] * ( fp1[0] - f0[0] )*0.5;

    ret.loadBlock(i, block);
  }
  ret.loadZeroBlock( bundle.size()-1 );

  return ret;
};


Bundle BundleInterpolator4th::interpolate( ) {

  // prepare target bundle
  Bundle ret;
  ret.resize( bundle.size() );

  // compute flux (inner region)
  vblock_t block, fp2, fp1, f0, fm1, Delta;

  ret.loadZeroBlock(0);
  ret.loadZeroBlock(1);
  for(size_t i=2; i<bundle.size()-2; i++) {
    fm1     = bundle.getSlice(i-1);
    f0      = bundle.getSlice(i  );
    fp1     = bundle.getSlice(i+1);
    fp2     = bundle.getSlice(i+2);

    // get shift 
    Delta     = getDeltaSlice(i);
    Delta[0] *= dt / bundle.getDx(i);

    // flux limiter
    if (Delta[0] > 0.90) { Delta[0] = 0.90; }


    // 4th order conservative Lagrangian interpolation
    block[0] =     Delta[0]    * (-fp2[0] + 7.0*fp1[0] + 7.0*f0[0] - fm1[0] )/12.0
              +pow(Delta[0],2) * ( fp2[0] -15.0*fp1[0] +15.0*f0[0] - fm1[0] )/24.0
              +pow(Delta[0],3) * ( fp2[0] -     fp1[0] -     f0[0] + fm1[0] )/12.0
              +pow(Delta[0],4) * (-fp2[0] + 3.0*fp1[0] - 3.0*f0[0] + fm1[0] )/24.0;

    ret.loadBlock(i, block);
  }
  ret.loadZeroBlock( bundle.size()-2 );
  ret.loadZeroBlock( bundle.size()-1 );

  return ret;
};


Bundle BundleInterpolator4PIC::interpolate( ) {

  // prepare target bundle
  Bundle ret;
  ret.resize( bundle.size() );

  // compute flux (inner region)
  vblock_t block;
  Realf fp2, fp1, fp0, fm1, fm2;


  // temporary variables for loop
  Realf aa, as, aa1, aa2, aa3;
  int fa, ss, j1;

  Realf hmax1, hmax2, hmin1, hmin2;
  Realf hmax, hmin;
  Realf ep1, ep2, ep3;
  Realf flux;


  ret.loadZeroBlock(0);
  ret.loadZeroBlock(1);
  for(int i=2; i<(int)bundle.size()-2; i++) {

    // get shift (CFL number) and its sign (for detecting upwind)
    aa = getDeltaSlice(i)[0] * dt/bundle.getDx(i);
    fa = (int) floor(aa);
    ss = sign(aa);
    as = aa * ss;

    // compute how many grid indices do we advect
    j1 = i - fa;
    if (j1 < 2) {
      ret.loadZeroBlock(i);
      continue;
    }
    if (j1 > (int)bundle.size()-2) {
      ret.loadZeroBlock(i);
      continue;
    }

    // cfl shifts for interpolation formula
    aa1= -(as+1) * (as-1) * (as-2);
    aa2=2*(as+3) * (as-1) * (as-2);
    aa3= -(as+2) * (as+1) * (as-1);

    // elements (shifted for upwind)
    fm2     = bundle.getSlice(j1-ss-ss)[0];
    fm1     = bundle.getSlice(j1-ss   )[0];
    fp0     = bundle.getSlice(j1      )[0];
    fp1     = bundle.getSlice(j1+ss   )[0];
    fp2     = bundle.getSlice(j1+ss+ss)[0];

    hmax1 = max( max(fm1,fp0), min(fm1*2-fm2,fp0*2-fp1) );
    hmin1 = min( min(fm1,fp0), max(fm1*2-fm2,fp0*2-fp1) );
    hmax2 = max( max(fp1,fp0), min(fp0*2-fm1,fp1*2-fp2) );
    hmin2 = min( min(fp1,fp0), max(fp0*2-fm1,fp1*2-fp2) );
        
    hmax = max(hmax1,hmax2);
    hmin = max( (Realf)0.0, min(hmin1,hmin2));

    ep3=-min(fm1-fp0, (Realf) 2.4*(fp0-hmin))*(fp0<=fm1)
        +min(fp0-fm1, (Realf) 4.0*(hmax-fp0))*(fp0> fm1);
    ep2= min(fp1-fp0, (Realf) 1.6*(fp0-hmin))*(fp1>=fp0)
        -min(fp0-fp1, (Realf) 1.6*(hmax-fp0))*(fp1< fp0);
    ep1= min(fp2-fp1, (Realf) 2.0*(fp0-hmin))*(fp2>=fp1)*(ep2>=ep3) 
        -min(fp1-fp2, (Realf) 1.1*(fp0-hmin))*(fp2< fp1)*(ep2>=ep3) 
        +min(fp2-fp1, (Realf) 0.8*(hmax-fp0))*(fp2>=fp1)*(ep2< ep3) 
        -min(fp1-fp2, (Realf) 1.9*(hmax-fp0))*(fp2< fp1)*(ep2< ep3);

    flux = aa*((ep1*aa1 + ep2*aa2 + ep3*aa3)/24+fp0);
    block[0] = flux;

    ret.loadBlock(i, block);
  }
  ret.loadZeroBlock( bundle.size()-2 );
  ret.loadZeroBlock( bundle.size()-1 );

  return ret;
};
