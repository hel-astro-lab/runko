
#include <array>
#include <vector>
#include <cmath>

#include "bundles.h"


using namespace bundles;



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

    // 4th order conservative Lagrangian interpolation
    block[0] = Delta[0]        * (-fp2[0] + 7.0*fp1[0] + 7.0*f0[0] - fm1[0] )/12.0
              +pow(Delta[0],2) * ( fp2[0] -15.0*fp1[0] +15.0*f0[0] - fm1[0] )/24.0
              +pow(Delta[0],3) * ( fp2[0] -     fp1[0] -     f0[0] + fm1[0] )/12.0
              +pow(Delta[0],4) * (-fp2[0] + 3.0*fp1[0] - 3.0*f0[0] + fm1[0] )/24.0;

    ret.loadBlock(i, block);
  }
  ret.loadZeroBlock( bundle.size()-2 );
  ret.loadZeroBlock( bundle.size()-1 );

  return ret;
};


