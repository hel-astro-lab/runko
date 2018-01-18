#pragma once

#include <vector>
#include <array>


#include "../definitions.h"


namespace toolbox {

  // -------------------------------------------------- 
  /// Sheet of velocities from vmesh
  // These are always slices of the full mesh along some dimension at some 
  // location.
  class Sheet {

    public:

      /// guiding grid along horizontal dimensions of the sheet (i.e., x)
      std::vector<Realf> iGrid;

      // sheet size in horizontal dim
      size_t Ni = 0;

      /// guiding grid along vertical dimensions of the sheet (i.e., y)
      std::vector<Realf> jGrid;

      // sheet size in horizontal dim
      size_t Nj = 0;

      /// coordinate of the sheet in the slicing dimension
      double sliceVal;

      /// Value storage of the sheet
      std::vector<Realf> values;


      void resize(size_t Ni_, size_t Nj_);

      size_t getIndex(size_t i, size_t j);

      void loadValue(size_t i, size_t j, Realf val);

      void loadZeroBlock(size_t i, size_t j);

      void loadBlock(size_t i, size_t j, vblock_t block);

      vblock_t getBlock(size_t i, size_t j);

      bool isNonZero(size_t i, size_t j);

      void checkSizes(const Sheet& s);

      Realf sum();

      // Sheet diff();


      // Sheet arithmetics
      Sheet& operator+=(const Sheet& rhs);

      Sheet& operator-=(const Sheet& rhs);

      Sheet& operator*=(const Sheet& rhs);

      Sheet& operator/=(const Sheet& rhs);

      template<typename T>
      Sheet& operator*=(const T rhs);

      template<typename T>
      Sheet& operator/=(const T rhs);
  };

  inline Sheet operator+(Sheet lhs, const Sheet& rhs) {
    lhs += rhs;
    return lhs;
  };

  inline Sheet operator-(Sheet lhs, const Sheet& rhs) {
    lhs -= rhs;
    return lhs;
  };

  inline Sheet operator*(Sheet lhs, const Sheet& rhs) {
    lhs *= rhs;
    return lhs;
  };

  template<typename T>
  inline Sheet operator*(Sheet lhs, const T rhs) {
    lhs *= (Realf) rhs;
    return lhs;
  };

  template<typename T>
  inline Sheet operator*(const T rhs, Sheet lhs)  {
    lhs *= (Realf) rhs;
    return lhs;
  };

  template<typename T>
  inline Sheet operator/(Sheet lhs, const T rhs) {
    lhs /= (Realf) rhs;
    return lhs;
  };

  template<typename T>
  inline Sheet operator/(const T rhs, Sheet lhs)  {
    lhs /= (Realf) rhs;
    return lhs;
  };


};


// --------------------------------------------------
// Implementations


/// Resize the sheet into correct size
inline void toolbox::Sheet::resize(size_t Ni_, size_t Nj_) {
  iGrid.resize(Ni_);
  jGrid.resize(Nj_);
  values.resize(Ni_*Nj_);

  Ni = Ni_;
  Nj = Nj_;
}

/// internal function to get general id from sheet indices
inline size_t toolbox::Sheet::getIndex(size_t i, size_t j) {
  return Ni*j + i;
}

/// Load scalar to the sheet
inline void toolbox::Sheet::loadValue(size_t i, size_t j, Realf val) {
  size_t indx = getIndex(i, j);
  values[indx] = val;
}

/// load zeros to location (i,j)
inline void toolbox::Sheet::loadZeroBlock(size_t i, size_t j) {
  size_t indx = getIndex(i, j);

  // TODO add block instead of scalar
  values[indx] = 0.0;
}

inline void toolbox::Sheet::loadBlock(size_t i, size_t j, vblock_t block) {
  size_t indx = getIndex(i, j);

  // TODO add block instead of scalar
  values[indx] = block[0];
}

/// return block at location (i,j)
inline vblock_t toolbox::Sheet::getBlock(size_t i, size_t j) {
  vblock_t ret;
  size_t indx = getIndex(i, j);
  ret[0] = values[indx]; //TODO return block instead element

  return ret;
}

/// check if block at location (i,j) is zero
inline bool toolbox::Sheet::isNonZero(size_t i, size_t j) {
  size_t indx = getIndex(i, j);
  if ( values[indx] == 0.0 ) { return false; };
  return true;
}

/// differential volume elements of sheet cells
/* TODO implement; is this even needed?
Sheet Sheet::diff() {

  // initialize new return sheet
  Sheet ret;
  ret.resize(Ni, Nj);

  // compute lengths


  return ret;
}
*/

/// 
inline Realf toolbox::Sheet::sum() {
  Realf sum = 0.0;
  for(auto v : values) { sum += v; }

  return sum;
}

/// Check that another sheet conforms to my size & dimensions
inline void toolbox::Sheet::checkSizes(const Sheet& s) {
  if(this->Ni != s.Ni) throw std::range_error ("i dimensions do not match"); 
  if(this->Nj != s.Nj) throw std::range_error ("j dimensions do not match"); 
}

// Sheet arithmetics
inline toolbox::Sheet& toolbox::Sheet::operator+=(const Sheet& rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] += rhs.values[q];
  return *this;
}

inline toolbox::Sheet& toolbox::Sheet::operator-=(const Sheet& rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] -= rhs.values[q];
  return *this;
}

template<typename T>
inline toolbox::Sheet& toolbox::Sheet::operator*=(const T rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] *= (Realf) rhs;
  return *this;
}

template<typename T>
inline toolbox::Sheet& toolbox::Sheet::operator/=(const T rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] /= (Realf) rhs;
  return *this;
}

inline toolbox::Sheet& toolbox::Sheet::operator*=(const Sheet& rhs) {
  this->checkSizes(rhs);

  for(size_t q=0; q<(this->Ni*this->Nj); q++) { 
    this->values[q] *= rhs.values[q];
  }
  return *this;
}



