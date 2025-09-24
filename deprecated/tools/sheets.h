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
      std::vector<float_m> iGrid;

      // sheet size in horizontal dim
      size_t Ni = 0;

      /// guiding grid along vertical dimensions of the sheet (i.e., y)
      std::vector<float_m> jGrid;

      // sheet size in horizontal dim
      size_t Nj = 0;

      /// coordinate of the sheet in the slicing dimension
      double slice_value;

      /// Value storage of the sheet
      std::vector<float_m> values;


      void resize(size_t Ni_, size_t Nj_);

      size_t get_index(size_t i, size_t j);

      void load_value(size_t i, size_t j, float_m val);

      void load_zero_block(size_t i, size_t j);

      void load_block(size_t i, size_t j, vblock_t block);

      vblock_t get_block(size_t i, size_t j);

      bool is_non_zero(size_t i, size_t j);

      void check_size(const Sheet& s);

      float_m sum();

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
    lhs *= (float_m) rhs;
    return lhs;
  };

  template<typename T>
  inline Sheet operator*(const T rhs, Sheet lhs)  {
    lhs *= (float_m) rhs;
    return lhs;
  };

  template<typename T>
  inline Sheet operator/(Sheet lhs, const T rhs) {
    lhs /= (float_m) rhs;
    return lhs;
  };

  template<typename T>
  inline Sheet operator/(const T rhs, Sheet lhs)  {
    lhs /= (float_m) rhs;
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
inline size_t toolbox::Sheet::get_index(size_t i, size_t j) {
  return Ni*j + i;
}

/// Load scalar to the sheet
inline void toolbox::Sheet::load_value(size_t i, size_t j, float_m val) {
  size_t indx = get_index(i, j);
  values[indx] = val;
}

/// load zeros to location (i,j)
inline void toolbox::Sheet::load_zero_block(size_t i, size_t j) {
  size_t indx = get_index(i, j);

  // TODO add block instead of scalar
  values[indx] = 0.0;
}

inline void toolbox::Sheet::load_block(size_t i, size_t j, vblock_t block) {
  size_t indx = get_index(i, j);

  // TODO add block instead of scalar
  values[indx] = block[0];
}

/// return block at location (i,j)
inline vblock_t toolbox::Sheet::get_block(size_t i, size_t j) {
  vblock_t ret;
  size_t indx = get_index(i, j);
  ret[0] = values[indx]; //TODO return block instead element

  return ret;
}

/// check if block at location (i,j) is zero
inline bool toolbox::Sheet::is_non_zero(size_t i, size_t j) {
  size_t indx = get_index(i, j);
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
inline float_m toolbox::Sheet::sum() {
  float_m sum = 0.0;
  for(auto v : values) { sum += v; }

  return sum;
}

/// Check that another sheet conforms to my size & dimensions
inline void toolbox::Sheet::check_size(const Sheet& s) {
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
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] *= (float_m) rhs;
  return *this;
}

template<typename T>
inline toolbox::Sheet& toolbox::Sheet::operator/=(const T rhs) {
  for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] /= (float_m) rhs;
  return *this;
}

inline toolbox::Sheet& toolbox::Sheet::operator*=(const Sheet& rhs) {
  this->check_size(rhs);

  for(size_t q=0; q<(this->Ni*this->Nj); q++) { 
    this->values[q] *= rhs.values[q];
  }
  return *this;
}



