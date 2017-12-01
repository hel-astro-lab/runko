#pragma once

#include <vector>
#include <stdexcept>


namespace toolbox {


/*! \brief simple 3D dense matrix class
 *
 * Internally this is just a thin wrapper around STL vector class.
 *
 */

template <class T> 
class DenseGrid {
  public:

    /// internal storage
    std::vector<T> mat;

    /// grid size along x
    size_t Nx;
      
    /// grid size along y
    size_t Ny;

    /// grid size along z
    size_t Nz;

    /// Internal indexing
    size_t indx(size_t i, size_t j, size_t k) {
      size_t indx = i + Ny*(j + Nz*k);
      return indx;
    }

    /// standard (i,j,k) syntax
    T& operator()(size_t i, size_t j, size_t k) { 
      return mat[ indx(i,j,k) ];
    }
    const T& operator()(size_t i, size_t j, size_t k) const { 
      return mat[ indx(i,j,k) ];
    }


    /// Default initialization
    DenseGrid(size_t Nx, size_t Ny, size_t Nz) : Nx(Nx), 
                                                 Ny(Ny),
                                                 Nz(Nz) {
      mat.resize(Nx*Ny*Nz);
      std::fill(mat.begin(), mat.end(), T() ); // fill with zeros
    };

    /// clear internal storage (overriding with zeros to avoid garbage)
    void clear() {
      std::fill(mat.begin(), mat.end(), T() ); // fill with zeros
    }


    // DenseGrid arithmetics
    DenseGrid& operator=(const DenseGrid& rhs);

    DenseGrid& operator=(const T& rhs);

    DenseGrid& operator+=(const DenseGrid& rhs);

    DenseGrid& operator-=(const DenseGrid& rhs);

    DenseGrid& operator*=(const T& rhs);

    DenseGrid& operator/=(const T& rhs);


    /// Validate that two grids match
    void validateDims(const DenseGrid<T>& rhs) {
      if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
      if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");
      if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
    };

};


template <class T>
DenseGrid<T>& DenseGrid<T>::operator=(const DenseGrid<T>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] = rhs.mat[i];
  }
  return *this;
}

template <class T>
DenseGrid<T>& DenseGrid<T>::operator=(const T& rhs) {
  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] = rhs;
  }
  return *this;
}

template <class T>
DenseGrid<T>& DenseGrid<T>::operator+=(const DenseGrid<T>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] += rhs.mat[i];
  }
  return *this;
}

template <class T>
DenseGrid<T>& DenseGrid<T>::operator-=(const DenseGrid<T>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] -= rhs.mat[i];
  }
  return *this;
}

template <class T>
DenseGrid<T>& DenseGrid<T>::operator*=(const T& rhs) {
  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] *= rhs;
  }
  return *this;
}

template <class T>
DenseGrid<T>& DenseGrid<T>::operator/=(const T& rhs) {
  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] /= rhs;
  }
  return *this;
}


// Array arithmetics 
//-------------------------------------------------- 
template <class T>
inline DenseGrid<T> operator+(DenseGrid<T> lhs, const DenseGrid<T>& rhs) {
  lhs += rhs;
  return lhs;
}

template <class T>
inline DenseGrid<T> operator-(DenseGrid<T> lhs, const DenseGrid<T>& rhs) {
  lhs -= rhs;
  return lhs;
}

// Single value operators
//-------------------------------------------------- 
template <class T>
inline DenseGrid<T> operator*(DenseGrid<T> lhs, const T& rhs) {
  lhs *= rhs;
  return lhs;
}

template <class T>
inline DenseGrid<T> operator/(DenseGrid<T> lhs, const T& rhs) {
  lhs /= rhs;
  return lhs;
}




} // end of namespace toolbox
