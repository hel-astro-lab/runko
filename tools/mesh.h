#pragma once

#include <vector>
#include <stdexcept>


namespace toolbox {


/*! \brief simple dense 3D simulation mesh class
 *
 * Internally this is just a thin wrapper around STL vector class.
 *
 */

template <class T, int H=0> 
class Mesh {
  public:

    /// internal storage
    std::vector<T> mat;

    /// grid size along x
    size_t Nx;
      
    /// grid size along y
    size_t Ny;

    /// grid size along z
    size_t Nz;


    /// Internal indexing with halo region padding of width H
    size_t indx(int i, int j, int k) const {
      int indx = (i + H) + (Ny + 2*H)*( (j + H) + (Nz + 2*H)*(k + H));
      return indx;
    }

    /// standard (i,j,k) syntax
    T& operator()(int i, int j, int k) { 
      return mat[ indx(i,j,k) ];
    }
    const T& operator()(int i, int j, int k) const { 
      return mat[ indx(i,j,k) ];
    }


    /// Default initialization
    Mesh(size_t Nx, size_t Ny, size_t Nz) : Nx(Nx), Ny(Ny), Nz(Nz) {
      mat.resize( (Nx + 2*H)*(Ny + 2*H)*(Nz + 2*H) );
      std::fill(mat.begin(), mat.end(), T() ); // fill with zeros
    };

    /// clear internal storage (overriding with zeros to avoid garbage)
    void clear() {
      std::fill(mat.begin(), mat.end(), T() ); // fill with zeros
    }


    // Mesh arithmetics
    Mesh& operator=(const Mesh& rhs);

    Mesh& operator=(const T& rhs);

    Mesh& operator+=(const Mesh& rhs);

    Mesh& operator-=(const Mesh& rhs);

    Mesh& operator*=(const T& rhs);

    Mesh& operator/=(const T& rhs);


    /// Validate that two grids match
    void validateDims(const Mesh& rhs) {
      if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
      if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");
      if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
    };

    void copyVert(Mesh& rhs, int lhsJ, int rhsJ);

    void copyHorz(Mesh& rhs, int lhsI, int rhsI);

};


template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] = rhs.mat[i];
  }
  return *this;
}

template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator=(const T& rhs) {
  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] = rhs;
  }
  return *this;
}

template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator+=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] += rhs.mat[i];
  }
  return *this;
}

template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator-=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] -= rhs.mat[i];
  }
  return *this;
}

template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator*=(const T& rhs) {
  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] *= rhs;
  }
  return *this;
}

template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator/=(const T& rhs) {
  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] /= rhs;
  }
  return *this;
}


// Array arithmetics 
//-------------------------------------------------- 
template <class T, int H>
inline Mesh<T,H> operator+(Mesh<T,H> lhs, const Mesh<T,H>& rhs) {
  lhs += rhs;
  return lhs;
}

template <class T, int H>
inline Mesh<T,H> operator-(Mesh<T,H> lhs, const Mesh<T,H>& rhs) {
  lhs -= rhs;
  return lhs;
}

// Single value operators
//-------------------------------------------------- 
template <class T, int H>
inline Mesh<T,H> operator*(Mesh<T,H> lhs, const T& rhs) {
  lhs *= rhs;
  return lhs;
}

template <class T, int H>
inline Mesh<T,H> operator/(Mesh<T,H> lhs, const T& rhs) {
  lhs /= rhs;
  return lhs;
}




// Copy operators for different slices
//-------------------------------------------------- 

/// Copy vertical slice
// TODO: assumes implicitly 2D (x-y) arrays only by setting k=0 and ignoring it
template <class T, int H>
void Mesh<T,H>::copyVert(Mesh<T,H>& rhs, int lhsI, int rhsI) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int j=0; j<this->Ny; j++) { 
    this->operator()(lhsI, j, 0) = rhs(rhsI, j, 0);
  }

  /*
  for(int k=-H; k<this->Nz+H; k++) {
    for(int j=0; j<this->Ny; j++) { 
      this->operator()(lhsI, j, k) = rhs(rhsI, j, k);
    }
  }
  */

}


/// Copy horizontal slice 
// TODO: assumes implicitly 2D (x-y) arrays only by setting k=0 and ignoring it
template <class T, int H>
void Mesh<T,H>::copyHorz(Mesh<T,H>& rhs, int lhsJ, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");

  for(int i=0; i<this->Nx; i++) { 
    this->operator()(i, lhsJ, 0) = rhs(i, rhsJ, 0);
  }

  /*
  for(int k=-H; k<this->Nz+H; k++) {
    for(int i=0; i<this->Nx; i++) { 
      this->operator()(i, lhsJ, k) = rhs(i, rhsJ, k);
    }
  }
  */
}


} // end of namespace toolbox
