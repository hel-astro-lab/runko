#pragma once

#include <vector>
#include <stdexcept>
#include <cassert>


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
      assert( (i >= -H) && (i <  (int)Nx + H)  );
      assert( (j >= -H) && (j <  (int)Ny + H)  );
      assert( (k >= -H) && (k <  (int)Nz + H)  );

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


    /// serialize 3D data cube into 1D vector
    std::vector<T> serialize() const {
      std::vector<T> ret;
      ret.reserve(Nx*Ny*Nz);

      for(int k=0; k<int(Nz); k++)
      for(int j=0; j<int(Ny); j++)
      for(int i=0; i<int(Nx); i++)
        ret.push_back(mat[ indx(i,j,k) ] );

      return ret;
    }

    /// load 3D data cube from 1D serial vector
    // TODO: vec or vec& ?
    void load_serial(
        std::vector<T> vec, 
        size_t Nx, size_t Ny, size_t Nz
        ) {
      mat.resize( (Nx + 2*H)*(Ny + 2*H)*(Nz + 2*H) );

      int q = 0;
      for(int k=0; k<int(Nz); k++)
      for(int j=0; j<int(Ny); j++)
      for(int i=0; i<int(Nx); i++) 
      {
        mat[ indx(i,j,k) ] = vec[q];
        q++;
      }

    }




    // Mesh arithmetics
    //Mesh& operator=(const Mesh& rhs);

    template<int H2>
    Mesh& operator=(const Mesh<T, H2>& rhs);

    Mesh& operator=(const T& rhs);

    template<int H2>
    Mesh& operator+=(const Mesh<T, H2>& rhs);

    template<int H2>
    Mesh& operator-=(const Mesh<T, H2>& rhs);

    Mesh& operator*=(const T& rhs);

    Mesh& operator/=(const T& rhs);


    /// Validate that two grids match
    // TODO: unify index testing; use assert() ?
    template<int H2>
    void validateDims(const Mesh<T,H2>& rhs) {
      if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
      if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");
      if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
    };

    template<int H2>
    void copyVert(Mesh<T, H2>& rhs, int lhsJ, int rhsJ);

    template<int H2>
    void addVert(Mesh<T, H2>& rhs, int lhsJ, int rhsJ);


    template<int H2>
    void copyHorz(Mesh<T, H2>& rhs, int lhsI, int rhsI);

    template<int H2>
    void addHorz(Mesh<T, H2>& rhs, int lhsI, int rhsI);


    template<int H2>
    void copyFace(Mesh<T, H2>& rhs, int lhsK, int rhsK);

    template<int H2>
    void addFace(Mesh<T, H2>& rhs, int lhsK, int rhsK);


    template<int H2>
    void copyZdirPencil(Mesh<T, H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ);

    template<int H2>
    void addZdirPencil(Mesh<T, H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ);
};


/*
template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] = rhs.mat[i];
  }
  return *this;
}
*/


template<typename T, int H>
template <int H2>
Mesh<T,H>& Mesh<T,H>::operator=(const Mesh<T,H2>& rhs) {
  validateDims(rhs);

  for(size_t k=0;  k<this->Nz; k++) {
    for(size_t j=0;  j<this->Ny; j++) {
      for(size_t i=0;  i<this->Nx; i++) {
        this->operator()(i,j,k) = rhs(i,j,k);
      }
    }
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


/*
template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator+=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] += rhs.mat[i];
  }
  return *this;
}
*/

template<typename T, int H>
template <int H2>
Mesh<T,H>& Mesh<T,H>::operator+=(const Mesh<T,H2>& rhs) {
  validateDims(rhs);

  for(size_t k=0;  k<this->Nz; k++) {
    for(size_t j=0;  j<this->Ny; j++) {
      for(size_t i=0;  i<this->Nx; i++) {
        this->operator()(i,j,k) += rhs(i,j,k);
      }
    }
  }

  return *this;
}


/*
template <class T, int H>
Mesh<T,H>& Mesh<T,H>::operator-=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] -= rhs.mat[i];
  }
  return *this;
}
*/

template<typename T, int H>
template <int H2>
Mesh<T,H>& Mesh<T,H>::operator-=(const Mesh<T,H2>& rhs) {
  validateDims(rhs);

  for(size_t k=0;  k<this->Nz; k++) {
    for(size_t j=0;  j<this->Ny; j++) {
      for(size_t i=0;  i<this->Nx; i++) {
        this->operator()(i,j,k) -= rhs(i,j,k);
      }
    }
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
/*
template <class T, int H>
inline Mesh<T,H> operator+(Mesh<T,H> lhs, const Mesh<T,H>& rhs) {
  lhs += rhs;
  return lhs;
}
*/

template <class T, int H, int H2>
inline Mesh<T,H> operator+(Mesh<T,H> lhs, const Mesh<T,H2>& rhs) {
  lhs += rhs;
  return lhs;
}


/*
template <class T, int H>
inline Mesh<T,H> operator-(Mesh<T,H> lhs, const Mesh<T,H>& rhs) {
  lhs -= rhs;
  return lhs;
}
*/

template <class T, int H, int H2>
inline Mesh<T,H> operator-(Mesh<T,H> lhs, const Mesh<T,H2>& rhs) {
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
template <int H2>
void Mesh<T,H>::copyVert(Mesh<T,H2>& rhs, int lhsI, int rhsI) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int j=0; j<(int)this->Ny; j++) { 
    for(int k=0; k<(int)this->Nz; k++) {
      this->operator()(lhsI, j, k) = rhs(rhsI, j, k);
    }
  }
}

/// Add vertical slice
template <class T, int H>
template <int H2>
void Mesh<T,H>::addVert(Mesh<T,H2>& rhs, int lhsI, int rhsI) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int j=0; j<(int)this->Ny; j++) { 
    for(int k=0; k<(int)this->Nz; k++) {
      this->operator()(lhsI, j, k) += rhs(rhsI, j, k);
    }
  }
}


/// Copy horizontal slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::copyHorz(Mesh<T,H2>& rhs, int lhsJ, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");

  for(int i=0; i<(int)this->Nx; i++) { 
    for(int k=0; k<(int)this->Nz; k++) {
      this->operator()(i, lhsJ, k) = rhs(i, rhsJ, k);
    }
  }
}
  
/// Add horizontal slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::addHorz(Mesh<T,H2>& rhs, int lhsJ, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");

  for(int i=0; i<(int)this->Nx; i++) { 
    for(int k=0; k<(int)this->Nz; k++) {
      this->operator()(i, lhsJ, k) += rhs(i, rhsJ, k);
    }
  }
}


/// Copy face slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::copyFace(Mesh<T,H2>& rhs, int lhsK, int rhsK) {
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int i=0; i<(int)this->Nx; i++) { 
    for(int j=0; j<(int)this->Ny; j++) {
      this->operator()(i, j, lhsK) = rhs(i, j, rhsK);
    }
  }
}
  
/// Add face slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::addFace(Mesh<T,H2>& rhs, int lhsK, int rhsK) {
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int i=0; i<(int)this->Nx; i++) { 
    for(int j=0; j<(int)this->Ny; j++) {
      this->operator()(i, j, lhsK) += rhs(i, j, rhsK);
    }
  }
}

// copy pencil pointing along Z
template <class T, int H>
template <int H2>
void Mesh<T,H>::copyZdirPencil(Mesh<T,H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) { 
    this->operator()(lhsI, lhsJ, k) = rhs(rhsI, rhsJ, k);
  }
}

// add pencil pointing along Z
template <class T, int H>
template <int H2>
void Mesh<T,H>::addZdirPencil(Mesh<T,H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) { 
    this->operator()(lhsI, lhsJ, k) += rhs(rhsI, rhsJ, k);
  }
}



} // end of namespace toolbox
