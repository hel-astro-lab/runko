#pragma once

#include <iostream>

#include<algorithm>
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
    inline size_t indx(int i, int j, int k) const {
      assert( (i >= -H) && (i <  (int)Nx + H)  );
      assert( (j >= -H) && (j <  (int)Ny + H)  );
      assert( (k >= -H) && (k <  (int)Nz + H)  );


      int indx = (i + H) + (Nx + 2*H)*( (j + H) + (Ny + 2*H)*(k + H));

      /*
      std::cout << "indx: " << indx;
      std::cout << "i: " << i << " j " << j << " k " << k << "\n";
      std::cout << "H: " << H << "\n";
      std::cout << "size: " << mat.size() << "\n";
      std::cout << "Nx: " << Nx << " Ny " << Ny << " Nz " << Nz << "\n";
      */

      assert( (indx >= 0) && (indx <  (int)mat.size() ) );

      return indx;
    }

    /// standard (i,j,k) syntax
    T& operator()(int i, int j, int k) { 
      return mat[ indx(i,j,k) ];
    }

    const T& operator()(int i, int j, int k) const { 
      return mat[ indx(i,j,k) ];
    }

    /// empty default constructor
    //Mesh() {};

    /// Default initialization
    Mesh(size_t Nx_in, size_t Ny_in, size_t Nz_in) : 
      Nx(Nx_in), Ny(Ny_in), Nz(Nz_in) 
    {
      mat.resize( (Nx + 2*H)*(Ny + 2*H)*(Nz + 2*H) );
      std::fill(mat.begin(), mat.end(), T() ); // fill with zeros
    };

    // 2D shortcut
    Mesh(size_t Nx_in, size_t Ny_in) :
      Nx(Nx_in), Ny(Ny_in), Nz(1) { Mesh(Nx, Ny, Nz); }

    // 1D shortcut
    Mesh(size_t Nx_in) : 
      Nx(Nx_in), Ny(1), Nz(1)     { Mesh(Nx, Ny, Nz); }


    // explicit default copy operator
    Mesh(Mesh& other) = default;
    Mesh(const Mesh& other) = default;
    
    //Mesh& operator=(const Mesh& other) = default;

    /// address to data
    T* data() { return mat.data(); }

    const T* data() const {return mat.data(); }

    /// internal storage size
    size_t size() { return mat.size(); }

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
    void unserialize(
        std::vector<T>& vec, 
        size_t Nx_in, size_t Ny_in, size_t Nz_in
        ) {

      Nx = Nx_in;
      Ny = Ny_in;
      Nz = Nz_in;
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
    //=
    Mesh& operator=(const Mesh<T, H>& rhs);

    template<int H2>
    Mesh& operator=(const Mesh<T, H2>& rhs);

    Mesh& operator=(const T& rhs);


    //+=
    Mesh& operator+=(const Mesh<T, H>& rhs);

    template<int H2>
    Mesh& operator+=(const Mesh<T, H2>& rhs);

    //-=
    Mesh& operator-=(const Mesh<T, H>& rhs);

    template<int H2>
    Mesh& operator-=(const Mesh<T, H2>& rhs);

    //*=
    Mesh& operator*=(const T& rhs);

    // /=
    Mesh& operator/=(const T& rhs);


    /// Validate that two grids match
    // TODO: unify index testing; use assert() ?
    template<int H2>
    void validateDims(const Mesh<T,H2>& rhs) {
      if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
      if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");
      if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
      //if(this->mat.size() != rhs.mat.size()) throw std::range_error ("container sizes do not match");
    }

    template<int H2>
    void copy_vert(Mesh<T, H2>& rhs, int lhsI, int rhsI);

    template<int H2>
    void add_vert(Mesh<T, H2>& rhs, int lhsI, int rhsI);


    template<int H2>
    void copy_horz(Mesh<T, H2>& rhs, int lhsJ, int rhsJ);

    template<int H2>
    void add_horz(Mesh<T, H2>& rhs, int lhsJ, int rhsJ);


    template<int H2>
    void copy_face(Mesh<T, H2>& rhs, int lhsK, int rhsK);

    template<int H2>
    void add_face(Mesh<T, H2>& rhs, int lhsK, int rhsK);


    template<int H2>
    void copy_z_pencil(Mesh<T, H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ);

    template<int H2>
    void add_z_pencil(Mesh<T, H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ);
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

/// = with any halo size
template<typename T, int H>
Mesh<T,H>& Mesh<T,H>::operator=(const Mesh<T,H>& rhs) {
  validateDims(rhs);
    
  // explicit deep copy
  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] = rhs.mat[i];
  }
    
  // implicit deep copy
  //this->mat = rhs.mat;
  return *this;
}


/// = with differing halo size
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


template<typename T, int H>
Mesh<T,H>& Mesh<T,H>::operator+=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] += rhs.mat[i];
  }

  // TODO: do not operate on halo regions
  //for(size_t k=0;  k<this->Nz; k++) {
  //  for(size_t j=0;  j<this->Ny; j++) {
  //    for(size_t i=0;  i<this->Nx; i++) {
  //      this->operator()(i,j,k) += rhs(i,j,k);
  //    }
  //  }
  //}

  return *this;
}

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


template<typename T, int H>
Mesh<T,H>& Mesh<T,H>::operator-=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  for(size_t i=0; i<this->mat.size(); i++) {
    this->mat[i] -= rhs.mat[i];
  }

  // TODO: do not operate on halo regions
  //for(size_t k=0;  k<this->Nz; k++) {
  //  for(size_t j=0;  j<this->Ny; j++) {
  //    for(size_t i=0;  i<this->Nx; i++) {
  //      this->operator()(i,j,k) -= rhs(i,j,k);
  //    }
  //  }
  //}
  return *this;
}

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
void Mesh<T,H>::copy_vert(Mesh<T,H2>& rhs, int lhsI, int rhsI) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) {
    for(int j=0; j<(int)this->Ny; j++) { 
      this->operator()(lhsI, j, k) = rhs(rhsI, j, k);
    }
  }
}

/// Add vertical slice
template <class T, int H>
template <int H2>
void Mesh<T,H>::add_vert(Mesh<T,H2>& rhs, int lhsI, int rhsI) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) {
    for(int j=0; j<(int)this->Ny; j++) { 
      this->operator()(lhsI, j, k) += rhs(rhsI, j, k);
    }
  }
}


/// Copy horizontal slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::copy_horz(Mesh<T,H2>& rhs, int lhsJ, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) {
    for(int i=0; i<(int)this->Nx; i++) { 
      this->operator()(i, lhsJ, k) = rhs(i, rhsJ, k);
    }
  }
}
  
/// Add horizontal slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::add_horz(Mesh<T,H2>& rhs, int lhsJ, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) {
    for(int i=0; i<(int)this->Nx; i++) { 
      this->operator()(i, lhsJ, k) += rhs(i, rhsJ, k);
    }
  }
}


/// Copy face slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::copy_face(Mesh<T,H2>& rhs, int lhsK, int rhsK) {
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int j=0; j<(int)this->Ny; j++) {
    for(int i=0; i<(int)this->Nx; i++) { 
      this->operator()(i, j, lhsK) = rhs(i, j, rhsK);
    }
  }
}
  
/// Add face slice 
template <class T, int H>
template <int H2>
void Mesh<T,H>::add_face(Mesh<T,H2>& rhs, int lhsK, int rhsK) {
  if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
  if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");

  for(int j=0; j<(int)this->Ny; j++) {
    for(int i=0; i<(int)this->Nx; i++) { 
      this->operator()(i, j, lhsK) += rhs(i, j, rhsK);
    }
  }
}

// copy pencil pointing along Z
template <class T, int H>
template <int H2>
void Mesh<T,H>::copy_z_pencil(Mesh<T,H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) { 
    this->operator()(lhsI, lhsJ, k) = rhs(rhsI, rhsJ, k);
  }
}

// add pencil pointing along Z
template <class T, int H>
template <int H2>
void Mesh<T,H>::add_z_pencil(Mesh<T,H2>& rhs, int lhsI, int lhsJ, int rhsI, int rhsJ) {
  if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");

  for(int k=0; k<(int)this->Nz; k++) { 
    this->operator()(lhsI, lhsJ, k) += rhs(rhsI, rhsJ, k);
  }
}



} // end of namespace toolbox
