#pragma once

#include <iostream>

#include<algorithm>
//#include <vector>
#include <stdexcept>
#include <cassert>
#include <exception>

#include "iter/allocator.h"
#include "iter/devcall.h"

#ifdef GPU
#include <cuda_runtime_api.h>
#else
#endif


namespace toolbox {


/*! \brief simple dense 3D simulation mesh class
 *
 * Internally this is just a thin wrapper around STL vector class.
 *
 */

template <typename T, int H> 
class Mesh 
{

  private:
    /// internal storage
    T *ptr;
    bool allocated{false};
    int count{0};
  public:

    /// grid size along x
    int Nx{0};
      
    /// grid size along y
    int Ny{0};

    /// grid size along z
    int Nz{0};

    /// Internal indexing with halo region padding of width H
    DEVCALLABLE
    inline size_t indx(int i, int j, int k) const {
      //if(not( (i >= -H) && (i <  (int)Nx + H)  )) std::cout << "XXX(" << i <<","<< j <<","<< k <<")\n";
      //if(not( (j >= -H) && (j <  (int)Ny + H)  )) std::cout << "XXX(" << i <<","<< j <<","<< k <<")\n";
      //if(not( (k >= -H) && (k <  (int)Nz + H)  )) std::cout << "XXX(" << i <<","<< j <<","<< k <<")\n";

      //assert( (i >= -H) && (i < (int)Nx + H)  );
      //assert( (j >= -H) && (j < (int)Ny + H)  );
      //assert( (k >= -H) && (k < (int)Nz + H)  );

      bool inx = (i >= -H) && (i < (int)Nx + H);
      bool iny = (j >= -H) && (j < (int)Ny + H);
      bool inz = (k >= -H) && (k < (int)Nz + H);

      if( !inx || !iny || !inz) {
          std::cerr << "MESH OUTSIDE TILE " << std::endl;
          std::cerr << i << " /" << Nx << " +H " << H << std::endl;
          std::cerr << j << " /" << Ny << " +H " << H << std::endl;
          std::cerr << k << " /" << Nz << " +H " << H << std::endl;

          if(i < -H) i = -H;
          if(j < -H) j = -H;
          if(k < -H) k = -H;

          if(i >= (int)Nx + H) i = Nx + H - 1;
          if(j >= (int)Ny + H) j = Ny + H - 1;
          if(k >= (int)Nz + H) k = Nz + H - 1;

          assert(false);
          //abort();
      }

      // this is true if above is true
      //int indx = (i + H) + (Nx + 2*H)*( (j + H) + (Ny + 2*H)*(k + H));
      //assert( (indx >= 0) && (indx < (int)count ) );

      assert(allocated);

      //return indx;
      return i + H + (Nx + 2*H)*( (j + H) + (Ny + 2*H)*(k + H));
    }

    /// 1D index 
    DEVCALLABLE
    inline T& operator()(size_t ind) { 
      return ptr[ ind ];
    }

    DEVCALLABLE
    inline const T& operator()(size_t ind) const { 
      return ptr[ ind ];
    }

    /// standard (i,j,k) syntax
    DEVCALLABLE
    inline T& operator()(int i, int j, int k) { 
      size_t ind = indx(i,j,k);
      assert(ind < count);
      return ptr[ind];
    }

    DEVCALLABLE
    inline const T& operator()(int i, int j, int k) const { 
      size_t ind = indx(i,j,k);
      assert(ind < count);
      return ptr[ind];
    }

    /// empty default constructor
    //Mesh() = default;
    Mesh() :
    allocated(false), count(0)
    { }


    /// standard initialization
    Mesh(int Nx, int Ny, int Nz) : 
      Nx(Nx), 
      Ny(Ny), 
      Nz(Nz),
      allocated(false), count(0)
      //mat( (Nx + 2*H)*(Ny + 2*H)*(Nz + 2*H) )
    {
      alloc( (Nx + 2*H)*(Ny + 2*H)*(Nz + 2*H));
      try {
        //mat.resize( (Nx + 2*H)*(Ny + 2*H)*(Nz + 2*H) ); //automatically done at construction
        //std::fill(ptr, ptr+count, T() ); // fill with zeros
        clear();
      } catch ( std::exception& e) {
        // whoops... if control reaches here, a memory allocation
        // failure occurred somewhere.
        std::cerr << "Standard exception: " << e.what() << std::endl;
        assert(false);
      }

      //if(mat.empty()) assert(false);
    };

    // explicit default copy operator
    //Mesh(Mesh& other) = default;
    Mesh(Mesh& other) :
      Nx(other.Nx),
      Ny(other.Ny),
      Nz(other.Nz),
      allocated(false), count(0)
      //mat(other.mat)
    { 
      Nx = other.Nx; 
      Ny = other.Ny; 
      Nz = other.Nz; 
      alloc(other.size());
      //mat.resize(other.mat.size());
      for(size_t i=0; i<other.size(); i++) ptr[i] = other.ptr[i];
    }

    // Mesh(const Mesh& other) = default;
    Mesh(const Mesh& other) :
      Nx(other.Nx),
      Ny(other.Ny),
      Nz(other.Nz),
      allocated(false), count(0)
      //mat(other.mat)
    { 
      Nx = other.Nx; 
      Ny = other.Ny; 
      Nz = other.Nz; 
      alloc(other.size());
      for(size_t i=0; i<other.size(); i++) ptr[i] = other.ptr[i];
    }
    
    // public swap for efficient memory management
    friend void swap(Mesh& first, Mesh& second)
    {
        using std::swap;
        swap(first.Nx, second.Nx);
        swap(first.Ny, second.Ny);
        swap(first.Nz, second.Nz);
        swap(first.ptr, second.ptr);
        swap(first.count, second.count);
        swap(first.allocated, second.allocated);
    }

    //Mesh& operator=(const Mesh& other) = default;
    // copy-and-swap algorithm
    //
    // NOTE: rhs is passed by value.
    // See: 
    // https://web.archive.org/web/20140113221447/http://cpp-next.com/archive/2009/08/want-speed-pass-by-value/
    // https://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
    Mesh& operator=(Mesh other) 
    {
      swap(*this, other); 
      return *this;
    }

    // move constructor
    Mesh(Mesh&& other) noexcept
        : Mesh() // initialize via default constructor, C++11 only
    {
        swap(*this, other);
    }

    ~Mesh()
    {
      // todo fix this 
      if(allocated) UniAllocator::deallocate(ptr);

      allocated = false;
      count = 0;
      
    }

    /// address to data
    DEVCALLABLE
    T* data() { return ptr; }
    
    DEVCALLABLE
    const T* data() const {return ptr; }

    /// internal storage size
    size_t size() const { return count; }

    /// clear internal storage (overriding with zeros to avoid garbage)
    void clear() {
      #ifdef GPU
        cudaMemset ( ptr, 0, count*sizeof(T) );
      #else
        std::fill(ptr, ptr+count, T() ); // fill with zeros
      #endif
      
    }

    /// fill halos with zeros
    void clear_halos() {
        for(int k=-H;  k<this->Nz+H; k++) {
        for(int j=-H;  j<this->Ny+H; j++) {
        for(int i=-H;  i<this->Nx+H; i++) {

            if(
                (i >= 0 && i < this->Nx) &&  
                (j >= 0 && j < this->Ny) &&  
                (k >= 0 && k < this->Nz) 
              ) { continue; }

            ptr[ indx(i,j,k) ] = 0.0;
        }}}

    }


    /// serialize 3D data cube into 1D vector
    std::vector<T> serialize() const {
      std::vector<T> ret;
      ret.reserve(Nx*Ny*Nz);

      for(int k=0; k<int(Nz); k++)
      for(int j=0; j<int(Ny); j++)
      for(int i=0; i<int(Nx); i++)
        ret.push_back(ptr[ indx(i,j,k) ] );

      return ret;
    }

    /// load 3D data cube from 1D serial vector
    void unserialize(
        std::vector<T>& vec, 
        int Nx_in, int Ny_in, int Nz_in
        ) {

      Nx = Nx_in;
      Ny = Ny_in;
      Nz = Nz_in;
      alloc( (Nx + 2*H)*(Ny + 2*H)*(Nz + 2*H) );

      int q = 0;
      for(int k=0; k<int(Nz); k++)
      for(int j=0; j<int(Ny); j++)
      for(int i=0; i<int(Nx); i++) 
      {
        ptr[ indx(i,j,k) ] = vec[q];
        q++;
      }
    }

    void alloc(int count_){
        if(allocated) UniAllocator::deallocate(ptr);

        ptr = UniAllocator::allocate<T>(count_);
        allocated = true;
        count = count_;

		return;
      }

    // Mesh arithmetics
    //=
    //Mesh& operator=(const Mesh<T, H>& rhs);

    template<int H2>
    Mesh& operator=(const Mesh<T, H2>& rhs);


    // scalar assignment
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
      //if(this->Nx != rhs.Nx) throw std::range_error ("x dimensions do not match");
      //if(this->Ny != rhs.Ny) throw std::range_error ("y dimensions do not match");
      //if(this->Nz != rhs.Nz) throw std::range_error ("z dimensions do not match");
      //if(this->mat.size() != rhs.mat.size()) throw std::range_error ("container sizes do not match");
      assert(this->Nx == rhs.Nx);
      assert(this->Ny == rhs.Ny);
      assert(this->Nz == rhs.Nz);
    }
};





//--------------------------------------------------

/// = with differing halo size
template<typename T, int H>
template <int H2>
inline Mesh<T,H>& Mesh<T,H>::operator=(const Mesh<T,H2>& rhs) {
  validateDims(rhs);
  //for(size_t i=0; i<this->mat.size(); i++) this->mat[i] = rhs.mat[i];

  for(int k=0;  k<this->Nz; k++) {
    for(int j=0;  j<this->Ny; j++) {
      for(int i=0;  i<this->Nx; i++) {
        this->operator()(i,j,k) = rhs(i,j,k);
      }
    }
  }

  return *this;
}

/// + with same halo size
template <class T, int H>
inline Mesh<T,H>& Mesh<T,H>::operator=(const T& rhs) {
  // overwriting internal container with a scalar
  for(size_t i=0; i<this->size(); i++) {
    this->ptr[i] = rhs;
  }
  return *this;
}


template<typename T, int H>
inline Mesh<T,H>& Mesh<T,H>::operator+=(const Mesh<T,H>& rhs) {
  validateDims(rhs);
  for(size_t i=0; i<this->size(); i++) this->ptr[i] += rhs.ptr[i];

  // TODO: do not operate on halo regions
  //for(int k=0;  k<this->Nz; k++) {
  //  for(int j=0;  j<this->Ny; j++) {
  //    for(int i=0;  i<this->Nx; i++) {
  //      this->operator()(i,j,k) += rhs(i,j,k);
  //    }
  //  }
  //}

  return *this;
}

template<typename T, int H>
template <int H2>
inline Mesh<T,H>& Mesh<T,H>::operator+=(const Mesh<T,H2>& rhs) {
  validateDims(rhs);
  //for(size_t i=0; i<this->mat.size(); i++) this->mat[i] += rhs.mat[i];

  for(int k=0;  k<this->Nz; k++) {
    for(int j=0;  j<this->Ny; j++) {
      for(int i=0;  i<this->Nx; i++) {
        this->operator()(i,j,k) += rhs(i,j,k);
      }
    }
  }

  return *this;
}


template<typename T, int H>
inline Mesh<T,H>& Mesh<T,H>::operator-=(const Mesh<T,H>& rhs) {
  validateDims(rhs);

  // purely vectorized version
  for(size_t i=0; i<this->size(); i++) this->ptr[i] -= rhs.ptr[i];

  // Version that does not operate on halo regions
  // this is more correct but every so slightly slower
  // because vectorization is interrupted.
  //for(int k=0;  k<this->Nz; k++) {
  //  for(int j=0;  j<this->Ny; j++) {
  //    for(int i=0;  i<this->Nx; i++) {
  //      this->operator()(i,j,k) -= rhs(i,j,k);
  //    }
  //  }
  //}
  return *this;
}

// -= for differing halos
template<typename T, int H>
template <int H2>
inline Mesh<T,H>& Mesh<T,H>::operator-=(const Mesh<T,H2>& rhs) {
  validateDims(rhs);
  //for(size_t i=0; i<this->mat.size(); i++) this->mat[i] -= rhs.mat[i];

  for(int k=0;  k<this->Nz; k++) {
    for(int j=0;  j<this->Ny; j++) {
      for(int i=0;  i<this->Nx; i++) {
        this->operator()(i,j,k) -= rhs(i,j,k);
      }
    }
  }

  return *this;
}

template <class T, int H>
inline Mesh<T,H>& Mesh<T,H>::operator*=(const T& rhs) {
  for(size_t i=0; i<this->size(); i++) {
    this->ptr[i] *= rhs;
  }
  return *this;
}

template <class T, int H>
inline Mesh<T,H>& Mesh<T,H>::operator/=(const T& rhs) {
  for(size_t i=0; i<this->size(); i++) {
    this->ptr[i] /= rhs;
  }
  return *this;
}


// Array arithmetics 
//-------------------------------------------------- 
template <class T, int H, int H2>
inline Mesh<T,H> operator+(Mesh<T,H> lhs, const Mesh<T,H2>& rhs) {
  lhs += rhs;
  return lhs;
}

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




} // end of namespace toolbox
