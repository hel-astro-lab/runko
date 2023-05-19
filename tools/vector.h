#pragma once

#include <iostream>
#include <array>
#include <cmath>

namespace toolbox {

  using std::ostream;

template <typename T> 
class Vec3 
{
  public:

    std::array<T, 3> data;

    inline T& operator()(size_t ind) { 
      return data[ ind ];
    }

    inline const T& operator()(size_t ind) const { 
      return data[ ind ];
    }

    Vec3() {};

    Vec3(T x, T y, T z) 
    {
      data[0] = x;  
      data[1] = y;  
      data[2] = z;  
    }

  friend ostream& operator<<(ostream& os, Vec3<T>& v)
  {
    os << "[";
    for(size_t i=0; i < 3; i++) os << v(i) << " ";
    os << "]";
    return os; 
  }

}; // end of Vec3 class


template <typename T> 
class Vec4 
{
  public:

    std::array<T, 4> data;

    inline T& operator()(size_t ind) { 
      return data[ ind ];
    }

    inline const T& operator()(size_t ind) const { 
      return data[ ind ];
    }

    Vec4() {};

    Vec4(T x, T y, T z, T w) 
    {
      data[0] = x;  
      data[1] = y;  
      data[2] = z;  
      data[3] = w;  
    }

  friend ostream& operator<<(ostream& os, const Vec4<T>& v)
  {
    os << "[";
    for(size_t i=0; i < 4; i++) os << v(i) << " ";
    os << "]";
    return os; 
  }
}; // end of Vec4 class


template <typename T> 
class Mat3
{
public:

  std::array<T, 9> data;

  inline T& operator()(size_t i, size_t j) { 
    return data[ j + 3*i ]; // NOTE: same ordering as in numpy
  }

  inline const T& operator()(size_t i, size_t j) const { 
    return data[ j + 3*i ]; 
  }

  Mat3() {};

  Mat3( Vec3<T>& x, Vec3<T>& y, Vec3<T>& z)
  {
    data[0] = x(0);  
    data[1] = x(1);  
    data[2] = x(2);  

    data[3] = y(0);  
    data[4] = y(1);  
    data[5] = y(2);  

    data[6] = z(0);  
    data[7] = z(1);  
    data[8] = z(2);  
  }

  friend ostream& operator<<(ostream& os, const Mat3<T>& M)
  {
    os << std::endl;
    os << "[";
    os << "[ "  << M(0, 0) << " , " << M(0, 1) << " , " << M(0, 2) << "]" << std::endl;
    os << " [ " << M(1, 0) << " , " << M(1, 1) << " , " << M(1, 2) << "]" << std::endl;
    os << " [ " << M(2, 0) << " , " << M(2, 1) << " , " << M(2, 2) << "]";
    os << "]" << std::endl;

    return os; 
  }

}; // end of 3 class


template <typename T> 
class Mat4
{
  public:

    std::array<T, 16> data;

    inline T& operator()(size_t i, size_t j) { 
      return data[ j + 4*i ]; // NOTE: same ordering as in numpy
    }

    inline const T& operator()(size_t i, size_t j) const { 
      return data[ j + 4*i ]; 
    }

    Mat4() {};

    Mat4( Vec4<T>& x, Vec4<T>& y, Vec4<T>& z, Vec4<T>& w)
    {
      data[0]  = x(0);  
      data[1]  = x(1);  
      data[2]  = x(2);  
      data[3]  = x(3);  

      data[4]  = y(0);  
      data[5]  = y(1);  
      data[6]  = y(2);  
      data[7]  = y(3);  

      data[8]  = z(0);  
      data[9]  = z(1);  
      data[10] = z(2);  
      data[11] = z(3);  

      data[12] = w(0);  
      data[13] = w(1);  
      data[14] = w(2);  
      data[15] = w(3);  
    }

  friend ostream& operator<<(ostream& os, const Mat4<T>& M)
  {
    os << std::endl;
    os << "[";
    os << "[ "  << M(0, 0) << " , " << M(0, 1) << " , " << M(0, 2) << " , " << M(0,3) << "]" << std::endl;
    os << " [ " << M(1, 0) << " , " << M(1, 1) << " , " << M(1, 2) << " , " << M(1,3) << "]" << std::endl;
    os << " [ " << M(2, 0) << " , " << M(2, 1) << " , " << M(2, 2) << " , " << M(2,3) << "]" << std::endl;
    os << " [ " << M(3, 0) << " , " << M(3, 1) << " , " << M(3, 2) << " , " << M(3,3) << "]";
    os << "]" << std::endl;

    return os; 
  }

}; // end of Mat3 class


//--------------------------------------------------
// arithmetics
// elementwise arithmetics
  
template <class T>
inline Vec3<T> operator/(Vec3<T>& v, T s) 
{
  Vec3<T> ret;
  for(size_t i=0; i<3; i++) ret(i) = v(i)/s;
  return ret;
}

template <class T>
inline Vec3<T> operator*(T s, Vec3<T>& v) 
{
  Vec3<T> ret;
  for(size_t i=0; i<3; i++) ret(i) = v(i)*s;
  return ret;
}

template <class T>
inline const Vec3<T> operator*(const T s, const Vec3<T>& v) 
{
  Vec3<T> ret;
  for(size_t i=0; i<3; i++) ret(i) = v(i)*s;
  return ret;
}


template <class T>
inline Vec3<T> operator+(Vec3<T>& v1, Vec3<T>& v2) 
{
  Vec3<T> ret;
  for(size_t i=0; i<3; i++) ret(i) = v1(i) + v2(i);
  return ret;
}


template <class T>
inline Vec3<T> operator-(Vec3<T>& v1, Vec3<T>& v2) 
{
  Vec3<T> ret;
  for(size_t i=0; i<3; i++) ret(i) = v1(i) - v2(i);
  return ret;
}


template <class T>
inline T dot(Vec3<T>& v1, Vec3<T>& v2) 
{
  double ret = 0.0;
  for(size_t i=0; i<3; i++) ret += v1(i)*v2(i);
  return ret;  
}

template <class T>
inline double sum(Vec3<T>& v)
{
  return v(0) + v(1) + v(2);
}

// matrix-vector product
template <class T>
inline Vec3<T> dot(Mat3<T>& M, Vec3<T>& v) 
{
  Vec3<T> ret;
  ret(0) = v(0)*M(0,0) + v(1)*M(0,1) + v(2)*M(0,2);
  ret(1) = v(0)*M(1,0) + v(1)*M(1,1) + v(2)*M(1,2);
  ret(2) = v(0)*M(2,0) + v(1)*M(2,1) + v(2)*M(2,2);
  return ret;
}

// matrix-vector product
template <class T>
inline Vec4<T> dot(Mat4<T>& M, Vec4<T>& v) 
{
  Vec4<T> ret;
  ret(0) = v(0)*M(0,0) + v(1)*M(0,1) + v(2)*M(0,2) + v(3)*M(0,3);
  ret(1) = v(0)*M(1,0) + v(1)*M(1,1) + v(2)*M(1,2) + v(3)*M(1,3);
  ret(2) = v(0)*M(2,0) + v(1)*M(2,1) + v(2)*M(2,2) + v(3)*M(2,3);
  ret(3) = v(0)*M(3,0) + v(1)*M(3,1) + v(2)*M(3,2) + v(3)*M(3,3);
  return ret;
}


template <class T>
inline Vec3<T> cross(Vec3<T>& v1, Vec3<T>& v2) 
{
  Vec3<T> ret;
  ret(0) =  v1(1)*v2(2) - v1(2)*v2(1);
  ret(1) = -v1(0)*v2(2) + v1(2)*v2(0);
  ret(2) =  v1(0)*v2(1) - v1(1)*v2(0);
  return ret;
}

// unit vector into direction of v1 x v2
template <class T>
inline Vec3<T> unit_cross(Vec3<T>& v1, Vec3<T>& v2) 
{
  Vec3<T> ret = cross(v1, v2);
  // TODO prevent dividing by zero
  return ret/norm(ret);
}


// matrix inversion
template <class T>
inline Mat3<T> inv(Mat3<T>& M) 
{
  Mat3<T> ret;

  double det = 0.0;
  for(int i=0; i<3; i++) det += (M(0,i) * ( M(1,(i+1)%3) * M(2,(i+2)%3) - M(1,(i+2)%3) * M(2,(i+1)%3) ));

  for(int i = 0; i<3; i++) {
  for(int j = 0; j<3; j++) {
    ret(i,j) = ( (M((j+1)%3,(i+1)%3) * M((j+2)%3,(i+2)%3) ) - (M((j+1)%3,(i+2)%3) * M( (j+2)%3,(i+1)%3)) )/ det;
  }}

  return ret;
}

//--------------------------------------------------
// vector norm
inline double norm(double x, double y, double z)
{
  return std::sqrt(x*x + y*y + z*z);
}

template <typename T> 
inline T norm(Vec3<T>& v)
{
  return std::sqrt( v(0)*v(0) + v(1)*v(1) + v(2)*v(2) );
}






} // enf of ns toolbox
