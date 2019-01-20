#pragma once

#include <tuple>


namespace hilbert {


unsigned int generalhilbertindex(unsigned int m0, unsigned int m1,  int x,  int y);
unsigned int generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int m2,  int x,  int y,  int z);
void generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h);
void generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h);



class Hilbert2D {
  unsigned int m0, m1;

  public:

  Hilbert2D(int m0i, int m1i)
  {
    m0 = static_cast<unsigned int>(m0i);
    m1 = static_cast<unsigned int>(m1i);

    //std::cout << " creating Hilb with " << m0 << " " << m1 << " " << m0i << " " << m1i << "\n";
  }


  unsigned int hindex(int x, int y)
  {
    //std::cout << " calling hindex with " << m0 << " " << m1 << " " << x << " " << y << "\n";
    return generalhilbertindex(m0, m1, x,  y);
  }

  std::tuple<int, int> inv(unsigned int h)
  {
    unsigned int x, y;
    generalhilbertindexinv(m0, m1, &x, &y, h);

    return std::make_tuple<int, int>(x,y);
  }
};


class Hilbert3D {

  unsigned int m0, m1, m2;

  public:

  Hilbert3D(int m0i, int m1i, int m2i) 
  { 
    m0 = static_cast<unsigned int>(m0i);
    m1 = static_cast<unsigned int>(m1i);
    m2 = static_cast<unsigned int>(m2i);
  }


  unsigned int hindex(int x, int y, int z)
  {
    return generalhilbertindex(m0, m1, m2, x,  y, z);
  }

  std::tuple<int, int, int> inv(unsigned int h)
  {
    unsigned int x, y, z;
    generalhilbertindexinv(m0, m1, m2, &x, &y, &z, h);

    return std::make_tuple<int, int, int>(x,y,z);
  }
};


}
