#include "skinny_yee.h"
#include "../emf/tile.h"


ffe::SkinnyYeeLattice& 
  ffe::SkinnyYeeLattice::operator +=(const ffe::SkinnyYeeLattice& rhs)
{
  ex += rhs.ex;
  ey += rhs.ey;
  ez += rhs.ez;
  bx += rhs.bx;
  by += rhs.by;
  bz += rhs.bz;

  return *this;
}

ffe::SkinnyYeeLattice& 
  ffe::SkinnyYeeLattice::operator -=(const ffe::SkinnyYeeLattice& rhs) 
{
  ex -= rhs.ex;
  ey -= rhs.ey;
  ez -= rhs.ez;
  bx -= rhs.bx;
  by -= rhs.by;
  bz -= rhs.bz;

  return *this;
}

ffe::SkinnyYeeLattice& 
  ffe::SkinnyYeeLattice::operator *=(double rhs) 
{
  ex *= static_cast<float_m>(rhs);
  ey *= static_cast<float_m>(rhs);
  ez *= static_cast<float_m>(rhs);
  bx *= static_cast<float_m>(rhs);
  by *= static_cast<float_m>(rhs);
  bz *= static_cast<float_m>(rhs);

  return *this;
}

ffe::SkinnyYeeLattice& 
  ffe::SkinnyYeeLattice::operator /=(double rhs) 
{
  ex /= static_cast<float_m>(rhs);
  ey /= static_cast<float_m>(rhs);
  ez /= static_cast<float_m>(rhs);
  bx /= static_cast<float_m>(rhs);
  by /= static_cast<float_m>(rhs);
  bz /= static_cast<float_m>(rhs);

  return *this;
}


/// copy yee grid to skinny yee
void ffe::SkinnyYeeLattice::set_yee(const emf::YeeLattice& yee)
{
  ex = yee.ex;
  ey = yee.ey;
  ez = yee.ez;
  bx = yee.bx;
  by = yee.by;
  bz = yee.bz;
}

ffe::SkinnyYeeLattice 
  ffe::operator +(ffe::SkinnyYeeLattice lhs, const ffe::SkinnyYeeLattice& rhs)
{
  lhs += rhs;
  return lhs;
}

ffe::SkinnyYeeLattice 
  ffe::operator -(ffe::SkinnyYeeLattice lhs, const ffe::SkinnyYeeLattice& rhs)
{
  lhs -= rhs;
  return lhs;
}

ffe::SkinnyYeeLattice 
  ffe::operator *(ffe::SkinnyYeeLattice lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}

ffe::SkinnyYeeLattice 
  ffe::operator /(ffe::SkinnyYeeLattice lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}

