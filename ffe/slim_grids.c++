#include "slim_grids.h"
#include "../emf/tile.h"


ffe::SlimGrids& 
  ffe::SlimGrids::operator +=(const ffe::SlimGrids& rhs)
{
  ex += rhs.ex;
  ey += rhs.ey;
  ez += rhs.ez;
  bx += rhs.bx;
  by += rhs.by;
  bz += rhs.bz;

  return *this;
}

ffe::SlimGrids& 
  ffe::SlimGrids::operator -=(const ffe::SlimGrids& rhs) 
{
  ex -= rhs.ex;
  ey -= rhs.ey;
  ez -= rhs.ez;
  bx -= rhs.bx;
  by -= rhs.by;
  bz -= rhs.bz;

  return *this;
}

ffe::SlimGrids& 
  ffe::SlimGrids::operator *=(double rhs) 
{
  ex *= static_cast<float_m>(rhs);
  ey *= static_cast<float_m>(rhs);
  ez *= static_cast<float_m>(rhs);
  bx *= static_cast<float_m>(rhs);
  by *= static_cast<float_m>(rhs);
  bz *= static_cast<float_m>(rhs);

  return *this;
}

ffe::SlimGrids& 
  ffe::SlimGrids::operator /=(double rhs) 
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
void ffe::SlimGrids::set_grids(const emf::Grids& yee)
{
  ex = yee.ex;
  ey = yee.ey;
  ez = yee.ez;
  bx = yee.bx;
  by = yee.by;
  bz = yee.bz;
}

ffe::SlimGrids 
  ffe::operator +(ffe::SlimGrids lhs, const ffe::SlimGrids& rhs)
{
  lhs += rhs;
  return lhs;
}

ffe::SlimGrids 
  ffe::operator -(ffe::SlimGrids lhs, const ffe::SlimGrids& rhs)
{
  lhs -= rhs;
  return lhs;
}

ffe::SlimGrids 
  ffe::operator *(ffe::SlimGrids lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}

ffe::SlimGrids 
  ffe::operator /(ffe::SlimGrids lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}

