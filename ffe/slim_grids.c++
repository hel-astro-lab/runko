#include "slim_grids.h"
#include "emf/tile.h"


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
  ex *= static_cast<float>(rhs);
  ey *= static_cast<float>(rhs);
  ez *= static_cast<float>(rhs);
  bx *= static_cast<float>(rhs);
  by *= static_cast<float>(rhs);
  bz *= static_cast<float>(rhs);

  return *this;
}

ffe::SlimGrids& 
  ffe::SlimGrids::operator /=(double rhs) 
{
  ex /= static_cast<float>(rhs);
  ey /= static_cast<float>(rhs);
  ez /= static_cast<float>(rhs);
  bx /= static_cast<float>(rhs);
  by /= static_cast<float>(rhs);
  bz /= static_cast<float>(rhs);

  return *this;
}


/// copy gs grid to skinny gs
void ffe::SlimGrids::set_grids(const emf::Grids& gs)
{
  ex = gs.ex;
  ey = gs.ey;
  ez = gs.ez;
  bx = gs.bx;
  by = gs.by;
  bz = gs.bz;
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

