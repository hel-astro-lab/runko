
#include "core/pic/boundaries/gap.h"
#include "tools/vector.h"
#include "tools/signum.h"
#include "tools/staggered_grid.h"

#include <cmath> 
#include <cassert>
#include <string>
#include <map>


using std::min;
using std::max;
using std::abs;
using std::sqrt;
using std::sin;
using std::cos;

using toolbox::Vec3;
using toolbox::norm;
using toolbox::norm1d;
using toolbox::norm2d;
using toolbox::StaggeredSphericalCoordinates;
using toolbox::shape; // tanh function



template<size_t D>
void pic::Gap<D>::insert_em(
    pic::Tile<D>& tile)
{

}



template<size_t D>
void pic::Gap<D>::update_b(
    pic::Tile<D>& tile)
{

}



template<size_t D>
void pic::Gap<D>::update_e(
    pic::Tile<D>& tile)
{

}


// current from moving frame
template<size_t D>
void pic::Gap<D>::update_j(
  pic::Tile<D>& tile)
{

}


template<size_t D>
void pic::Gap<D>::solve(
    pic::Tile<D>& tile)
{


}






//--------------------------------------------------
// explicit template instantiation
template class pic::Gap<1>; // 1D only
