
#include "particle.h"
#include "../tools/wrap.h"
#include "../tools/iter/devcall.h"
#include "../tools/iter/iter.h"

#include <algorithm>
#include <map>
#include <utility>
#include <mpi.h>
#include <functional>

#include <cuda_runtime_api.h>
#include <nvtx3/nvToolsExt.h> 


namespace pic {

template<>
void ParticleContainer<3>::check_outgoing_particles(
    std::array<double,3>& mins,
    std::array<double,3>& maxs)
{
nvtxRangePush(__PRETTY_FUNCTION__);

  to_other_tiles.clear();
  outgoing_count = 0;

  // unpack limits
  double 
    xmin = mins[0],
    ymin = mins[1],
    zmin = mins[2],

    xmax = maxs[0],
    ymax = maxs[1],
    zmax = maxs[2];

  int lenx = static_cast<int>( xmax - xmin );
  int leny = static_cast<int>( ymax - ymin );
  int lenz = static_cast<int>( zmax - zmin );

  //DEV_REGISTER

  // shortcut for particle locations
  real_prtcl* locn[3];
  for( int i=0; i<3; i++) 
    locn[i] = &( loc(i,0) );

  #ifdef GPU_no

  //getErrorCuda((cudaMallocManaged((void**)&count, sizeof(int))));
  int maxCap = to_other_tiles.capacity();
  // 
  to_other_tiles.resize(to_other_tiles.capacity());
  auto listPtr = to_other_tiles.data();
  
  UniIter::iterate(
    [=] DEVFUNIFGPU (int n, int &count)
    {
      int i = 0;
      int j = 0;
      int k = 0;

      int i0, j0, k0;
      i0 = static_cast<int>( floor(locn[0][n] - mins[0]) );
      j0 = static_cast<int>( floor(locn[1][n] - mins[1]) );
      k0 = static_cast<int>( floor(locn[2][n] - mins[2]) );

      if(i0 <  0)    i--; // left wrap
      if(i0 >= lenx) i++; // right wrap

      if(j0 <  0)    j--; // bottom wrap
      if(j0 >= leny) j++; // top wrap

      if(k0 <  0)    k--; // back
      if(k0 >= lenz) k++; // front

      if ( (i != 0) || (j != 0) || (k != 0) ) 
      {
        //to_other_tiles.push_back( {i,j,k,n} );
        #ifdef GPU
          int pos = atomicAdd(&count, 1);
        #else
        int pos;
          #pragma omp atomic capture 
          {
            pos = count;
            count++;
          }
        #endif
        if(pos < maxCap)
          listPtr[pos] = {i,j,k,n};
      }
    },size(), outgoing_count);
    UniIter::sync();
    // check outgoing_count and react to it...
    if(outgoing_count > maxCap)
    {
      UniIter::sync();
      //std::cout << outgoing_count << " outgoing count over cap " << maxCap << std::endl;
      to_other_tiles.clear();
      to_other_tiles.reserve((outgoing_count));

      // clear, realloc and recall
      check_outgoing_particles(mins, maxs);
    }
    else{
      to_other_tiles.resize(outgoing_count);
    }

   #else
  for(size_t n=0; n<size(); n++) {
    int i,j,k; // relative indices
    i = 0;
    j = 0;
    k = 0;

    int i0, j0, k0;
    i0 = static_cast<int>( floor(locn[0][n] - mins[0]) );
    j0 = static_cast<int>( floor(locn[1][n] - mins[1]) );
    k0 = static_cast<int>( floor(locn[2][n] - mins[2]) );

    if(i0 <  0)    i--; // left wrap
    if(i0 >= lenx) i++; // right wrap

    if(j0 <  0)    j--; // bottom wrap
    if(j0 >= leny) j++; // top wrap

    if(k0 <  0)    k--; // back
    if(k0 >= lenz) k++; // front

    if ( (i != 0) || (j != 0) || (k != 0) ) 
      to_other_tiles.push_back( {i,j,k,n} );
  }
  #endif
  

nvtxRangePop();

}


}

