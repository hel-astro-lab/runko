
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

  DEV_REGISTER

  // shortcut for particle locations
  real_prtcl* locn[3];
  for( int i=0; i<3; i++) 
    locn[i] = &( loc(i,0) );


  //getErrorCuda((cudaMallocManaged((void**)&count, sizeof(int))));
  int maxCap = to_other_tiles.capacity();
  // 
  to_other_tiles.resize(to_other_tiles.capacity());
  auto listPtr = to_other_tiles.data();
  UniIter::sync();
  
  //#ifdef GPU
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
/*    
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
  *
  */

nvtxRangePop();

}


template<>
void ParticleContainer<3>::transfer_and_wrap_particles( 
    ParticleContainer& neigh,
    std::array<int,3>    dirs, 
    std::array<double,3>& global_mins, 
    std::array<double,3>& global_maxs
    )
{
nvtxRangePush(__PRETTY_FUNCTION__);

  // particle overflow from tiles is done in shortest precision
  // to avoid rounding off errors and particles left in a limbo
  // between tiles.
  real_prtcl locx, locy, locz, velx, vely, velz, wgt;
  int id, proc;

  int i;
  //for (auto&& elem : neigh.to_other_tiles) {
  for (size_t ii = 0; ii < neigh.to_other_tiles.size(); ii++)
  {
    const auto &elem = neigh.to_other_tiles[ii];

      
    if(elem.i == 0 && 
       elem.j == 0 &&
       elem.k == 0) continue; 

    // NOTE: directions are flipped (- sign) so that they are
    // in directions in respect to the current tile

    if (elem.i == -dirs[0] &&
        elem.j == -dirs[1] &&
        elem.k == -dirs[2] ) {

      i = elem.n;

      locx = wrap( neigh.loc(0, i), static_cast<real_prtcl>(global_mins[0]), static_cast<real_prtcl>(global_maxs[0]) );
      locy = wrap( neigh.loc(1, i), static_cast<real_prtcl>(global_mins[1]), static_cast<real_prtcl>(global_maxs[1]) );
      locz = wrap( neigh.loc(2, i), static_cast<real_prtcl>(global_mins[2]), static_cast<real_prtcl>(global_maxs[2]) );

      velx = neigh.vel(0, i);
      vely = neigh.vel(1, i);
      velz = neigh.vel(2, i);

      wgt  = neigh.wgt(i);

      id   = neigh.id(0,i);
      proc = neigh.id(1,i);

      //add_identified_particle({locx,locy,locz}, {velx,vely,velz}, wgt, id, proc);
      
        //assert(prtcl_loc.size() == 3);
        //assert(prtcl_vel.size() == 3);

        //for (size_t i=0; i<3; i++) 
        locArr[0].push_back(locx);
        locArr[1].push_back(locy);
        locArr[2].push_back(locz);
        //for (size_t i=0; i<3; i++) 
        velArr[0].push_back(velx);
        velArr[1].push_back(vely);
        velArr[2].push_back(velz);
        wgtArr.push_back(wgt);

        indArr[0].push_back(id);
        indArr[1].push_back(proc);

        Nprtcls++;
      
    }
  }
nvtxRangePop();

  }

}


