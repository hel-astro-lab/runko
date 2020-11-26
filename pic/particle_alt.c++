
#include "particle.h"
#include "../tools/wrap.h"
#include "../tools/iter/devcall.h"
//#include "../tools/iter/iter.h"

#include <algorithm>
#include <map>
#include <utility>
#include <mpi.h>
#include <functional>

#include <cuda_runtime_api.h>
//#include <nvtx3/nvToolsExt.h> 

#include <nvtx3/nvToolsExt.h> 

namespace pic {

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
  for (auto&& elem : neigh.to_other_tiles) {
//  for (size_t ii = 0; ii < neigh.to_other_tiles.size(); ii++)
//  {
//    const auto &elem = neigh.to_other_tiles[ii];

      
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

      //std::cout << locx << " " << locy << " " << locz << " " <<  velx << " " << vely << " " << velz << " " <<  wgt << " " <<  id << " " <<  proc << std::endl;
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

