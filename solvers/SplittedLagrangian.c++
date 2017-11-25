#include <vector>
#include <cmath>

#include "../solvers.h"
#include "../bundles.h"


namespace vlasov {

/* \brief Splitted Lagrangian Velocity solver for Vlasov fluids
 *
 * Full solution is obtained by splitting each velocity dimension
 * vx/vy/vz into separate solutions that are then cycled over to
 * obtain full 3D solution.
 *
 * Uses BundeInterpolator that can be set to whatever order 
 * interpolator one wishes to use.
 *
 *
 * TODO: dimension rotation is static (i.e., always the same order)
 *       but it should be trivial to update it into Strang splitting 
 *       to get 2nd order accuracy.
 *
 */
class SplittedLagrangian : public VlasovVelocitySolver {

  public:
    void solve() {

      // operate on a pointer to the velomesh
      vmesh::VeloMesh* vmesh = cell->getDataPtr();

      // setup force
      for (size_t dim=0; dim<3; dim++) {

        size_t Nb1, Nb2;
        switch(dim) {
          case 0: Nb1 = vmesh->Nblocks[1];
                  Nb2 = vmesh->Nblocks[2];
                  break;
          case 1: Nb1 = vmesh->Nblocks[0];
                  Nb2 = vmesh->Nblocks[2];
                  break;
          case 2: Nb1 = vmesh->Nblocks[0];
                  Nb2 = vmesh->Nblocks[1];
                  break;
        }

        // fmt::print("Solving for dim {} with {} {}\n", dim, Nb1, Nb2);

        bundles::Bundle delta; 
        delta.resize(vmesh->Nblocks[dim]); 
        vblock_t block;

        /*
           for (size_t q=0; q<vmesh.Nblocks[dim]; q++) {
           block[0] = force[dim];
           delta.loadBlock(q, block);
           }
           intp->setDelta( delta );
           */

        double Bx = 0.0;
        double By = 0.0;
        double Bz = 0.001;

        intp->dt = 1.0;

        // loop over other directions
        double vx, vy, vz;
        double force;

        std::array<double, 3> vel; 

        for(size_t i2=0; i2<Nb2; i2++) {
          for(size_t i1=0; i1<Nb1; i1++) {

            // compute cross product against B field
            switch(dim) {
              case 0: vel = vmesh->getCenterIndx( {{0, i1, i2}} );
                      vy = vel[1];
                      vz = vel[2];
                      force = (vy*Bz - vz*By);
                      break;
              case 1: vel = vmesh->getCenterIndx( {{i1, 0, i2}} );
                      vx = vel[0];
                      vz = vel[2];
                      force = (vz*Bx - vx*Bz);
                      break;
              case 2: vel = vmesh->getCenterIndx( {{i1, i2, 0}} );
                      vx = vel[0];
                      vy = vel[1];
                      force = (vx*By - vy*Bx);
                      break;
            }

            // add some extra force from E field
            force += 0.01;

            // create force bundle to act on the distribution
            for (size_t q=0; q<vmesh->Nblocks[dim]; q++) {
              block[0] = force;
              delta.loadBlock(q, block);
            }


            intp->setDelta( delta );

            // get bundle at the location
            bundles::Bundle vbundle = vmesh->getBundle(dim, i1, i2);

            // interpolate numerical flux
            intp->setBundle(vbundle);
            bundles::Bundle U0 = intp->interpolate();

            // apply flux to the mesh
            vmesh->addBundle(dim, i1, i2, U0);
          }
        }




      }

    };

};

} // end of namespace vlasov
