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
class MomentumLagrangianSolver : public VlasovVelocitySolver {

  public:
    void solve() {

      // get reference to the Vlasov fluid that we are solving
      VlasovFluid& gr          = cell->getPlasmaGrid();

      // get reference to the Yee grid 
      maxwell::YeeLattice& yee = cell->getYee();


      // loop over the cell's internal grid
      for(size_t q=0; q<gr.Nx; q++) {
        for(size_t r=0; r<gr.Ny; r++) {

          //--------------------------------------------------
          // Initialize
          vmesh::VeloMesh& vmesh = gr.electrons(q,r,0);


          // Get local field components
          double Bx = yee.bx(q,r,0);
          double By = yee.by(q,r,0);
          double Bz = yee.bz(q,r,0);

          double Ex = yee.ex(q,r,0);
          double Ey = yee.ey(q,r,0);
          double Ez = yee.ez(q,r,0);

          intp->dt = cell->dt;


          double vx, vy, vz;
          double force;

          std::array<double, 3> vel; 


          // Roll over all dimensions
          //--------------------------------------------------
          for (size_t dim=0; dim<3; dim++) {

            size_t Nb1, Nb2;
            switch(dim) {
              case 0: Nb1 = vmesh.Nblocks[1];
                      Nb2 = vmesh.Nblocks[2];
                      break;
              case 1: Nb1 = vmesh.Nblocks[0];
                      Nb2 = vmesh.Nblocks[2];
                      break;
              case 2: Nb1 = vmesh.Nblocks[0];
                      Nb2 = vmesh.Nblocks[1];
                      break;
            }

            // fmt::print("Solving for dim {} with {} {}\n", dim, Nb1, Nb2);

            bundles::Bundle delta; 
            delta.resize(vmesh.Nblocks[dim]); 
            vblock_t block;


            for(size_t i2=0; i2<Nb2; i2++) {
              for(size_t i1=0; i1<Nb1; i1++) {

                // compute cross product against B field
                // then add the proper E field component
                //
                // i.e., Lorentz Force
                switch(dim) {
                  case 0: vel = vmesh.getCenterIndx( {{0, i1, i2}} );
                          vy = vel[1];
                          vz = vel[2];
                          force = (vy*Bz - vz*By);

                          force += Ex;
                          break;
                  case 1: vel = vmesh.getCenterIndx( {{i1, 0, i2}} );
                          vx = vel[0];
                          vz = vel[2];
                          force = (vz*Bx - vx*Bz);

                          force += Ey;
                          break;
                  case 2: vel = vmesh.getCenterIndx( {{i1, i2, 0}} );
                          vx = vel[0];
                          vy = vel[1];
                          force = (vx*By - vy*Bx);

                          force += Ez;
                          break;
                }


                // create force bundle to act on the distribution
                for (size_t q=0; q<vmesh.Nblocks[dim]; q++) {
                  block[0] = force;
                  delta.loadBlock(q, block);
                }


                intp->setDelta( delta );

                // get bundle at the location
                bundles::Bundle vbundle = vmesh.getBundle(dim, i1, i2);

                // interpolate numerical flux
                intp->setBundle(vbundle);
                bundles::Bundle U0 = intp->interpolate();

                // apply flux to the mesh
                vmesh.addBundle(dim, i1, i2, U0);


              }
            }

          }

        }
      }
    }

};

} // end of namespace vlasov
