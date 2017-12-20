#include <vector>
#include <cmath>
#include <stdexcept>

#include "../solvers.h"
#include "../sheets.h"


namespace vlasov {

  /* \brief Splitted Lagrangian spatial solver for Vlasov fluids
   *
   * Full solution is obtained by splitting each spatial dimension
   * x/y/z into separate solutions that are then cycled over to
   * obtain full 3D solution.
   *
   * Uses Sheets to vectorize calculations.
   *
   *
   */
  class SpatialLagrangianSolver2nd : public VlasovSpatialSolver {

    public:
      void solve() {

        // size_t Nintp = 4;

        // get target cell that we operate on
        // TODO fix solver interface to take cell AND node as inputs; 
        // this will avoid pointer casting?
        uint64_t cid = grid->cellId(targeti, targetj);
        Grid::CellPtr cellPtr = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid));


        // Get reference to the velocity mesh we are solving
        VlasovFluid& gr0 = cellPtr->getPlasmaGrid();

        // get reference to a new mesh that we are updating into
        VlasovFluid& gr1 = cellPtr->getNewPlasmaGrid();


        // numerical flux from moving the Vlasov fluid around
        Realf fluxX, fluxY, fluxZ;
        fluxX = fluxY = fluxZ = 0.0;


        // TODO only x dir update is done here; multidimensionalize
        for (size_t dim=0; dim<1; dim++) {
          // fmt::print("Solving for dim {}\n", dim);


          // get neighbors for interpolation
          auto nindx_p1    = cellPtr->neighs(+1, 0); // i+1 neighbor
          uint64_t cid_p1  = grid->cellId( std::get<0>(nindx_p1), std::get<1>(nindx_p1) );
          Grid::CellPtr cellPtr_p1 = 
            std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_p1));

          VlasovFluid& gr_p1 = cellPtr_p1->getPlasmaGrid();


          auto nindx_m1    = cellPtr->neighs(-1, 0); // i-1 neighbor
          uint64_t cid_m1  = grid->cellId( std::get<0>(nindx_m1), std::get<1>(nindx_m1) );
          Grid::CellPtr cellPtr_m1 = 
            std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_m1));

          VlasovFluid&     gr_m1 = cellPtr_m1->getPlasmaGrid();

          if (gr0.Nx == 1 && gr0.Ny == 1 && gr0.Nz == 1) {

            // No-block formatting at all
            throw std::range_error("not implemented");

            vmesh::VeloMesh& v0 = gr0.electrons(0,0,0);
            vmesh::VeloMesh& v1 = gr1.electrons(0,0,0);

            vmesh::VeloMesh& vp1 = gr_p1.electrons(0, 0, 0); // xxx | 0, 1, 2, 
            vmesh::VeloMesh& vm1 = gr_m1.electrons(0, 0, 0); // 0, 1, 2, .. | xxx

            fluxX += solve1d(v0, vm1, vp1, v1, dim, cellPtr);

          } else {

            size_t first = 0;
            size_t last  = gr0.Nx-1;

            // Block structured cells
            for(size_t k = 0; k<gr0.Nz; k++) {
              for(size_t j = 0; j<gr0.Ny; j++) {

                // leftmost side blocks (-1 value from left neighbor)
                fluxX += solve1d(
                    gr0.  electrons(first,   j,k),
                    gr_m1.electrons(last,    j,k),
                    gr0.  electrons(first+1, j,k),
                    gr1.  electrons(first,   j,k),
                    dim, cellPtr);

                // inner blocks
                for(size_t i=1; i<gr0.Nx-1; i++) {
                  
                  fluxX += solve1d(
                      gr0.electrons(i,   j,k),
                      gr0.electrons(i-1, j,k),
                      gr0.electrons(i+1, j,k),
                      gr1.electrons(i,   j,k),
                      dim, cellPtr);

                }

                // rightmost side blocks (+1 value from right neighbor)
                fluxX += solve1d(
                    gr0.  electrons(last,   j,k),
                    gr0.  electrons(last-1, j, k),
                    gr_p1.electrons(first,  j, k),
                    gr1.  electrons(last,   j,k),
                    dim, cellPtr);

              }
            }
          
          }



          //-------------------------------------------------- 
          // TODO collect and deposit flux from the sheets to the Yee lattice
            



        } // end of dimension cycle

        return;
      }


      /*! \brief Solve 1D advection problem 
       *
       * Internally we use 2D Lagrangian interpolation.
       */

      Realf solve1d(vmesh::VeloMesh& v0,
                   vmesh::VeloMesh& vm1,
                   vmesh::VeloMesh& vp1,
                   vmesh::VeloMesh& v1,
                   size_t dim,
                   Grid::CellPtr cellPtr) {

        Realf numerical_flux = 0.0;

        // initialize cell and step size
        Realf dt = cellPtr->dt;
        Realf dx = cellPtr->dx;


        size_t Nb = v0.Nblocks[dim]; // number of slices to loop over

        // loop over every sheet in the mesh
        sheets::Sheet Up, Um, flux;
        sheets::Sheet s0, sp1, sm1;
        Realf aa;
        for(size_t i=0; i<Nb; i++) {
          sm1 = vm1.getSheet(dim, i);
          s0  = v0 .getSheet(dim, i);
          sp1 = vp1.getSheet(dim, i);

          aa = (Realf) s0.sliceVal * (dt/dx);

          // U_+1/2
          Up = (sp1 + s0)*0.5* aa   
            -(sp1 - s0)*0.5* aa*aa;

          // U_-1/2
          Um = (s0 + sm1)*0.5* aa
            -(s0 - sm1)*0.5* aa*aa;

          // dF = U_-1/2 - U_+1/2
          flux = s0 + (Um - Up);

          v1.addSheet(dim, i,  flux);
        }

        return numerical_flux;
      }



    };




} // end of namespace vlasov

