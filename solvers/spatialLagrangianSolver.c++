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
        // TODO fix solver interface to take cell AND node as inputs; this will avoid pointer casting?
        uint64_t cid = grid->cellId(targeti, targetj);
        Grid::CellPtr cellPtr = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid));


        // Get reference to the velocity mesh we are solving
        VlasovFluid& gr0 = cellPtr->getPlasmaGrid();

        // get reference to a new mesh that we are updating into
        VlasovFluid& gr1 = cellPtr->getNewPlasmaGrid();

        // initialize cell and step size
        Realf dt = cellPtr->dt;
        Realf dx = cellPtr->dx;


        // TODO only x dir update is done here
        for (size_t dim=0; dim<1; dim++) {
          // fmt::print("Solving for dim {}\n", dim);


          // get neighbors for interpolation
          auto nindx_p1    = cellPtr->neighs(+1, 0); // i+1 neighbor
          uint64_t cid_p1  = grid->cellId( std::get<0>(nindx_p1), std::get<1>(nindx_p1) );
          Grid::CellPtr cellPtr_p1 = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_p1));

          VlasovFluid& gr_p1 = cellPtr_p1->getPlasmaGrid();


          auto nindx_m1    = cellPtr->neighs(-1, 0); // i+1 neighbor
          uint64_t cid_m1  = grid->cellId( std::get<0>(nindx_m1), std::get<1>(nindx_m1) );
          Grid::CellPtr cellPtr_m1 = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_m1));

          VlasovFluid&     gr_m1 = cellPtr_m1->getPlasmaGrid();

          // catch non-implemented usage
          if (gr0.Nx != 1 && gr0.Ny != 1 && gr0.Nz != 1) throw std::range_error("not implemetned");

          // TODO make this loop over all elements instead
          //-------------------------------------------------- 
          vmesh::VeloMesh& v0 = gr0.electrons(0,0,0);
          vmesh::VeloMesh& v1 = gr1.electrons(0,0,0);

          vmesh::VeloMesh& vp1 = gr_p1.electrons(0, 0, 0); // xxx | 0, 1, 2, 
          vmesh::VeloMesh& vm1 = gr_m1.electrons(0, 0, 0); // 0, 1, 2, .. | xxx

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
          //-------------------------------------------------- 
            
            
        } // end of dimension cycle


        return;
      }

    };




} // end of namespace vlasov

