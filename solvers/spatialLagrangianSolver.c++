#include <vector>
#include <cmath>

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
        uint64_t cid = grid->cellId(targeti, targetj);
        Grid::CellPtr cellPtr = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid));


        // Get pointer to the velocity mesh we are solving
        vmesh::VeloMesh* v0 = cellPtr->data.get();
        // v0->clip();

        // get pointer to a new mesh 
        vmesh::VeloMesh* v0new = cellPtr->data.getNew();


        Realf dt = 1.0e-2;
        Realf dx = 1.0;


        // for (size_t dim=0; dim<3; dim++) {
        // XXX only x dir update is done here
        for (size_t dim=0; dim<1; dim++) {
          // fmt::print("Solving for dim {}\n", dim);

          size_t Nb = v0->Nblocks[dim]; // number of slices to loop over

          // get neighbors for interpolation
          auto nindx_p1    = cellPtr->neighs(+1, 0); // i+1 neighbor
          uint64_t cid_p1  = grid->cellId( std::get<0>(nindx_p1), std::get<1>(nindx_p1) );
          Grid::CellPtr cellPtr_p1 = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_p1));
          vmesh::VeloMesh* vp1       = cellPtr_p1->data.get();

          auto nindx_m1    = cellPtr->neighs(-1, 0); // i+1 neighbor
          uint64_t cid_m1  = grid->cellId( std::get<0>(nindx_m1), std::get<1>(nindx_m1) );
          Grid::CellPtr cellPtr_m1 = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_m1));
          vmesh::VeloMesh* vm1       = cellPtr_m1->data.get();


          // loop over every sheet in the mesh
          sheets::Sheet Up, Um, flux;
          sheets::Sheet s0, sp1, sm1;
          Realf aa;
          for(size_t i=0; i<Nb; i++) {
            sm1 = vm1->getSheet(dim, i);
            s0  = v0 ->getSheet(dim, i);
            sp1 = vp1->getSheet(dim, i);

            aa = (Realf) s0.sliceVal * (dt/dx);

            // U_+1/2
            Up = (sp1 + s0)*0.5* aa   
              -(sp1 - s0)*0.5* aa*aa;

            // U_-1/2
            Um = (s0 + sm1)*0.5* aa
              -(s0 - sm1)*0.5* aa*aa;

            // dF = U_-1/2 - U_+1/2
            flux = s0 + (Um - Up);

            v0new ->addSheet(dim, i,  flux);
          }
        } // end of dimension cycle

        return;
      }

    };




} // end of namespace vlasov

