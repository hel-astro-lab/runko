#pragma once

#include "bundles.h"
#include "sheets.h"
#include "velomesh.h"
#include "cell.h"
#include "grid.h"


namespace vlasov {

/* -------------------------------------------------- 
 *
 * Generic interfaces for the Vlasov solvers are defined here. Actual
 * solver implementations are located in solvers/ directory.
 *
 */


/*! \brief General Vlasov velocity solver
 *
 * Defines the general interface for VELOCITY Vlasov solvers.
 * Actual solver should be implemented in solve()
 */
class VlasovVelocitySolver {

  public:
    /// Bundle interpolator pointer
    bundles::BundleInterpolator *intp;

    /// Velocity mesh to solve
    vlasov::VlasovCell* cell;

    // construct solver with a strict reference to the cell
    // VlasovVelocitySolver(vlasov::VlasovCell& cell) : cell(cell) {}

    /// Set internal mesh
    // void setMesh(vmesh::VeloMesh _vmesh) { vmesh = &_vmesh; };
    void setCell(vlasov::VlasovCell& c) { cell = &c; };

    /// Set internal interpolator
    void setInterpolator( bundles::BundleInterpolator &_intp ) { intp = &_intp; };

    /// actual solver implementation
    virtual void solve() = 0;


};




/*! \brief General Vlasov location solver
 *
 * Defines the general interface for SPATIAL Vlasov solvers.
 * Actual solver should be implemented in solve()
 *
 * Notes: Gets pointed to a focus cell
 *        Fetches neighbors using cells own interface
 *        Solves the advection equation locally
 *
 * TODO check and ensure correct initialization
 */
class VlasovSpatialSolver {

  public:

    // reference to the node
    vlasov::Grid* grid;

    /// Construct solver always with information of the node
    VlasovSpatialSolver(vlasov::Grid& grid) : grid(&grid) {}

    /// default construct 
    VlasovSpatialSolver() : grid(nullptr) {}

    /// Target cell 
    size_t targeti = -1;
    size_t targetj = -1;

    /// Set grid address so we can probe neighbors for halo data
    void setGrid(vlasov::Grid& newGrid) {
      grid = &newGrid;
    };

    /// set internal cell that is solved
    void setTargetCell(size_t i, size_t j) {
      targeti = i;
      targetj = j;
    };

    /// actual solver implementation
    virtual void solve() = 0;

};


} // end of namespace vlasov
