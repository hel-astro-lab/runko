#ifndef COMMON_H
#define COMMON_H

#include <limits>
#include <vector>



#define sqr(x) ((x)*(x))
#define pow2(x) sqr(x)
#define pow3(x) ((x)*(x)*(x))
#define MASTER_RANK 0


// Define uint64_t for MPI
#ifndef MPI_UINT64_T
#define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
#endif

// Define size_t for MPI
#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "Error while initializing size_t type for MPI"
#endif


// TODO remove double definition in Node (done for python bindings)
// struct globalflags {
//     static bool master;
// };


namespace conf {

    /// Grid dimensions
    const size_t Nx = 10;
    const size_t Ny = 10;

    /// block size inside spatial cell
    const size_t NxCell = 2;
    const size_t NyCell = 2;

    /// physical grid dimensions
    const double xmin = 0.0;
    const double xmax = 1.0;

    const double ymin = 0.0;
    const double ymax = 1.0;

}


namespace BC {

    /// Periodic x boundary condition
    size_t xwrap( int i ) {
        while (i < 0) {
            i += conf::Nx;
        }
        while (i >= conf::Nx) {
            i -= conf::Nx;
        }
        return size_t(i);
    }

    /// Periodic y boundary condition
    size_t ywrap( int j ) {
        while (j < 0) {
            j += conf::Ny;
        }
        while (j >= conf::Ny) {
            j -= conf::Ny;
        }
        return size_t(j);
    }

}


/*! Contains different cell boundary types. Type of computation will depend on 
 * what this type is.
 */
namespace cellType {
    enum {
        LOCAL,    //! Default type indicating that cell is owned by the current process
        VIRTUAL,  //! virtual cell (owned by other process)
        OUTFLOW,  //! outflow cell
        N_CELLTYPES
    };
}

namespace commType {
    enum {
        NCELLS,   //! Number of incoming cells,
        CELLDATA, //! Cell data array
        N_COMMTYPES
    };
}




#endif
