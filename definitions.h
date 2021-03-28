#pragma once

// math constants
#define PI 3.14159265358979323846 

// General double and single precision lengths used in the modules.
//
// By default, most of the arrays are stored with float_m. Some quantities
// that need more precision are temporarily casted to double (and 
// back to float_m in the end). I.e., most of the algorithms are mixed-precision.
// It is also worth mentioning that most IO files are already manually set 
// to float type to reduce storage requirements.
using float_m = float;  /// mesh/grid floating point accuracy


// nan prevention tolerance
#define EPS 1e-7f

// Separate short float type for particles. In most applications it is enough
// to keep this as float. However, some special cases might require extra
// resolution (e.g. with large Lorentz factors) so fall-back to double might 
// be needed.
//
// NOTE: Some calculations are still mixed-precision, as before, because critical
// numerical parts are still done via double precision.
//
using float_p = float; /// particle storage type; NOTE: don't forget MPI type below

// corresponding MPI datatype for particle communications; MPI_FLOAT or MPI_DOUBLE
#define MPI_FLOAT_TP MPI_FLOAT



// FIXME: These should be cleaned away also at some point and moved to 
// relevant modules instead (like vlv).
//#define NBLOCKS 20     /// default cube size
//#define BLOCK_WID 4    /// block width

//using indices_t = std::array<size_t, 3>;
//using vblock_t = std::array<float, BLOCK_WID>;


