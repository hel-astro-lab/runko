#pragma once


// Old type definitions
// FIXME: remove completely and switch to the new one below.
using Realf = double;
using Real = double;
using creal = const double;


// General double and single precision lengths used in the modules.
//
// By default, most of the arrays are stored with real_short. Some quantities
// that need more precision are temporarily casted to real_long (and often
// back to real_short in the end). I.e., most of the algorithms are mixed-precision.
//
// For extreme simulations the memory footprint can be reduced further by 
// setting real_long = float. In this case, however, you must be sure that
// physics is not affected by insufficient floating point accuracy.
//
// It is also worth mentioning that most IO files are already manually set 
// to float type to reduce storage requirements.
using real_long  = double; /// some spesific algorithmic manipulations are done via this
using real_short = float;  /// mainly used for mesh storage types

// nan prevention tolerance
#define EPS 1e-7f

// Separate short float type for particles. In most applications it is enough
// to keep this as float. However, some special cases might require extra
// resolution so fall-back to double might be needed.
//
// NOTE: Some calculations are still mixed-precision, as before, because critical
// numerical parts are still done via real_long precision.
//
// NOTE: To force a purely single-precision particle algorithms make 
// real_long also float.
using real_prtcl = float; /// particle storage type; NOTE: don't forget MPI type below

// corresponding MPI datatype for particle communications; MPI_FLOAT or MPI_DOUBLE
#define MPI_FLOAT_TP MPI_FLOAT



// FIXME: These should be cleaned away also at some point and moved to 
// relevant modules instead (like vlv).
#define NBLOCKS 20     /// default cube size
#define BLOCK_WID 4    /// block width

using indices_t = std::array<size_t, 3>;
using vblock_t = std::array<Real, BLOCK_WID>;


