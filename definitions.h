#pragma once


// Old type definitions
// FIXME: remove completely and switch to the new one below.
using Realf = double;
using Real = double;
using creal = const double;


// General double and single precision lengths used in the modules.
//
// By default, most of the arrays are stored with float_t. Some quantities
// that need more precision are temporarily casted to double_t (and often
// back to float_t in the end). I.e., most of the algorithms are mixed-precision.
//
// For extreme simulations the memory footprint can be reduced further by 
// setting double_t = float. In this case, however, you must be sure that
// physics is not affected by insufficient floating point accuracy.
//
// It is also worth mentioning that most IO files are already manually set 
// to float type to reduce storage requirements.
using double_t = double;
using float_t = float;

// Separate short float type for particles. In most applications it is enough
// to keep this as float. However, some special cases might require extra
// resolution so fall-back to double might be needed.
//
// NOTE: Some calculations are still mixed-precision, as before, because critical
// numerical parts are still done via double_t precision.
//
// NOTE: To force a purely single-precision particle algorithms make 
// double_t also float.
using float_tp = float;

// corresponding MPI datatype for particle communications; MPI_FLOAT or MPI_DOUBLE
#define MPI_FLOAT_TP MPI_FLOAT



// FIXME: These should be cleaned away also at some point and moved to 
// relevant modules instead (like vlv).
#define NBLOCKS 20     /// default cube size
#define BLOCK_WID 4    /// block width

typedef std::array<size_t, 3> indices_t;

typedef std::array<Real, BLOCK_WID> vblock_t;


