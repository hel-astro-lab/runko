#pragma once


typedef float Realf;
typedef float Real;
typedef const float creal;

#define NBLOCKS 20     /// default cube size
#define BLOCK_WID 4    /// block width

typedef std::array<size_t, 3> indices_t;

typedef std::array<Real, BLOCK_WID> vblock_t;


