
#include <array>
#include <vector>

#include "bundles.h"


using namespace bundles;



/// resize the container 
void Bundle::resize( size_t N) { 
    pencil.resize(N);
    grid.resize(N);
};

size_t Bundle::size() {
    return grid.size();
};

/// load bundle full of zeros
void Bundle::loadZeroBlock(size_t q) {
    pencil[q] = 0.0;
};

/// load block in x/y/z order 
void Bundle::loadBlock(size_t q, vblock_t block) {
    pencil[q] = block[0]; 
};

/// load values to the grid and transform the incoming cube according to dim
void Bundle::loadGrid(size_t q, Realf val) {
    grid[q] = val;
};

/// return the guiding grid
std::vector<Realf> Bundle::getGrid() {
    return grid;
};

/// return the pencil values
std::vector<Realf> Bundle::getPencil() {
    return pencil;
};

/// If current bundle slice is non-zero
bool Bundle::isNonZero(size_t q) {
    if ( pencil[q] == 0.0 ) { return false; };
    
    return true;
};

/// return q:th slice of the bundle
vblock_t Bundle::getSlice(size_t q) {
    vblock_t ret;
    ret[0] = pencil[q];

    return ret;
};

/// get grid size
Realf Bundle::getDx(size_t q) {
    return std::abs( grid[q+1] - grid[q] );
};
