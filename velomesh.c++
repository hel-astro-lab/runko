#include <vector>
#include <cmath>
#include <unordered_map>


#include "bundles.h"
#include "velomesh.h"
#include "definitions.h"

using namespace vmesh;
using namespace sheets;
using namespace bundles;


/// Morton Z-fill velocity mesh with blocks
void VeloMesh::zFill( 
        std::array<double, 3> mins_,
        std::array<double, 3> maxs_){

    // Compute grid geometry
    mins   = {{ mins_[0], mins_[1], mins_[2] }};
    maxs   = {{ maxs_[0], maxs_[1], maxs_[2] }};

    for (size_t i=0; i<3; i++) {
        lens[i] = (maxs[i] - mins[i])/( (double)Nblocks[i] - 1.0);
    }

    
   // fmt::print("z-filling from x:{} {} d: {}| y:{} {} d:{}| z:{} {} d:{}\n",
   //     mins[0], maxs[0], lens[0], 
   //     mins[1], maxs[1], lens[1], 
   //     mins[2], maxs[2], lens[2]);


    // fill mesh in Morton z-order
    uint64_t cellID = 1;
    for(size_t zi=0; zi<Nblocks[2]; zi++) {
        for(size_t yi=0; yi<Nblocks[1]; yi++) {
            for(size_t xi=0; xi<Nblocks[0]; xi++) {
                // fmt::print("({},{},{}\n",xi,yi,zi);
                vblock_t vblock;
                blockContainer.insert( std::make_pair(cellID, vblock ) );

                number_of_blocks++;
                cellID++;
            }
        }
    }
    

}



/// Get block of cell id based on the global cellID
vblock_t VeloMesh::getBlock( const uint64_t cellID ) const {
    typename std::unordered_map<uint64_t, vblock_t>::const_iterator it 
      = blockContainer.find(cellID);

    return it->second;
};


/*! Transform (i,j,k) indices (in z-ordering) to unique global IDs on top level 
 *  of refinement.
 */
uint64_t VeloMesh::getBlockID( const indices_t index ) const {

    // check for bad input
    //
    // these are ignored because we use size_t type that cant be negative
    // if (index[0] < 0)          {return error_block;};
    // if (index[1] < 0)          {return error_block;};
    // if (index[2] < 0)          {return error_block;};
    if (index[0] >= Nblocks[0]) {return error_block;};
    if (index[1] >= Nblocks[1]) {return error_block;};
    if (index[2] >= Nblocks[2]) {return error_block;};

    uint64_t GID = 1; // we start cell order from 1; 0 is error cell
    GID += index[2] * Nblocks[1]*Nblocks[0];
    GID += index[1] * Nblocks[0];
    GID += index[0];

    return GID;
};

/// Get indices from cellid number
indices_t VeloMesh::getIndices( uint64_t cellID ) {
    if (cellID == error_block) { 
        indices_t indx = {{ error_index, error_index, error_index }};
        return indx; 
    };

    // TODO get refinement level
    // TODO subtract larger cells

    cellID -= 1; // numbering starts from zero

    indices_t indx = {{ 
         cellID %  Nblocks[0], 
        (cellID /  Nblocks[0]) % Nblocks[1] ,
         cellID / (Nblocks[0]  * Nblocks[1] )
         }};


    /* (cell % (this->length.get()[0] * (uint64_t(1) << refinement_level)))
		* (uint64_t(1) << (max_refinement_level - refinement_level))

       ((cell / (this->length.get()[0] * (uint64_t(1) << refinement_level)))
			% (this->length.get()[1] * (uint64_t(1) << refinement_level)))
		* (uint64_t(1) << (max_refinement_level - refinement_level)),
       (cell / ( this->length.get()[0] this->length.get()[1] (uint64_t(1) << 
       (2 * refinement_level))
		)) * (uint64_t(1) << (max_refinement_level - refinement_level))
    */


    return indx;
};

/// get cell size from cell id
std::array<double, 3> VeloMesh::getSize( const uint64_t cellID ) {
    // TODO: check which refinement level we are on
    int refLevel = 0; 

    std::array<double, 3> wid;
    for (int i=0; i<3; i++) { wid[i] = wid[i] / std::pow(2.0, refLevel); };

    return wid;
};


/// get cell center from cell id
std::array<double, 3> VeloMesh::getCenter( const uint64_t cellID ) {
    // TODO check for out-of-bounds ID
    indices_t indx = getIndices( cellID );

    std::array<double, 3> center;

    // TODO add refinement
    center[0] = mins[0] + lens[0]/2.0 + (double)indx[0] * lens[0];
    center[1] = mins[1] + lens[1]/2.0 + (double)indx[1] * lens[1];
    center[2] = mins[2] + lens[2]/2.0 + (double)indx[2] * lens[2];
    // fmt::print("mins {} lens {} indx {}\n", mins[0], lens[0], indx[0]);

return center;
};


// get cell center from cell index
// TODO use operator overloading instead (in pybind section)
std::array<double, 3> VeloMesh::getCenterIndx( const indices_t indx ) {
    std::array<double, 3> center;

    // TODO add refinement
    center[0] = mins[0] + lens[0]/2.0 + (double)indx[0] * lens[0];
    center[1] = mins[1] + lens[1]/2.0 + (double)indx[1] * lens[1];
    center[2] = mins[2] + lens[2]/2.0 + (double)indx[2] * lens[2];
    // fmt::print("mins {} lens {} indx {}\n", mins[0], lens[0], indx[0]);

return center;
};


/// grid along x dir
std::vector<double> VeloMesh::getXGrid() {
    std::vector<double> ret;
    ret.resize(Nblocks[0]);
    for(size_t i=0; i<Nblocks[0]; i++) {
        auto ceni = getCenterIndx( {{i, 0, 0}} );
        ret[i] = ceni[0];
    }
    return ret;
};


/// grid along y dir
std::vector<double> VeloMesh::getYGrid() {
    std::vector<double> ret;
    ret.resize(Nblocks[1]);
    for(size_t i=0; i<Nblocks[1]; i++) {
        auto ceni = getCenterIndx( {{0, i, 0}} );
        ret[i] = ceni[1];
    }
    return ret;
};


/// grid along z dir
std::vector<double> VeloMesh::getZGrid() {
    std::vector<double> ret;
    ret.resize(Nblocks[2]);
    for(size_t i=0; i<Nblocks[2]; i++) {
        auto ceni = getCenterIndx( {{0, 0, i}} );
        ret[i] = ceni[2];
    }
    return ret;
};


/// return a list of all blocks
std::vector<uint64_t> VeloMesh::allBlocks( bool sorted ) {

    std::vector<uint64_t> ret_val;

    for (auto item: blockContainer) {
        uint64_t cellID = item.first;

        ret_val.push_back( cellID );
    };
        

    if (sorted && ret_val.size() > 0) {
		std::sort(ret_val.begin(), ret_val.end());
	}

    return ret_val;
};



/// Clip all the blocks below threshold
bool VeloMesh::clip() {

    std::vector<uint64_t> below_threshold;

    for (const uint64_t block: this->allBlocks() ){
        vblock_t& blockData = blockContainer.at( block );
        // fmt::print("block: {} with data {} (len {})\n", 
        // block, blockData[0], blockData.size() );

        if (blockData[0] < threshold) { below_threshold.push_back( block ); }
    };

    // TODO collect / bookkeep how much was lost
    for (uint64_t block: below_threshold) { 
        // fmt::print("Erasing {}\n", block);
        blockContainer.erase( block ); 
        number_of_blocks -= 1;
    };


    return true;
};


/// Bytes actually used by the container
size_t VeloMesh::sizeInBytes() const {
    return blockContainer.size()*sizeof(vblock_t);
};


// Capacity of the container because of hash map complexity and bucket division
size_t VeloMesh::capacityInBytes() const {
    return blockContainer.bucket_count() * (sizeof(vblock_t));
};


/*! returns a sheet from the mesh that is oriented perpendicular to dim at 
 * location i
 */
Sheet VeloMesh::getSheet(size_t dim, size_t i) {

    // get i,j,k elements rotated along the correct dimension
    size_t x,y,z;
    std::vector<double> horz, vert;
    double sliceVal;
    switch(dim) {
        case 0:  // x
            x = 0;
            y = 1;
            z = 2;
            horz = getYGrid();
            vert = getZGrid();
            sliceVal = getCenterIndx({{i, 0, 0}})[0];
            break;
        case 1:  // y
            x = 1;
            y = 0;
            z = 2;
            horz = getXGrid();
            vert = getZGrid();
            sliceVal = getCenterIndx({{0, i, 0}})[1];
            break;
        case 2:  // z
            z = 2;
            y = 0;
            z = 1;
            horz = getXGrid();
            vert = getYGrid();
            sliceVal = getCenterIndx({{0, 0, i}})[2];
            break;
    }

    // create & initialize sheet
    Sheet sheet;
    sheet.resize(Nblocks[y], Nblocks[z]);
    for(size_t q=0; q<horz.size(); q++) sheet.iGrid[q] = horz[q];
    for(size_t q=0; q<vert.size(); q++) sheet.jGrid[q] = vert[q]; 
    sheet.sliceVal = sliceVal;


    // get values
    uint64_t cid;
    for(size_t k=0; k<Nblocks[z]; k++) {
        for(size_t j=0; j<Nblocks[y]; j++) {
            switch(dim) {
                case 0: cid = getBlockID( {{i, j, k}} ); // x
                        break;
                case 1: cid = getBlockID( {{j, i, k}} ); // y 
                        break;
                case 2: cid = getBlockID( {{j, k, i}} ); // z
                        break;
            }

            // lets get actual values 
            auto it = blockContainer.find(cid);

            // if block does not exist, fill with zero
            if( it == blockContainer.end() ) {
                sheet.loadZeroBlock(j,k);
                continue;
            }

            // TODO transform block to correct order
            vblock_t block = it->second;
            sheet.loadBlock(j, k, block);
        }
    }


    return sheet;
};



/*! return full bundle of pencils penetrating the box at i1 & i2 coordinates 
 * along dim
 */
Bundle VeloMesh::getBundle(size_t dim, size_t i1, size_t i2) {

    size_t Nb = Nblocks[dim];

    // target bundle
    Bundle vbundle;
    vbundle.resize( Nb );

    uint64_t cid;
    for (size_t q=0; q<Nb; q++) {

        switch(dim) {
            case 0: cid = getBlockID( {{q, i1, i2}} ); // x pencil
                    break;
            case 1: cid = getBlockID( {{i1, q, i2}} ); // y pencil
                    break;
            case 2: cid = getBlockID( {{i1, i2, q}} ); // z pencil
                    break;
        }

        // add guiding grid
        auto center = getCenter(cid);
        vbundle.loadGrid(q, center[dim] );

        // next lets get actual values 
        auto it = blockContainer.find(cid);

        // if block does not exist, fill with zero
        if( it == blockContainer.end() ) {
            vbundle.loadZeroBlock(q);
            continue;
        }

        // TODO transform block to correct order
        vblock_t block = it->second;

        vbundle.loadBlock(q, block);
    }

    return vbundle;
};

/// Add given sheet to the mesh
void VeloMesh::addSheet(size_t dim, size_t i, Sheet sheet) {

    // rotate dimensions to match incoming sheet
    size_t x,y,z;
    switch(dim) {
        case 0:  // x
            x = 0;
            y = 1;
            z = 2;
            break;
        case 1:  // y
            x = 1;
            y = 0;
            z = 2;
            break;
        case 2:  // z
            z = 2;
            y = 0;
            z = 1;
            break;
    }


    uint64_t cid;
    for(size_t k=0; k<Nblocks[z]; k++) {
        for(size_t j=0; j<Nblocks[y]; j++) {

            // check if there is something coming to this block
            if(!sheet.isNonZero(j,k) ) continue; 

            // non-zero block; lets add it
            switch(dim) {
                case 0: cid = getBlockID( {{i, j, k}} ); // x
                        break;
                case 1: cid = getBlockID( {{j, i, k}} ); // y 
                        break;
                case 2: cid = getBlockID( {{j, k, i}} ); // z
                        break;
            }

            vblock_t vb = sheet.getBlock(j, k);

            // next lets get correct block
            auto it = blockContainer.find(cid);

            // if block does not exist, create it 
            if( it == blockContainer.end() ) {
                blockContainer.insert( std::make_pair(cid, vb ) );
                continue;
            }

            // non-zero address block
            // TODO: real, full block, addition
            vblock_t targetBlock = it->second;
            targetBlock[0] = vb[0];

            it->second = targetBlock;
        } // end of loop over sheet dimensions
    }

};


/// Add given bundle to the right blocks along the corresponding dimension
void VeloMesh::addBundle(size_t dim, size_t i1, size_t i2, Bundle vbundle) {

    size_t Nb = Nblocks[dim];

    uint64_t cid=0;
    for (size_t q=1; q<Nb; q++) {

        // check if there is something coming to this block
        if (!vbundle.isNonZero(q-1) && !vbundle.isNonZero(q) ){ continue; };

        // non-zero bundle; lets add it
        switch(dim) {
            case 0: cid = getBlockID( {{q, i1, i2}} ); // x pencil
                    break;
            case 1: cid = getBlockID( {{i1, q, i2}} ); // y pencil
                    break;
            case 2: cid = getBlockID( {{i1, i2, q}} ); // z pencil
                    break;
        }

        // get slice and compute flux *difference*
        // TODO rotate the block to correct dimensions
        vblock_t vb0  = vbundle.getSlice(q);
        vblock_t vbm1 = vbundle.getSlice(q-1);
        vblock_t vb;
        vb[0] = vbm1[0] - vb0[0];



        // next lets get correct block
        auto it = blockContainer.find(cid);
        
        // if block does not exist, create it 
        if( it == blockContainer.end() ) {
            blockContainer.insert( std::make_pair(cid, vb ) );
            continue;
        }

        // non-zero address block
        // TODO: real, full block, addition
        vblock_t targetBlock = it->second;
          
        targetBlock[0] += vb[0];

        it->second = targetBlock;
    }

};



