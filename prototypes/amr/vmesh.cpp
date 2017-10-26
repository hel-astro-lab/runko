#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <fmt/format.h>
#include <fmt/format.cc>
#include <fmt/string.h>
#include <fmt/ostream.h>

#include <vector>
#include <cmath>
#include <unordered_map>


#include "definitions.h"



// --------------------------------------------------
namespace conf {

    /// Grid dimensions
    const size_t Nx = 50;
    const size_t Ny = 1;

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
// --------------------------------------------------

namespace vmesh {


    static const uint64_t error_block = 0;
    static const uint64_t error_index = 0xFFFFFFFFFFFFFFFF;


    //-------------------------------------------------- 
    class vBlock {
        public:
            Real data = 0.0;

            std::array<Real, 3> loc;
            std::array<Real, 3> dls;

            std::array<Real, 3> get_loc(){ return loc; };
            std::array<Real, 3> get_dls(){ return dls; };

            int refLevel = 0;

            vBlock( Real vx, Real vy, Real vz, Real dx, Real dy, Real dz);


    }; // end of vBlock class header

    // default constructor
    vBlock::vBlock( Real vx, Real vy, Real vz, Real dx, Real dy, Real dz )
    {
        loc[0] = vx;
        loc[1] = vy;
        loc[2] = vz;

        dls[0] = dx;
        dls[1] = dy;
        dls[2] = dz;
    };

    //-------------------------------------------------- 
    /// Bundle of pencils
    class Bundle {

        /// values along the pencil
        std::vector<Realf> pencil;

        /// guiding grid for the pencil
        std::vector<Realf> grid;

        public:

        /// resize the container 
        void resize( size_t N) { 
            pencil.resize(N);
            grid.resize(N);
        };

        size_t size() {
            return grid.size();
        };

        /// load bundle full of zeros
        void loadZeroBlock(size_t q) {
            pencil[q] = 0.0;
        };

        /// load block in x/y/z order 
        void loadBlock(size_t q, vblock_t block) {
            pencil[q] = block[0]; 
        };

        /// load values to the grid and transform the incoming cube according to dim
        void loadGrid(size_t q, Realf val) {
            grid[q] = val;
        };

        /// return the guiding grid
        std::vector<Realf> getGrid() {
            return grid;
        };

        /// return the pencil values
        std::vector<Realf> getPencil() {
            return pencil;
        };

        /// If current bundle slice is non-zero
        bool isNonZero(size_t q) {
            if ( pencil[q] == 0.0 ) { return false; };
            
            return true;
        };

        /// return q:th slice of the bundle
        vblock_t getSlice(size_t q) {
            vblock_t ret;
            ret[0] = pencil[q];

            return ret;
        };

        /// get grid size
        Realf getDx(size_t q) {
            return std::abs( grid[q+1] - grid[q] );
        };


    }; // end of Bundle class


    // -------------------------------------------------- 
    /// Sheet of velocities from vmesh
    // These are always slices of the full mesh along some dimension at some location
    class Sheet {

        public:

        /// guiding grid along horizontal dimensions of the sheet (i.e., x)
        std::vector<Realf> iGrid;

        // sheet size in horizontal dim
        size_t Ni = 0;

        /// guiding grid along vertical dimensions of the sheet (i.e., y)
        std::vector<Realf> jGrid;

        // sheet size in horizontal dim
        size_t Nj = 0;

        /// coordinate of the sheet in the slicing dimension
        double sliceVal;

        /// Value storage of the sheet
        std::vector<Realf> values;

        /// Resize the sheet into correct size
        void resize(size_t Ni_, size_t Nj_) {
            iGrid.resize(Ni_);
            jGrid.resize(Nj_);
            values.resize(Ni_*Nj_);

            Ni = Ni_;
            Nj = Nj_;
        }

        /// internal function to get general id from sheet indices
        size_t getIndex(size_t i, size_t j) {
            return Ni*j + i;
        }

        /// Load scalar to the sheet
        void loadValue(size_t i, size_t j, Realf val) {
            size_t indx = getIndex(i, j);
            values[indx] = val;
        }

        void loadZeroBlock(size_t i, size_t j) {
            size_t indx = getIndex(i, j);

            // FIXME add block instead of scalar
            values[indx] = 0.0;
        }

        void loadBlock(size_t i, size_t j, vblock_t block) {
            size_t indx = getIndex(i, j);

            // FIXME add block instead of scalar
            values[indx] = block[0];
        }

        // return block at location (i,j)
        vblock_t getBlock(size_t i, size_t j) {
            vblock_t ret;
            size_t indx = getIndex(i, j);
            ret[0] = values[indx]; //FIXME return block instead element

            return ret;
        }

        bool isNonZero(size_t i, size_t j) {
            size_t indx = getIndex(i, j);
            if ( values[indx] == 0.0 ) { return false; };
            return true;
        };

        // Sheet arithmetics
        Sheet& operator+=(const Sheet& rhs) {
            for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] += rhs.values[q];
            return *this;
        }

        Sheet& operator-=(const Sheet& rhs) {
            for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] -= rhs.values[q];
            return *this;
        }

        Sheet& operator*=(const Realf rhs) {
            for(size_t q=0; q<(this->Ni*this->Nj); q++) this->values[q] *= rhs;
            return *this;
        }
    };

    inline Sheet operator+(Sheet lhs, const Sheet& rhs) {
        lhs += rhs;
        return lhs;
    };

    inline Sheet operator-(Sheet lhs, const Sheet& rhs) {
        lhs -= rhs;
        return lhs;
    };

    inline Sheet operator*(Sheet lhs, const Realf rhs) {
        lhs *= rhs;
        return lhs;
    };





   // -------------------------------------------------- 
    class vMesh {

        public:
            size_t number_of_blocks = 0;
    
            /// main Hashmap block container
            std::unordered_map<uint64_t, vblock_t> blockContainer;
              
            /// returns a pointer to the data of given block
            vblock_t* operator [] (const uint64_t cellID) const
	        {
                if (this->blockContainer.count(cellID) > 0) {
	        		return (vblock_t*) &(this->blockContainer.at(cellID));
	        	} else {
	        		return NULL;
	        	}
	        }

            /* For python interface we need to define __getitem__ and __setitem__
               float operator[](size_t index) const { return m_data[index]; }
               float &operator[](size_t index) { return m_data[index]; }
            */
            vblock_t __getitem__(const uint64_t cellID) const {
                return this->blockContainer.at( cellID );
            };
            vblock_t __getitem2__(const size_t i, const size_t j, const size_t k) const {
                uint64_t cellID = this->get_block_ID( {{i,j,k}} );
                // fmt::print("({},{},{}) = {}\n",i,j,k,cellID);
                return this->__getitem__(cellID);
            };

            void __setitem__(const uint64_t cellID, const vblock_t vals) {
                blockContainer[cellID] = vals;
            }
            void __setitem2__(const size_t i, const size_t j, const size_t k, vblock_t vals) {
                uint64_t cellID = this->get_block_ID( {{i,j,k}} );
                blockContainer[cellID] = vals;
            };

            std::array<double, 3> mins, maxs, lens; // Geometry parameters
            indices_t Nblocks = {{ NBLOCKS, NBLOCKS, NBLOCKS }};

            /// Clipping threshold
            Real threshold = 1.0e-3;


            void zFill( std::array<double, 3> mins_,
                        std::array<double, 3> maxs_);

            vblock_t get_block( const uint64_t cellID ) const;

            uint64_t get_block_ID( const indices_t index ) const;

            indices_t get_indices( uint64_t cellID );

            std::array<double, 3> get_size( const uint64_t cellID );

            std::array<double, 3> get_center( const uint64_t cellID );
            std::array<double, 3> get_center_indx( const indices_t indx );

            std::vector<double> getXGrid();
            std::vector<double> getYGrid();
            std::vector<double> getZGrid();

            std::vector<uint64_t> all_blocks( bool sorted = false);

            vmesh::Bundle get_bundle(size_t, size_t, size_t);

            void add_bundle(size_t, size_t, size_t, vmesh::Bundle);

            vmesh::Sheet getSheet(size_t, size_t);
            void addSheet(size_t, size_t, vmesh::Sheet);

            bool clip( );

            size_t sizeInBytes() const;
            size_t capacityInBytes() const;



    }; // end of vMesh class header


    void vmesh::vMesh::zFill( 
            std::array<double, 3> mins_,
            std::array<double, 3> maxs_){

        // Compute grid geometry
        mins   = {{ mins_[0], mins_[1], mins_[2] }};
        maxs   = {{ maxs_[0], maxs_[1], maxs_[2] }};

        for (size_t i=0; i<3; i++) {
            lens[i] = (maxs[i] - mins[i])/( (double)Nblocks[i] - 1.0);
        }

        
       fmt::print("z-filling vel-space from x:{} {} d: {}| y:{} {} d:{}| z:{} {} d:{}\n",mins[0], maxs[0], lens[0], mins[1], maxs[1], lens[1], mins[2], maxs[2], lens[2]);

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


    // --------------------------------------------------

    // Get block of cell id based on the global cellID
    vblock_t vmesh::vMesh::get_block( const uint64_t cellID ) const {
        typename std::unordered_map<uint64_t, vblock_t>::const_iterator it = blockContainer.find(cellID);
        return it->second;
    };

    // Transform (i,j,k) indices (in z-ordering) to unique global IDs on top level of refinement
    uint64_t vmesh::vMesh::get_block_ID( const indices_t index ) const {

        // check for bad input
        // if (index[0] < 0)          {return vmesh::error_block;};
        // if (index[1] < 0)          {return vmesh::error_block;};
        // if (index[2] < 0)          {return vmesh::error_block;};
        if (index[0] >= Nblocks[0]) {return vmesh::error_block;};
        if (index[1] >= Nblocks[1]) {return vmesh::error_block;};
        if (index[2] >= Nblocks[2]) {return vmesh::error_block;};

        uint64_t GID = 1; // we start cell order from 1; 0 is error cell
        GID += index[2] * Nblocks[1]*Nblocks[0];
        GID += index[1] * Nblocks[0];
        GID += index[0];

        return GID;
    };

    indices_t vmesh::vMesh::get_indices( uint64_t cellID ) {
        if (cellID == vmesh::error_block) { 
            indices_t indx = {{ vmesh::error_index, vmesh::error_index, vmesh::error_index }};
            return indx; 
        };

        // TODO get refinement level
        // TODO substract larger cells

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
           (cell / ( this->length.get()[0] this->length.get()[1] (uint64_t(1) << (2 * refinement_level))
			)) * (uint64_t(1) << (max_refinement_level - refinement_level))
        */


        return indx;
    };


    std::array<double, 3> vmesh::vMesh::get_size( const uint64_t cellID ) {
        // TODO: check which refinement level we are on
        int refLevel = 0; 

        std::array<double, 3> wid;
        for (int i=0; i<3; i++) { wid[i] = wid[i] / std::pow(2.0, refLevel); };

        return wid;
    };


    std::array<double, 3> vmesh::vMesh::get_center( const uint64_t cellID ) {
        // TODO check for out-of-bounds ID
        indices_t indx = get_indices( cellID );

        std::array<double, 3> center;

        // TODO add refinement
        center[0] = mins[0] + lens[0]/2.0 + (double)indx[0] * lens[0];
        center[1] = mins[1] + lens[1]/2.0 + (double)indx[1] * lens[1];
        center[2] = mins[2] + lens[2]/2.0 + (double)indx[2] * lens[2];
        // fmt::print("mins {} lens {} indx {}\n", mins[0], lens[0], indx[0]);

    return center;
    };
    

    std::array<double, 3> vmesh::vMesh::get_center_indx( const indices_t indx ) {
        std::array<double, 3> center;

        // TODO add refinement
        center[0] = mins[0] + lens[0]/2.0 + (double)indx[0] * lens[0];
        center[1] = mins[1] + lens[1]/2.0 + (double)indx[1] * lens[1];
        center[2] = mins[2] + lens[2]/2.0 + (double)indx[2] * lens[2];
        // fmt::print("mins {} lens {} indx {}\n", mins[0], lens[0], indx[0]);

    return center;
    };

    /// grid along x dir
    std::vector<double> vmesh::vMesh::getXGrid() {
        std::vector<double> ret;
        ret.resize(Nblocks[0]);
        for(size_t i=0; i<Nblocks[0]; i++) {
            auto ceni = get_center_indx( {{i, 0, 0}} );
            ret[i] = ceni[0];
        }
        return ret;
    };

    /// grid along y dir
    std::vector<double> vmesh::vMesh::getYGrid() {
        std::vector<double> ret;
        ret.resize(Nblocks[1]);
        for(size_t i=0; i<Nblocks[1]; i++) {
            auto ceni = get_center_indx( {{0, i, 0}} );
            ret[i] = ceni[1];
        }
        return ret;
    };

    /// grid along z dir
    std::vector<double> vmesh::vMesh::getZGrid() {
        std::vector<double> ret;
        ret.resize(Nblocks[2]);
        for(size_t i=0; i<Nblocks[2]; i++) {
            auto ceni = get_center_indx( {{0, 0, i}} );
            ret[i] = ceni[2];
        }
        return ret;
    };
    



    /// return a list of all blocks
    std::vector<uint64_t> vmesh::vMesh::all_blocks( bool sorted ) {

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
    bool vmesh::vMesh::clip() {

        std::vector<uint64_t> below_threshold;

        for (const uint64_t block: this->all_blocks() ){
            vblock_t& blockData = blockContainer.at( block );
            // fmt::print("block: {} with data {} (len {})\n", block, blockData[0], blockData.size() );

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
    size_t vmesh::vMesh::sizeInBytes() const {
        return blockContainer.size()*sizeof(vblock_t);
    };

    // Capacity of the container because of hash map complexity and bucket division
    size_t vmesh::vMesh::capacityInBytes() const {
        return blockContainer.bucket_count() * (sizeof(vblock_t));
    };


    /// returns a sheet from the mesh that is oriented perpendicular to dim at location i
    vmesh::Sheet vmesh::vMesh::getSheet(size_t dim, size_t i) {

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
                sliceVal = get_center_indx({{i, 0, 0}})[0];
                break;
            case 1:  // y
                x = 1;
                y = 0;
                z = 2;
                horz = getXGrid();
                vert = getZGrid();
                sliceVal = get_center_indx({{0, i, 0}})[1];
                break;
            case 2:  // z
                z = 2;
                y = 0;
                z = 1;
                horz = getXGrid();
                vert = getYGrid();
                sliceVal = get_center_indx({{0, 0, i}})[2];
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
                    case 0: cid = get_block_ID( {{i, j, k}} ); // x
                            break;
                    case 1: cid = get_block_ID( {{j, i, k}} ); // y 
                            break;
                    case 2: cid = get_block_ID( {{j, k, i}} ); // z
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



    /// return full bundle of pencils penetrating the box at i1 & i2 coordinates along dim
    vmesh::Bundle vmesh::vMesh::get_bundle(size_t dim, size_t i1, size_t i2) {

        size_t Nb = Nblocks[dim];

        // target bundle
        Bundle vbundle;
        vbundle.resize( Nb );

        uint64_t cid;
        for (size_t q=0; q<Nb; q++) {

            switch(dim) {
                case 0: cid = get_block_ID( {{q, i1, i2}} ); // x pencil
                        break;
                case 1: cid = get_block_ID( {{i1, q, i2}} ); // y pencil
                        break;
                case 2: cid = get_block_ID( {{i1, i2, q}} ); // z pencil
                        break;
            }

            // add guiding grid
            auto center = get_center(cid);
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
    void vmesh::vMesh::addSheet(size_t dim, size_t i, vmesh::Sheet sheet) {

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
                    case 0: cid = get_block_ID( {{i, j, k}} ); // x
                            break;
                    case 1: cid = get_block_ID( {{j, i, k}} ); // y 
                            break;
                    case 2: cid = get_block_ID( {{j, k, i}} ); // z
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
    void vmesh::vMesh::add_bundle(size_t dim, size_t i1, size_t i2, vmesh::Bundle vbundle) {

        size_t Nb = Nblocks[dim];

        uint64_t cid=0;
        for (size_t q=1; q<Nb; q++) {

            // check if there is something coming to this block
            if (!vbundle.isNonZero(q-1) && !vbundle.isNonZero(q) ){ continue; };

            // non-zero bundle; lets add it
            switch(dim) {
                case 0: cid = get_block_ID( {{q, i1, i2}} ); // x pencil
                        break;
                case 1: cid = get_block_ID( {{i1, q, i2}} ); // y pencil
                        break;
                case 2: cid = get_block_ID( {{i1, i2, q}} ); // z pencil
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


    /// Abstract base class for bundle interpolator
    class BundleInterpolator {
        public:
            /// internal bundle that we interpolate
            Bundle bundle;

            /// force acting on the fluid
            Bundle delta;

            /// time step
            Realf dt = 0.0;

            /// numerical zero
            Realf nzero = 1.0e-4;


            virtual ~BundleInterpolator() { };

            void setBundle(Bundle _bundle) {
                bundle = _bundle;
            };

            Bundle getBundle( ) {
                return bundle;
            };

            void setDelta( Bundle _delta ) {
                delta = _delta;
            };

            vblock_t getDeltaSlice(size_t i) {
                return delta.getSlice(i);
            };

            virtual Bundle interpolate() = 0;
    };


    /// Second order Lagrangian interpolator
    class BundleInterpolator2nd : public BundleInterpolator {
        public:
            Bundle interpolate( ) {

                // prepare target bundle
                Bundle ret;
                ret.resize( bundle.size() );

                // compute flux (inner region)
                vblock_t block, fp1, f0, Delta;

                ret.loadZeroBlock(0);
                for(size_t i=1; i<bundle.size()-1; i++) {
                    fp1     = bundle.getSlice(i+1);
                    f0      = bundle.getSlice(i  );

                    // get shift 
                    Delta     = getDeltaSlice(i);
                    Delta[0] *= dt / bundle.getDx(i);

                    // 2nd order conservative Lagrangian interpolation
                    block[0] = Delta[0]          * ( fp1[0] + f0[0] )*0.5 
                             - Delta[0]*Delta[0] * ( fp1[0] - f0[0] )*0.5;

                    ret.loadBlock(i, block);
                }
                ret.loadZeroBlock( bundle.size()-1 );

                return ret;
            };
    };


    /// Fourth order Lagrangian interpolator
    class BundleInterpolator4th : public BundleInterpolator {
        public:
            Bundle interpolate( ) {

                // prepare target bundle
                Bundle ret;
                ret.resize( bundle.size() );

                // compute flux (inner region)
                vblock_t block, fp2, fp1, f0, fm1, Delta;

                ret.loadZeroBlock(0);
                ret.loadZeroBlock(1);
                for(size_t i=2; i<bundle.size()-2; i++) {
                    fm1     = bundle.getSlice(i-1);
                    f0      = bundle.getSlice(i  );
                    fp1     = bundle.getSlice(i+1);
                    fp2     = bundle.getSlice(i+2);

                    // get shift 
                    Delta     = getDeltaSlice(i);
                    Delta[0] *= dt / bundle.getDx(i);

                    // 4th order conservative Lagrangian interpolation
                    block[0] = Delta[0]        * (-fp2[0] + 7.0*fp1[0] + 7.0*f0[0] - fm1[0] )/12.0
                              +pow(Delta[0],2) * ( fp2[0] -15.0*fp1[0] +15.0*f0[0] - fm1[0] )/24.0
                              +pow(Delta[0],3) * ( fp2[0] -     fp1[0] -     f0[0] + fm1[0] )/12.0
                              +pow(Delta[0],4) * (-fp2[0] + 3.0*fp1[0] - 3.0*f0[0] + fm1[0] )/24.0;

                    ret.loadBlock(i, block);
                }
                ret.loadZeroBlock( bundle.size()-2 );
                ret.loadZeroBlock( bundle.size()-1 );

                return ret;
            };
    };




    /// General Vlasov velocity solver
    class vSolver {


        /// Bundle interpolator pointer
        BundleInterpolator *intp;

        public:
            /// Velocity mesh to solve
            vmesh::vMesh vmesh;

            // vSolver( vmesh::vMesh _vmesh ){ vmesh = _vmesh; };

            /// Set internal mesh
            void setMesh(vmesh::vMesh _vmesh) { vmesh = _vmesh; };

            /// Set internal interpolator
            void setInterpolator( BundleInterpolator &_intp ) { intp = &_intp; };


            //--------------------------------------------------
            /// actual solver
            // splitted x/y/z direction solve with static dimension rotation
            void solve( ) {

                //std::array<Realf, 3> force = {{0.2, 0.2, 0.2}};
                
                // setup force
                for (size_t dim=0; dim<3; dim++) {

                    size_t Nb1, Nb2;
                    switch(dim) {
                        case 0: Nb1 = vmesh.Nblocks[1];
                                Nb2 = vmesh.Nblocks[2];
                                break;
                        case 1: Nb1 = vmesh.Nblocks[0];
                                Nb2 = vmesh.Nblocks[2];
                                break;
                        case 2: Nb1 = vmesh.Nblocks[0];
                                Nb2 = vmesh.Nblocks[1];
                                break;
                    }

                    // fmt::print("Solving for dim {} with {} {}\n", dim, Nb1, Nb2);

                    Bundle delta; 
                    delta.resize(vmesh.Nblocks[dim]); 
                    vblock_t block;

                    /*
                    for (size_t q=0; q<vmesh.Nblocks[dim]; q++) {
                        block[0] = force[dim];
                        delta.loadBlock(q, block);
                    }
                    intp->setDelta( delta );
                    */

                    double Bx = 0.0;
                    double By = 0.0;
                    double Bz = 0.001;

                    intp->dt = 1.0;

                    // loop over other directions
                    double vx, vy, vz;
                    double force;

                    std::array<double, 3> vel; 

                    for(size_t i2=0; i2<Nb2; i2++) {
                        for(size_t i1=0; i1<Nb1; i1++) {

                            // compute cross product against B field
                            switch(dim) {
                                case 0: vel = vmesh.get_center_indx( {{0, i1, i2}} );
                                        vy = vel[1];
                                        vz = vel[2];
                                        force = (vy*Bz - vz*By);
                                        break;
                                case 1: vel = vmesh.get_center_indx( {{i1, 0, i2}} );
                                        vx = vel[0];
                                        vz = vel[2];
                                        force = (vz*Bx - vx*Bz);
                                        break;
                                case 2: vel = vmesh.get_center_indx( {{i1, i2, 0}} );
                                        vx = vel[0];
                                        vy = vel[1];
                                        force = (vx*By - vy*Bx);
                                        break;
                            }

                            // create force bundle to act on the distribution
                            for (size_t q=0; q<vmesh.Nblocks[dim]; q++) {
                                block[0] = force;
                                delta.loadBlock(q, block);
                            }


                            intp->setDelta( delta );

                            // get bundle at the location
                            Bundle vbundle = vmesh.get_bundle(dim, i1, i2);

                            // interpolate numerical flux
                            intp->setBundle(vbundle);
                            Bundle U0 = intp->interpolate();

                            // apply flux to the mesh
                            vmesh.add_bundle(dim, i1, i2, U0);
                        }
                    }

                }

            };



    }; // end of vSolver


    /// Container class for dealing with *actual* simulation data
    class DataContainer {
        std::vector<vmesh::vMesh> container;

        public:
              
            size_t currentStep = 0;

            /// method to add data into the container
            void push_back(vmesh::vMesh vm) {
                container.push_back(vm);
            }

            /// Get current element
            vmesh::vMesh* get() {
                // fmt::print("getting from DataContainer with {}\n", currentStep);
                return (vmesh::vMesh*) &(container[ currentStep ]);
            }

            vmesh::vMesh* getNew() {
                if (currentStep == 0) return (vmesh::vMesh*) &(container[1]);
                if (currentStep == 1) return (vmesh::vMesh*) &(container[0]);
            }

            vmesh::vMesh* getAll(size_t cs) {
                // fmt::print("pulling from DataContainer with {}\n", cs);
                return (vmesh::vMesh*) &(container[cs]);
            }
                

            // FIXME raw cycling for time step index
            void cycle() {
                // fmt::print(" calling cycle (originally {})\n", currentStep);
                currentStep++;

                // check bounds and cycle back
                if (currentStep > 1) currentStep = 0;
            }

    };

     

    /// XXX this is just a simplified hollow copy from LoGi
    class Cell {
        public:

            uint64_t cid;
            size_t i, j;
            int owner;

            /// Data container
            vmesh::DataContainer data;


            //-------------------------------------------------- 
            // initalize cell according to its location (i,j) and owner (o)
            Cell(size_t i, size_t j, int o) {
                this->i     = i;
                this->j     = j;
                this->owner = o;
            }

            const std::tuple<size_t, size_t> index() {
                return std::make_tuple( i, j );
            }

            const std::tuple<size_t, size_t> neighs(int ir, int jr) {
                size_t ii = BC::xwrap( (int)this->i + ir );
                size_t jj = BC::ywrap( (int)this->j + jr );
                return std::make_tuple( ii, jj );
            }


            //-------------------------------------------------- 
            // XXX new sugar on top of the logi interface
              
            void addData(vmesh::vMesh m) {
                data.push_back(m);
            }

            /*
            std::vector<vmesh::vMesh> getData() {
                std::vector<vmesh::vMesh> ret;
                ret.push_back( data.get() );
                return ret;
            }
            */

            // XXX defined only for python API
            vmesh::vMesh getData() {
                return *data.get();
                // return *data.getAll(cstep);
            };
              



    };

    // Dummy class for temporarily handling Cells
    // XXX this is just a hollow copy of the Node class in LoGi
    class Node {
        public:

            //-------------------------------------------------- 
            //-------------------------------------------------- 
            // old inherited stuff from LoGi
            std::unordered_map<uint64_t, vmesh::Cell> cells;

            /// Create unique cell ids based on Morton z ordering
            uint64_t cell_id(size_t i, size_t j) {
                return uint64_t( j*conf::Nx + i );
            }
              
            void add_local_cell( vmesh::Cell c ) {
                c.cid   = cell_id(c.i, c.j);
                c.owner = 0;
                
                cells.insert( std::make_pair(c.cid, c) );
            }

            vmesh::Cell* get_cell_data(const uint64_t cid) const {
            	if (this->cells.count(cid) > 0) {
            		return (vmesh::Cell*) &(this->cells.at(cid));
            	} else {
            		return NULL;
            	}
            }
              

            vmesh::Cell get_cell( uint64_t cid ) {
                return cells.at(cid);
            }

            vmesh::Cell get_cell_index(size_t i, size_t j) {
                uint64_t cid = cell_id(i,j);
                return cells.at(cid);
            }

            //-------------------------------------------------- 
            //-------------------------------------------------- 
            // XXX new sugar on top of the logi interface

            void cycle() {
                // for (auto it: cells) it.second.data.cycle();

                std::unordered_map<uint64_t, vmesh::Cell>::iterator it = cells.begin();
                while (it != cells.end()) {
                    it->second.data.cycle();
                    it++;
                }

            }
            


    };



    /// General Vlasov spatial solver
    class sSolver {
        // NOTES:
        // gets pointed to a focus cell;
        // fetches neighbors using cells interface
        // solves locally

        /// Bundle interpolator pointer
        // BundleInterpolator *intp;


        public:
            /// Spatial cell to solve
            // vmesh::Cell cell;

            // reference to the node
            vmesh::Node& node;
            
            /// Construct solver always with information of the node
            sSolver(vmesh::Node& node) : node(node) {}

            /// Target cell 
            size_t targeti, targetj;

            /// Set node address so we can probe neighbors for halo data
            /*
            void setNode(vmesh::Node n) {
                node = n;
            };
            */

            /// set internal cell that is solved
            void setTargetCell(size_t i, size_t j) {
                targeti = i;
                targetj = j;
            };

            /// Set internal interpolator
            // void setInterpolator( BundleInterpolator &_intp ) { intp = &_intp; };


            //--------------------------------------------------
            /// actual solver
            // splitted x/y/z direction solve with static dimension rotation
            // loop over y and z and solve for each x
            //
            // get interpolator minimum pencil width
            // create minimum pencil for each x value
            // put into bundle 
            // apply interpolator to these bundles
            // unload from return bundle and apply to *current* cell 
            void solve( ) {

                // size_t Nintp = 4;

                // get target cell that we operate on
                uint64_t cid = node.cell_id(targeti, targetj);
                Cell* cellPtr = node.get_cell_data(cid);
                  
                // Get pointer to the velocity mesh we are solving
                vMesh* v0 = cellPtr->data.get();
                // v0->clip();

                // get pointer to a new mesh 
                vMesh* v0new = cellPtr->data.getNew();


                Realf dt = 1.0e-2;
                Realf dx = 1.0;


                // for (size_t dim=0; dim<3; dim++) {
                // XXX only x dir update is done here
                for (size_t dim=0; dim<1; dim++) {
                    // fmt::print("Solving for dim {}\n", dim);

                    size_t Nb = v0->Nblocks[dim]; // number of slices to loop over
                  
                    // get neighbors for interpolation
                    auto nindx_p1    = cellPtr->neighs(+1, 0); // i+1 neighbor
                    uint64_t cid_p1  = node.cell_id( std::get<0>(nindx_p1), std::get<1>(nindx_p1) );
                    Cell* cellPtr_p1 = node.get_cell_data( cid_p1 );
                    vMesh* vp1       = cellPtr_p1->data.get();


                    auto nindx_m1    = cellPtr->neighs(-1, 0); // i-1 neighbor
                    uint64_t cid_m1  = node.cell_id( std::get<0>(nindx_m1), std::get<1>(nindx_m1) );
                    Cell* cellPtr_m1 = node.get_cell_data( cid_m1 );
                    vMesh* vm1       = cellPtr_m1->data.get();

                    // loop over every sheet in the mesh
                    Sheet Up, Um, flux;
                    Sheet s0, sp1, sm1;
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
            }

            /*
            void update() {
                for (auto it: node.cells) it->second.data.cycle();
            }
            */


    }; // end of sSolver





}






// --------------------------------------------------
PYBIND11_MODULE(vmesh, m) {


    py::class_<vmesh::vBlock>(m, "vBlock" )
        .def(py::init<Real, Real, Real, Real, Real, Real>())
        .def_readwrite("data",      &vmesh::vBlock::data)
        .def_readwrite("refLevel",  &vmesh::vBlock::refLevel)
        .def("loc", &vmesh::vBlock::get_loc)
        .def("dls", &vmesh::vBlock::get_dls);


    py::class_<vmesh::Bundle>(m, "Bundle" )
        .def(py::init<>())
        .def("getGrid",   &vmesh::Bundle::getGrid)
        .def("getPencil", &vmesh::Bundle::getPencil);


    py::class_<vmesh::Sheet>(m, "Sheet" )
        .def(py::init<>())
        .def_readwrite("iGrid",      &vmesh::Sheet::iGrid)
        .def_readwrite("jGrid",      &vmesh::Sheet::jGrid)
        .def_readwrite("Ni",         &vmesh::Sheet::Ni)
        .def_readwrite("Nj",         &vmesh::Sheet::Nj)
        .def_readwrite("values",     &vmesh::Sheet::values);



    // Trampoline class for BundleInterpolator (because of inheritance from Base class)
    class PyBundleInterpolator : public vmesh::BundleInterpolator {
    public:
        using vmesh::BundleInterpolator::BundleInterpolator;
        using vmesh::BundleInterpolator::setBundle;
        using vmesh::BundleInterpolator::getBundle;
        vmesh::Bundle interpolate() override {
            PYBIND11_OVERLOAD_PURE(vmesh::Bundle, vmesh::BundleInterpolator, interpolate, );
        }
    };

    py::class_<vmesh::BundleInterpolator, PyBundleInterpolator> bintp(m, "BundleInterpolator" );
    bintp
        .def(py::init<>())
        .def("setBundle",   &vmesh::BundleInterpolator::setBundle)
        .def("getBundle",   &vmesh::BundleInterpolator::getBundle)
        .def("interpolate", &vmesh::BundleInterpolator::interpolate);

    py::class_<vmesh::BundleInterpolator2nd>(m, "BundleInterpolator2nd", bintp)
        .def(py::init<>());

    py::class_<vmesh::BundleInterpolator4th>(m, "BundleInterpolator4th", bintp)
        .def(py::init<>());




    py::class_<vmesh::vMesh>(m, "vMesh" )
        .def(py::init<>())
        .def_readwrite("number_of_blocks", &vmesh::vMesh::number_of_blocks)
        .def_readwrite("mins", &vmesh::vMesh::mins)
        .def_readwrite("maxs", &vmesh::vMesh::maxs)
        .def_readwrite("lens", &vmesh::vMesh::lens)
        .def_readwrite("Nblocks", &vmesh::vMesh::Nblocks)
        .def("zFill", &vmesh::vMesh::zFill)
        .def("get_block", &vmesh::vMesh::get_block)
        .def("get_block_ID", &vmesh::vMesh::get_block_ID)
        .def("get_indices", &vmesh::vMesh::get_indices)
        .def("all_blocks", &vmesh::vMesh::all_blocks)
        .def("get_size", &vmesh::vMesh::get_size)
        .def("get_center", &vmesh::vMesh::get_center)
        .def("get_bundle", &vmesh::vMesh::get_bundle)
        .def("add_bundle", &vmesh::vMesh::add_bundle)
        .def("getSheet",   &vmesh::vMesh::getSheet)
        // Bare bones array interface
        /*
        .def("__getitem__", [](const Sequence &s, size_t i) {
            if (i >= s.size()) throw py::index_error();
            return s[i];
        })
        .def("__setitem__", [](Sequence &s, size_t i, float v) {
            if (i >= s.size()) throw py::index_error();
            s[i] = v;
        })
        // Slices [optional]
        .def("__getitem__", [](const Sequence &s, py::slice slice) -> Sequence* {
            size_t start, stop, step, slicelength;
            if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();
            Sequence *seq = new Sequence(slicelength);
            for (size_t i = 0; i < slicelength; ++i) {
                (*seq)[i] = s[start]; start += step;
            }
            return seq;
        })
        .def("__setitem__", [](Sequence &s, py::slice slice, const Sequence &value) {
            size_t start, stop, step, slicelength;
            if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
                throw py::error_already_set();
            if (slicelength != value.size())
                throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
            for (size_t i = 0; i < slicelength; ++i) {
                s[start] = value[i]; start += step;
            }
        })        */
        .def("__getitem__", [](const vmesh::vMesh &s, uint64_t i) {
                return s.__getitem__(i);
        })
        .def("__setitem__", [](vmesh::vMesh &s, uint64_t i, vblock_t v) {
                s.__setitem__(i, v);
        })
        // i,j,k indexing based interface
        .def("__getitem__", [](const vmesh::vMesh &s, py::tuple indx) {
                size_t i = indx[0].cast<size_t>();
                size_t j = indx[1].cast<size_t>();
                size_t k = indx[2].cast<size_t>();
                return s.__getitem2__( i,j,k );
        })
        .def("__setitem__", [](vmesh::vMesh &s, py::tuple indx, vblock_t v) {
                size_t i = indx[0].cast<size_t>();
                size_t j = indx[1].cast<size_t>();
                size_t k = indx[2].cast<size_t>();
                return s.__setitem2__( i,j,k, v);
        })

        // .def("__setitem__", [](vmesh::vMesh &s, uint64_t i, vblock_t v) {
        //         s.__setitem__(i, v);
        // })

        // other more advanced mesh manipulation functions
        .def("sizeInBytes", &vmesh::vMesh::sizeInBytes)
        .def("capacityInBytes", &vmesh::vMesh::capacityInBytes)
        .def("clip", &vmesh::vMesh::clip);


    py::class_<vmesh::vSolver>(m, "vSolver" )
        .def(py::init<>())
        .def_readwrite("vmesh",  &vmesh::vSolver::vmesh)
        .def("setMesh" ,         &vmesh::vSolver::setMesh)
        .def("setInterpolator",  &vmesh::vSolver::setInterpolator)
        .def("solve",            &vmesh::vSolver::solve);

    py::class_<vmesh::sSolver>(m, "sSolver" )
        .def(py::init<vmesh::Node&>())
        // .def_readwrite("node",   &vmesh::sSolver::node)
        // .def("setNode" ,         &vmesh::sSolver::setNode)
        .def("setTargetCell" ,   &vmesh::sSolver::setTargetCell)
        .def("solve",            &vmesh::sSolver::solve);
        // .def("update",           &vmesh::sSolver::update);


    py::class_<vmesh::Cell>(m, "Cell" )
        .def(py::init<size_t, size_t, int >())
        .def_readwrite("owner",                       &vmesh::Cell::owner)
        .def_readwrite("i",                           &vmesh::Cell::i)
        .def_readwrite("j",                           &vmesh::Cell::j)
        .def("addData",                               &vmesh::Cell::addData)
        .def("getData",                               &vmesh::Cell::getData)
        .def("index",                                 &vmesh::Cell::index)
        .def("neighs",                                &vmesh::Cell::neighs);

    py::class_<vmesh::Node>(m, "Node" )
        .def(py::init<>())
        .def("cell_id",           &vmesh::Node::cell_id)
        .def("get_cell",          &vmesh::Node::get_cell)
        .def("get_cell_index",    &vmesh::Node::get_cell_index)
        .def("add_local_cell",    &vmesh::Node::add_local_cell)
        .def("cycle",             &vmesh::Node::cycle);


}
