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

            std::vector<uint64_t> all_blocks( bool sorted = false);

            vmesh::Bundle get_bundle(size_t, size_t, size_t);

            void add_bundle(size_t, size_t, size_t, vmesh::Bundle);


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


    /// XXX this is just a simplified hollow copy from LoGi
    class Cell {
        public:

            uint64_t cid;
            size_t i, j;
            int owner;

            /// Data container
            std::vector<vmesh::vMesh> dataContainer;


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

            void addData(vmesh::vMesh m) {
                dataContainer.push_back(m);
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


    };












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


    py::class_<vmesh::Cell>(m, "Cell" )
        .def(py::init<size_t, size_t, int >())
        .def_readwrite("owner",                       &vmesh::Cell::owner)
        .def_readwrite("i",                           &vmesh::Cell::i)
        .def_readwrite("j",                           &vmesh::Cell::j)
        .def_readwrite("dataContainer",               &vmesh::Cell::dataContainer)
        .def("addData",                               &vmesh::Cell::addData)
        .def("index",                                 &vmesh::Cell::index)
        .def("neighs",                                &vmesh::Cell::neighs);

    py::class_<vmesh::Node>(m, "Node" )
        .def(py::init<>())
        .def("cell_id",           &vmesh::Node::cell_id)
        .def("get_cell",          &vmesh::Node::get_cell)
        .def("get_cell_index",    &vmesh::Node::get_cell_index)
        .def("add_local_cell",    &vmesh::Node::add_local_cell);


}
