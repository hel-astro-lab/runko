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



namespace vmesh {


    static const uint64_t error_block = 0;
    static const uint64_t error_index = 0xFFFFFFFFFFFFFFFF;

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


   // -------------------------------------------------- 
    class vMesh {

        public:
            size_t number_of_blocks = 0;
    
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
            Real threshold = 1.0e-2;


            void zFill( std::array<double, 3> mins_,
                        std::array<double, 3> maxs_);

            vblock_t get_block( const uint64_t cellID ) const;

            uint64_t get_block_ID( const indices_t index ) const;

            indices_t get_indices( uint64_t cellID );

            std::array<double, 3> get_size( const uint64_t cellID );

            std::array<double, 3> get_center( const uint64_t cellID );

            std::vector<uint64_t> all_blocks( bool sorted = false);


            std::array<Real, PENCIL_WID> get_x_pencil(size_t j, size_t k);
            std::array<Real, PENCIL_WID> get_y_pencil(size_t i, size_t k);
            std::array<Real, PENCIL_WID> get_z_pencil(size_t i, size_t j);


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

        std::array<double, 3> center = {{0.0, 0.0, 0.0}};

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




}






// --------------------------------------------------
PYBIND11_MODULE(vmesh, m) {


    py::class_<vmesh::vBlock>(m, "vBlock" )
        .def(py::init<Real, Real, Real, Real, Real, Real>())
        .def_readwrite("data",      &vmesh::vBlock::data)
        .def_readwrite("refLevel",  &vmesh::vBlock::refLevel)
        .def("loc", &vmesh::vBlock::get_loc)
        .def("dls", &vmesh::vBlock::get_dls);


    /*
            block_t get_block( uint64_t cellID );
            indices_t get_indices( uint64_t cellID );
            std::array<double, 3> get_size( uint64_t cellID );
            std::array<double, 3> get_center( uint64_t cellID );
    */


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


}
