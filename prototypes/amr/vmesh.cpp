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


    const uint64_t error_block = 0;
    const uint64_t error_index = 0;


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
            uint64_t nBlocks = 0;
    
            std::unordered_map<uint64_t, vmesh::vBlock> meshBlocks; // XXX soon obsolete
            std::unordered_map<uint64_t, vblock_t> blockContainer;

            std::array<double, 3> mins, maxs, dvs; // Geometry parameters

            indices_t nCells;

            void zFill( std::array<double, 3> mins_,
                        std::array<double, 3> maxs_,
                        std::array<double, 3> dvs_ );

            vblock_t get_block( uint64_t cellID );

            uint64_t get_block_ID( indices_t index );

            indices_t get_indices( uint64_t cellID );

            std::array<double, 3> get_size( uint64_t cellID );

            std::array<double, 3> get_center( uint64_t cellID );




    }; // end of vMesh class header


    void vmesh::vMesh::zFill( 
            std::array<double, 3> mins_,
            std::array<double, 3> maxs_,
            std::array<double, 3> dvs_ ) {

        double xi, yi, zi;

        // fix grid geometry
        mins   = {{ mins_[0], mins_[1], mins_[2] }};
        maxs   = {{ maxs_[0], maxs_[1], maxs_[2] }};
        dvs    = {{  dvs_[0],  dvs_[1],  dvs_[2] }};



        fmt::print("z-filling vel-space from x:{} {} | y:{} {} | z:{} {}\n",mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2]);

        // fill mesh in Morton z-order
        uint64_t indx   = 0;
        uint64_t cellID = 1;

        zi = mins[2] + dvs[2]/2.0;
        uint64_t nz = 0;
        while( zi <= maxs[2] - dvs[2]/2.0 ){
            yi = mins[1] + dvs[1]/2.0;

            uint64_t ny = 0;
            while( yi <= maxs[1] - dvs[1]/2.0 ){
                xi = mins[0] + dvs[0]/2.0;

                uint64_t nx = 0;
                while( xi <= maxs[0] - dvs[0]/2.0 ){
                    fmt::print("({},{},{})\n", xi, yi, zi);

                    vblock_t vblock = {{0.0, 0.0, 0.0, 0.0}};
                    blockContainer.insert( std::make_pair(cellID, vblock ) );

                    nBlocks++;
                    cellID++;
                    indx++;

                    xi += dvs[0];
                    nx++;
                };
                yi += dvs[1];
                ny++;

                nCells[0] = nx;
            };
            zi += dvs[2];
            nz++;

            nCells[1] = ny;
        };
        nCells[2] = nz;
        
    };


    // --------------------------------------------------

    // Get block of cell id based on the global cellID
    vblock_t vmesh::vMesh::get_block( uint64_t cellID ) {
        typename std::unordered_map<uint64_t, vblock_t>::const_iterator it = blockContainer.find(cellID);
        return it->second;
    };

    // Transform (i,j,k) indices (in z-ordering) to unique global IDs on top level of refinement
    uint64_t vmesh::vMesh::get_block_ID( indices_t index ) {

        // check for bad input
        // if (index[0] < 0)          {return vmesh::error_block;};
        // if (index[1] < 0)          {return vmesh::error_block;};
        // if (index[2] < 0)          {return vmesh::error_block;};
        if (index[0] >= nCells[0]) {return vmesh::error_block;};
        if (index[1] >= nCells[1]) {return vmesh::error_block;};
        if (index[2] >= nCells[2]) {return vmesh::error_block;};

        uint64_t GID = 1; // we start cell order from 1; 0 is error cell
        GID += index[2] * nCells[1];
        GID += index[1] * nCells[0];
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
             cellID % nCells[0], 
            (cellID /  nCells[0]) % nCells[1] ,
             cellID / (nCells[0] * nCells[1] )
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


    std::array<double, 3> vmesh::vMesh::get_size( uint64_t cellID ) {
        // TODO: check which refinement level we are on
        int refLevel = 0; 

        std::array<double, 3> len;
        for (int i=0; i<3; i++) { len[i] = dvs[i] / std::pow(2.0, refLevel); };

        return len;
    };


    std::array<double, 3> vmesh::vMesh::get_center( uint64_t cellID ) {
        // TODO check for out-of-bounds ID
        indices_t indx = get_indices( cellID );

        std::array<double, 3> center = {{0.0, 0.0, 0.0}};

        // TODO add refinement
        center = {{
				mins[0] + dvs[0]/2.0 + double(indx[0]) * dvs[0],
				mins[1] + dvs[1]/2.0 + double(indx[1]) * dvs[1],
				mins[2] + dvs[2]/2.0 + double(indx[2]) * dvs[2]
			    }};

    return center;
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
        .def_readwrite("nBlocks", &vmesh::vMesh::nBlocks)
        .def_readwrite("mins", &vmesh::vMesh::mins)
        .def_readwrite("maxs", &vmesh::vMesh::maxs)
        .def_readwrite("dvs", &vmesh::vMesh::dvs)
        .def_readwrite("nCells", &vmesh::vMesh::nCells)
        .def_readwrite("meshBlocks", &vmesh::vMesh::meshBlocks)
        .def("zFill", &vmesh::vMesh::zFill)
        .def("get_block", &vmesh::vMesh::get_block)
        .def("get_block_ID", &vmesh::vMesh::get_block_ID)
        .def("get_indices", &vmesh::vMesh::get_indices)
        .def("get_size", &vmesh::vMesh::get_size)
        .def("get_center", &vmesh::vMesh::get_center);

        // .def("get_vBlock_from_ID", &vmesh::vMesh::get_vBlock_from_ID)
        // .def("get_vBlock_from_index", &vmesh::vMesh::get_vBlock_from_index);


}
