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
    
            std::unordered_map<uint64_t, vmesh::vBlock> meshBlocks;

            std::array<Real, 3> mins, maxs, dvs; // Geometry parameters


            void zFill( std::array<Real, 3>  mins_,
                         std::array<Real, 3> maxs_,
                         std::array<Real, 3> dvs_ );

            vmesh::vBlock get_vBlock_from_ID( uint64_t cellID );



    }; // end of vMesh class header


    void vmesh::vMesh::zFill( 
            std::array<Real, 3> mins_,
            std::array<Real, 3> maxs_,
            std::array<Real, 3> dvs_ ) {

        Real xi, yi, zi;

        // fix grid geometry
        mins = {{mins_[0], mins_[1], mins_[2]}};
        maxs = {{maxs_[0], maxs_[1], maxs_[2]}};
        dvs  = {{ dvs_[0],  dvs_[1],  dvs_[2]}};

        

        zi = mins[2] + dvs[2]/2.0;
        while( zi <= maxs[2] - dvs[2]/2.0 ){
            yi = mins[1] + dvs[1]/2.0;
            while( yi <= maxs[1] - dvs[1]/2.0 ){
                xi = mins[0] + dvs[0]/2.0;
                while( xi <= maxs[0] - dvs[0]/2.0 ){
                    fmt::print("({},{},{})\n", xi, yi, zi);

                    vmesh::vBlock vb( xi, yi, zi, dvs[0], dvs[1], dvs[2] );
                    meshBlocks.insert( std::make_pair(nBlocks, vb) );
                    nBlocks++;

                    xi += dvs[0];
                };
                yi += dvs[1];
            };
            zi += dvs[2];
        };
        

    };

    vmesh::vBlock vmesh::vMesh::get_vBlock_from_ID( uint64_t cellID ) {
        typename std::unordered_map<uint64_t, vmesh::vBlock >::const_iterator it = meshBlocks.find(cellID);
        
        return it->second;
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



    py::class_<vmesh::vMesh>(m, "vMesh" )
        .def(py::init<>())
        .def_readwrite("nBlocks", &vmesh::vMesh::nBlocks)
        .def_readwrite("meshBlocks", &vmesh::vMesh::meshBlocks)
        .def("zFill", &vmesh::vMesh::zFill)
        .def("get_vBlock_from_ID", &vmesh::vMesh::get_vBlock_from_ID);


}
