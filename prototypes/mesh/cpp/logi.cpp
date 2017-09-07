#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <fmt/format.h>
#include <fmt/format.cc>
#include <fmt/string.h>
#include <fmt/ostream.h>

#include <vector>
#include <cmath>

// #include "definitions.h"


/*
    path = 'out'
    Nx = 20
    Ny = 20
    NxCell = 2
    NyCell = 2
    xmin = ymin = 0.0
    xmax = ymax = 1.0
*/


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

};


namespace BC {

    /// Periodic x boundary condition
    size_t xwrap( int i ) {
        while (i < 0) {
            i += conf::Nx;
        }
        while (i >= conf::Nx) {
            i -= conf::Nx;
        };
        return size_t(i);
    };

    /// Periodic y boundary condition
    size_t ywrap( int j ) {
        while (j < 0) {
            j += conf::Ny;
        }
        while (j >= conf::Ny) {
            j -= conf::Ny;
        };
        return size_t(j);
    };

};

namespace logi {

    

    class Cell {

        public:

            /// unique cell ID
            uint64_t cid;

            /// MPI rank of who owns me
            size_t owner;

            /// coarse mpiGrid grid indices
            size_t i, j;


            Cell(size_t i, size_t j, size_t o) {
                this->i     = i;
                this->j     = j;
                this->owner = 0;
            };

            /// return mpiGrid index
            const std::tuple<size_t, size_t> index() {
                return std::make_tuple( i, j );
            };


            // return index of cells in relative to my position
            const std::tuple<size_t, size_t> neighs(int ir, int jr) {
                size_t ii = BC::xwrap( int(this->i) + ir );
                size_t jj = BC::ywrap( int(this->j) + jr );
                return std::make_tuple( ii, jj );
            };


            /// Return full neighborhood around me
            std::vector< std::tuple<size_t, size_t> > nhood() {
                std::vector< std::tuple<size_t, size_t> > nh;
                for (int ir=-1; ir<=1; ir++) {
                    for (int jr=-1; jr<=1; jr++) {
                        if (ir != 0 && jr != 0) {
                            nh.push_back( neighs(ir, jr) );
                        };
                    };
                };
                return nh;
            };


    };







}




// --------------------------------------------------
PYBIND11_MODULE(logi, m) {


    py::class_<logi::Cell>(m, "Cell" )
        .def(py::init<size_t, size_t, size_t >())
        .def_readwrite("cid",   &logi::Cell::cid)
        .def_readwrite("owner", &logi::Cell::owner)
        .def_readwrite("i",     &logi::Cell::i)
        .def_readwrite("j",     &logi::Cell::j)
        .def("index",  &logi::Cell::index)
        .def("neighs", &logi::Cell::neighs)
        .def("nhood",  &logi::Cell::nhood);




}
