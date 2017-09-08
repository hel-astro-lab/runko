#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <fmt/format.h>
#include <fmt/format.cc>
#include <fmt/string.h>
#include <fmt/ostream.h>

#include <vector>
#include <cmath>

#include "mpi.h"

// #include "definitions.h"
#include "common.h"










namespace logi {


    class Cell {

        public:

            /// unique cell ID
            uint64_t cid;

            /// MPI rank of who owns me
            size_t owner;

            /// coarse mpiGrid grid indices
            size_t i, j;


            /// initalize cell according to its location (i,j) and owner (o)
            Cell(size_t i, size_t j, size_t o) {
                this->i     = i;
                this->j     = j;
                this->owner = 0;
            };

            /// return mpiGrid index
            const std::tuple<size_t, size_t> index() {
                return std::make_tuple( i, j );
            };


            /// return index of cells in relative to my position
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

    }; // end of Cell class



    class Node {

        public:

            /// MPI stuff
            int rank  = 0;
            int Nrank = 0;
            MPI_Comm comm;

            //TODO double definition for python debugging
            bool master = false;


            Node() {
                fmt::print("initializing node...");

            };


            /// Initialize MPI and related auxiliary variables
            void initMPI() {

                // TODO do this in main program arg and argv
                MPI_Init(NULL, NULL);

                comm = MPI_COMM_WORLD;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &Nrank);

                // detect master
                if (rank == MASTER_RANK) { master = true; };

                fmt::print("Hi from rank {}\n", rank);
                if (master) { fmt::print("master is {}\n", rank); };
            };


            /// Finalize MPI environment 
            void finalizeMPI() {
                MPI_Finalize();
            };



    }; // end of Node class




} // end of logi namespace




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

    py::class_<logi::Node>(m, "Node" )
        .def(py::init<>())
        .def_readwrite("rank",   &logi::Node::rank)
        .def_readwrite("Nrank",  &logi::Node::Nrank)
        .def("initMPI",          &logi::Node::initMPI)
        .def("finalizeMPI",          &logi::Node::finalizeMPI);



}
