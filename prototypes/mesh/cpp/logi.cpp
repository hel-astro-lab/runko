#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <fmt/format.h>
#include <fmt/format.cc>
#include <fmt/string.h>
#include <fmt/ostream.h>

#include <vector>
#include <algorithm>
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

            /// Cell type listing
            std::vector<int> types = { cellType::LOCAL };


            /// initalize cell according to its location (i,j) and owner (o)
            Cell(size_t i, size_t j, size_t o) {
                this->i     = i;
                this->j     = j;
                this->owner = 0;
            }

            /// return mpiGrid index
            const std::tuple<size_t, size_t> index() {
                return std::make_tuple( i, j );
            }


            /// return index of cells in relative to my position
            const std::tuple<size_t, size_t> neighs(int ir, int jr) {
                size_t ii = BC::xwrap( int(this->i) + ir );
                size_t jj = BC::ywrap( int(this->j) + jr );
                return std::make_tuple( ii, jj );
            }


            /// Return full neighborhood around me
            std::vector< std::tuple<size_t, size_t> > nhood() {
                std::vector< std::tuple<size_t, size_t> > nh;
                for (int ir=-1; ir<=1; ir++) {
                    for (int jr=-1; jr<=1; jr++) {
                        if (ir != 0 && jr != 0) {
                            nh.push_back( neighs(ir, jr) );
                        }
                    }
                }
                return nh;
            }


            /// Check if cell fulfills a single criteria
            bool is_type( int criteria ) {
                if( std::find(
                            types.begin(), 
                            types.end(), 
                            criteria) 
                        == types.end() 
                  ) {
                    return false;
                } 
                return true;
            }
                
            /// Vectorized version requiring cell to fulfill every criteria
            bool is_types( std::vector<int> criteria ) {
                for (auto crit: criteria) {
                    if (is_type(crit))  {
                            continue;
                    } else {
                        return false;
                    }
                }

                // passed all criteria
                return true;
            }

    }; // end of Cell class



    class Node {


        /// Global large scale grid where information
        // of all the mpi processes is stored
        int _mpiGrid[conf::Nx][conf::Ny];

        /// Map with cellID & cell data
        std::unordered_map<uint64_t, logi::Cell> cells;


        public:

            /// get mpi process for whatever location
            const int mpiGrid(const int i, const int j) {
                return _mpiGrid[i][j];
            }

            /// set new mpi process for some cell
            void set_mpiGrid(const int i, const int j, int val) {
                _mpiGrid[i][j] = val;
            }


            /// Create unique cell ids based on Morton z ordering
            uint64_t cell_id(size_t i, size_t j) {
                return uint64_t( j*conf::Nx + i );
            }
            

            /// Add local cell to the node
            void add_cell( logi::Cell c ) {

                // calculate unique global cell ID
                uint64_t cid = cell_id(c.i, c.j);

                //TODO Catch error if cell is not already mine?
                c.owner = rank;

                cells.insert( std::make_pair(cid, c) );
            }

            logi::Cell get_cell( uint64_t cid ) {
                return cells.at(cid);
            }


            /*! Return a vector of cell indices that fulfill a given criteria.
             *  By default all local cells are returned.
             */
            std::vector<uint64_t> get_cells(
                    const std::vector<int>& criteria = std::vector<int>(),
                    const bool sorted=false ) {
                std::vector<uint64_t> ret;

                for (auto it: cells) {
                    if (criteria.size() == 0) {
                        ret.push_back( it.first );
                        continue;
                    }

                    // criteria checking
                    auto c = it.second;
                    if (!c.is_types( criteria ) ) {
                        continue;
                    }
                      
                    ret.push_back( it.first );
                }


                // optional sort based on the cell id
                if (sorted && ret.size() > 0) {
			        std::sort(ret.begin(), ret.end());
		        }

                return ret;
            }

            /// Return all cells that are of VIRTUAL type.
            std::vector<uint64_t> get_virtuals(
                    const std::vector<int>& criteria = std::vector<int>(),
                    const bool sorted=false ) {

                std::vector<int> new_criteria = criteria;
                new_criteria.push_back( cellType::VIRTUAL );
                return get_cells(new_criteria, sorted);
            }
            


        public:
            // -------------------------------------------------- 
            /// MPI communication related stuff
            int rank  = 0;
            int Nrank = 0;
            MPI_Comm comm;

            //TODO double definition for python debugging
            bool master = false;


            Node() {
                fmt::print("initializing node...");

            }


            /// Initialize MPI and related auxiliary variables
            void init_mpi() {

                // TODO do this in main program with arg and argv
                MPI_Init(NULL, NULL);

                comm = MPI_COMM_WORLD;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &Nrank);

                // detect master
                if (rank == MASTER_RANK) { master = true; };

                fmt::print("Hi from rank {}\n", rank);
                if (master) { fmt::print("master is {}\n", rank); };
            }


            /// Finalize MPI environment 
            void finalize_mpi() {
                MPI_Finalize();
            }


            /// Broadcast master ranks mpiGrid to everybody
            void bcast_mpiGrid() {

                MPI_Bcast(&(_mpiGrid[0][0]), 
                        conf::Nx*conf::Ny, 
                        MPI_INT, 
                        MASTER_RANK, 
                        MPI_COMM_WORLD);

            }


            // -------------------------------------------------- 
            





    }; // end of Node class




} // end of logi namespace




// --------------------------------------------------
PYBIND11_MODULE(logi, m) {

    m.attr("Nx")     = conf::Nx;
    m.attr("Ny")     = conf::Ny;
    m.attr("NxCell") = conf::NxCell;
    m.attr("NyCell") = conf::NyCell;
    m.attr("xmin")   = conf::xmin;
    m.attr("xmax")   = conf::xmax;
    m.attr("ymin")   = conf::ymin;
    m.attr("ymax")   = conf::ymax;

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
        .def_readwrite("rank",    &logi::Node::rank)
        .def_readwrite("Nrank",   &logi::Node::Nrank)
        .def_readwrite("master",  &logi::Node::master)
        .def("mpiGrid",           &logi::Node::mpiGrid)
        .def("add_cell",          &logi::Node::add_cell)
        .def("get_cell",          &logi::Node::get_cell)
        .def("get_cells",         &logi::Node::get_cells,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)
        .def("get_virtuals",      &logi::Node::get_virtuals,
                py::arg("criteria") = std::vector<int>(),
                py::arg("sorted") = true)

        // communication wrappers
        .def("set_mpiGrid",       &logi::Node::set_mpiGrid)
        .def("init_mpi",          &logi::Node::init_mpi)
        .def("bcast_mpiGrid",     &logi::Node::bcast_mpiGrid)
        .def("finalize_mpi",      &logi::Node::finalize_mpi);



}
