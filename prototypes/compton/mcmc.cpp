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



namespace mcmc {

    class photon {

        public:
            std::array<double, 4> data;
            photon( double, double, double, double );

            /// photon four-velocity components
            const double hv(){ return data[0]; };
            const double vx(){ return data[1]; };
            const double vy(){ return data[2]; };
            const double vz(){ return data[3]; };

    };


    photon::photon(double E, double vx, double vy, double vz) {
        this->data = {{E, vx, vy, vz}};
    };


    class electron {

        public:
            std::array<double, 4> data; 
            electron( double, double, double, double );
           
            // four-velocity of electron
            const double v0(){ return data[0]; };
            const double vx(){ return data[1]; };
            const double vy(){ return data[2]; };
            const double vz(){ return data[3]; };

            /// 3-velocity in units of c, i.e. v/c
            const double v() {
                return std::sqrt(
                          std::pow(vx(), 2.0) 
                        + std::pow(vy(), 2.0) 
                        + std::pow(vz(), 2.0) 
                                ); };

            /// gamma factor
            const double gamma() { return 1.0 / std::sqrt( 1.0 - std::pow( v(), 2.0 ) ); };

            /// proper momentum
            const double p() { return gamma() * v(); };

    };

    electron::electron(double v0, double vx, double vy, double vz) {
        this->data = {{v0, vx, vy, vz}};
    };



    // --------------------------------------------------
    class photonBucket {

        /// size of the bucket (in terms of photons)
        size_t nPhotons = 0;

        /// Photon container
        std::vector<std::array<double, 4>> bucket;

        public:
            void push_back( photon ph );
            const size_t size( ) { return this->nPhotons; };

            void replace( size_t indx, photon ph );

            photon get( size_t indx );

    };

    void photonBucket::push_back( photon ph ) {
        bucket.push_back( ph.data );
        nPhotons++;
    };

    void photonBucket::replace( size_t indx, photon ph ) {
        bucket[indx] = ph.data;
    };

    photon photonBucket::get( size_t indx ) {
        std::array<double, 4> data = bucket[indx];
        photon ph(data[0], data[1], data[2], data[3]);

        return ph;
    };



}






// --------------------------------------------------
PYBIND11_MODULE(mcmc, m) {


    py::class_<mcmc::photon>(m, "photon" )
        .def(py::init<double, double, double, double >())
        .def_readwrite("data",      &mcmc::photon::data)
        .def("hv",  &mcmc::photon::hv)
        .def("vx",  &mcmc::photon::vx)
        .def("vy",  &mcmc::photon::vy)
        .def("vz",  &mcmc::photon::vz);

    py::class_<mcmc::electron>(m, "electron" )
        .def(py::init<double, double, double, double >())
        .def_readwrite("data",      &mcmc::electron::data)
        .def("v0",  &mcmc::electron::v0)
        .def("vx",  &mcmc::electron::vx)
        .def("vy",  &mcmc::electron::vy)
        .def("vz",  &mcmc::electron::vz)
        .def("v",  &mcmc::electron::v)
        .def("gamma", &mcmc::electron::gamma);

    py::class_<mcmc::photonBucket>(m, "photonBucket" )
        .def(py::init<>())
        .def("size", &mcmc::photonBucket::size)
        .def("replace", &mcmc::photonBucket::replace)
        .def("get", &mcmc::photonBucket::get)
        .def("push_back", &mcmc::photonBucket::push_back);



        // -------------------------------------------------- 
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
        })
        .def("__getitem__", [](const mcmc::vMesh &s, uint64_t i) {
                return s.__getitem__(i);
        })
        .def("__setitem__", [](mcmc::vMesh &s, uint64_t i, vblock_t v) {
                s.__setitem__(i, v);
        })
        // i,j,k indexing based interface
        .def("__getitem__", [](const mcmc::vMesh &s, py::tuple indx) {
                size_t i = indx[0].cast<size_t>();
                size_t j = indx[1].cast<size_t>();
                size_t k = indx[2].cast<size_t>();
                return s.__getitem2__( i,j,k );
        })
        .def("__setitem__", [](mcmc::vMesh &s, py::tuple indx, vblock_t v) {
                size_t i = indx[0].cast<size_t>();
                size_t j = indx[1].cast<size_t>();
                size_t k = indx[2].cast<size_t>();
                return s.__setitem2__( i,j,k, v);
        })

        // .def("__setitem__", [](mcmc::vMesh &s, uint64_t i, vblock_t v) {
        //         s.__setitem__(i, v);
        // })
        */


}
