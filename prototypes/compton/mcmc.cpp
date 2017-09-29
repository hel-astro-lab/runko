#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <fmt/format.h>
#include <fmt/format.cc>
#include <fmt/string.h>
#include <fmt/ostream.h>

#include <vector>
#include <cmath>
#include <random>

#include <math.h>

// #include "definitions.h"
const double pi = M_PI;

typedef std::array<double, 3> vec;


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

            /// Location 
            // const double x(){ return data[4]; };
            // const double y(){ return data[5]; };
            // const double z(){ return data[6]; };

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

            void replace( const size_t indx, photon ph );

            photon get( const size_t indx );
            
            std::array<double, 4> get_data( const size_t indx );

            void resize(const size_t N);

            void swap(const size_t i1, const size_t i2);

    };

    void photonBucket::push_back( photon ph ) {
        bucket.push_back( ph.data );
        nPhotons++;
    };

    void photonBucket::replace( const size_t indx, photon ph ) {
        bucket[indx] = ph.data;
    };

    void photonBucket::swap(const size_t i1, const size_t i2) {
        std::swap( bucket[i1], bucket[i2] );
    };

    photon photonBucket::get( const size_t indx ) {
        std::array<double, 4> data = bucket[indx];
        photon ph(data[0], data[1], data[2], data[3]);

        return ph;
    };

    std::array<double, 4> photonBucket::get_data( const size_t indx ) {
        return bucket[indx];
    };

    void photonBucket::resize(const size_t N) {
        // TODO error check if N < N
        bucket.resize(N);
        nPhotons = N;

        return;
    };



    class Slab {

        /// photon container
        photonBucket bucket;

        /// Box sizes
        double xmin, xmax, ymin, ymax, zmin, zmax;

        /// RNG seed
        uint32_t rngSeed = 1; 

        /// RNG engine (Mersenne twister)
        std::mt19937 rng;

        std::uniform_real_distribution<double> randPhi{0.0, 2.0*pi};
        std::uniform_real_distribution<double> randmu{0.0, 1.0};

        /// Random spherical direction (r, theta, phi)
        vec randSph() {
            return {{ 1.0, std::acos( randmu(rng) ), randPhi(rng) }};
        };

        /// Random direction in cartesian (vx, vy, vx) coordinates
        // done via spherical coordinates (r, theta, phi) and then transforming back
        vec randHalfSphere() {
            vec sphDirs = randSph();
            vec ret = 
            {{
                sphDirs[0] * std::sin(sphDirs[1]) * std::cos(sphDirs[2]),
                sphDirs[0] * std::sin(sphDirs[1]) * std::sin(sphDirs[2]),
                sphDirs[0] * std::cos(sphDirs[1])
            }};

            return ret;
        };



        public:

            /// Simulation time step (in units of c)
            double dt = 0.1;

            /// Slab height
            double height = 1.0;

            /// electron number density
            double ne = 0.0;

            /// Thomson optical depth
            double tau = 0.0;


            /// location containers
            std::vector<double> xloc, yloc, zloc;


            /// Constructor
            Slab(photonBucket b) {
                this->bucket = b;

                // prepare location containers
                xloc.resize( bucket.size() );
                yloc.resize( bucket.size() );
                zloc.resize( bucket.size() );

                // Finally seed the rng
                rng.seed( rngSeed );
            };
              
            /// Number of photons in the slab
            const size_t size( ) {return this->bucket.size(); };


            /// Set slab dimensions; z is implicitly assumed as the height
            void set_dimensions(double _xmin, double _xmax,
                                double _ymin, double _ymax,
                                double _zmin, double _zmax) {

                this->xmin = _xmin;
                this->xmax = _xmax;
                this->ymin = _ymin;
                this->ymax = _ymax;
                this->zmin = _zmin;
                this->zmax = _zmax;

                this->height = zmax - zmin;
            };


            /// Set number density and compute Thomson optical depth based on it
            void set_numberDensity(double _ne) {
                this->ne = _ne;

                // compute Thomson tau = sigma_T * n_e * H
                tau = 1.0 * ne * height;
            };


            /// Step in time performing the full radiation interactions
            void step() {

                push();
                // wrap();
                // emergingFlux();
                
                // check_scatter();
                // scatter();
                // inject();

            };


            /// Push photons
            // TODO: Properly vectorize although this probably implicitly works 
            // already on compiler level
            void push() {

                size_t N = this->size();
                std::vector<double> vx, vy, vz;
                vx.resize(N);
                vy.resize(N);
                vz.resize(N);

                // get velocities from bucket
                for (size_t i=0; i<N; i++) {
                    auto vel = this->bucket.get_data(i);
                    vx[i] = vel[1];
                    vy[i] = vel[2];
                    vz[i] = vel[3];
                }

                // step forward in time
                for (size_t i=0; i<N; i++) {
                    xloc[i] += vx[i]*dt;
                    yloc[i] += vy[i]*dt;
                    zloc[i] += vz[i]*dt;
                }
            };


            /// Inject more from the floor
            void inject(double flux) {

                // size of the floor
                double area = (xmax-xmin)*(ymax-ymin);

                // how many to inject based on the flux
                size_t Ninj = (size_t)flux*area*dt;

                // fmt::print("Injecting {} photons...\n", Ninj);

                // resize beforehand 
                size_t Ns = size();
                bucket.resize(Ns + Ninj);
                xloc.resize(Ns + Ninj);
                yloc.resize(Ns + Ninj);
                zloc.resize(Ns + Ninj);


                for (size_t i=0; i<Ninj; i++) {

                    // create random photon
                    double E = 0.01;
                    vec dir = randHalfSphere();
                    photon ph(E, dir[0], dir[1], dir[2]);

                    bucket.replace( Ns + i, ph );

                    // set location
                    xloc[Ns + i] = 0.0; // TODO random loc
                    yloc[Ns + i] = 0.0; // TODO random loc
                    zloc[Ns + i] = 0.0; // bottom
                }

            };


            /// Scrape photons that are overflowing from the slab
            void scrape() {

                size_t i = 0;
                size_t Nspills = 0;

                while ( i<size()-Nspills) {

                    if (zloc[i] >= height) {
                        // swamp everything to the end and remove later
                        std::swap( xloc[i], xloc[size() - Nspills] );
                        std::swap( yloc[i], yloc[size() - Nspills] );
                        std::swap( zloc[i], zloc[size() - Nspills] );
                        bucket.swap( i, size()-Nspills );

                        Nspills++;
                    } else { 
                        i++;
                    }
                }

                // collecting spills
                fmt::print("Scraping spills: {}\n", Nspills);
                for (size_t i=size()-Nspills; i<size(); i++) {
                    photon ph = bucket.get(i);
                    fmt::print("  spill {}\n", ph.hv() );
                }

                // and finally remove 
                xloc.resize(size() - Nspills); 
                yloc.resize(size() - Nspills); 
                zloc.resize(size() - Nspills); 
                bucket.resize(size() - Nspills);

            };

            /// inject everything from point in bottom
            void floor() {
                for (size_t i=0; i<xloc.size(); i++) {
                    zloc[i] = zmin;

                    // center of floor
                    xloc[i] = 0.0;
                    yloc[i] = 0.0;
                }
            };


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
        .def("v0",    &mcmc::electron::v0)
        .def("vx",    &mcmc::electron::vx)
        .def("vy",    &mcmc::electron::vy)
        .def("vz",    &mcmc::electron::vz)
        .def("v",     &mcmc::electron::v)
        .def("gamma", &mcmc::electron::gamma);

    py::class_<mcmc::photonBucket>(m, "photonBucket" )
        .def(py::init<>())
        .def("size",      &mcmc::photonBucket::size)
        .def("replace",   &mcmc::photonBucket::replace)
        .def("get",       &mcmc::photonBucket::get)
        .def("push_back", &mcmc::photonBucket::push_back);


    py::class_<mcmc::Slab>(m, "Slab" )
        .def(py::init< mcmc::photonBucket >())
        .def_readwrite("xloc",    &mcmc::Slab::xloc)
        .def_readwrite("yloc",    &mcmc::Slab::yloc)
        .def_readwrite("zloc",    &mcmc::Slab::zloc)
        .def_readwrite("tau",     &mcmc::Slab::tau)
        .def_readwrite("height",  &mcmc::Slab::height)
        .def_readwrite("ne",      &mcmc::Slab::ne)
        .def("size",              &mcmc::Slab::size)
        .def("push",              &mcmc::Slab::push)
        .def("inject",            &mcmc::Slab::inject)
        .def("set_dimensions",    &mcmc::Slab::set_dimensions)
        .def("set_numberDensity", &mcmc::Slab::set_numberDensity)
        .def("scrape",            &mcmc::Slab::scrape)
        .def("floor",             &mcmc::Slab::floor);



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
