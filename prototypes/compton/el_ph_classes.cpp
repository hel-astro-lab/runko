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

#include <Eigen/Dense>
using namespace Eigen;


// --------------------------------------------------
// #include "definitions.h"
const double pi = M_PI;
typedef std::array<double, 3> vec;

using std::sin;
using std::cos;
using std::exp;
using std::acos;
using std::asin;
using std::pow;
using std::sqrt;
using std::log;



// --------------------------------------------------
namespace el_ph {



    class photon {

        public:
            std::array<double, 8> data;
            photon( double, double, double, double, double, double, double, double );

            /// photon four-velocity components
            const double hv(){ return data[0]; };
            const double vx(){ return data[1]; };
            const double vy(){ return data[2]; };
            const double vz(){ return data[3]; };
            const double x(){ return data[4]; };
            const double y(){ return data[5]; };
            const double z(){ return data[6]; };
            const double w(){ return data[7]; };

            /// Location 
            // const double x(){ return data[4]; };
            // const double y(){ return data[5]; };
            // const double z(){ return data[6]; };

    };


    photon::photon(double E, double vx, double vy, double vz, double x, double y, double z, double w) {
        this->data = {{E, vx, vy, vz, x, y, z, w}};
    };



    // --------------------------------------------------
    class photonBucket {

        /// size of the bucket (in terms of photons)
        size_t nPhotons = 0;

        /// Photon container
        std::vector<std::array<double, 8>> bucket;

        public:
            void push_back( photon ph );

            //const size_t size( ) { return this->nPhotons; };
            const size_t size( ) { return bucket.size(); };

            void replace( const size_t indx, photon ph );

            photon get( const size_t indx );
            
            std::array<double, 8> get_data( const size_t indx );

            void resize(const size_t N);

            void swap(const size_t i1, const size_t i2);

            std::vector<double> vecE();  
            std::vector<double> vecZ();  
    };

    /// Append photon into bucket
    void photonBucket::push_back( photon ph ) {
        bucket.push_back( ph.data );
        nPhotons++;
    };

    /// Replace i:th photon with the given one
    void photonBucket::replace( const size_t indx, photon ph ) {
        bucket[indx] = ph.data;
    };

    /// Swap i1 and i2 photons in the bucket
    void photonBucket::swap(const size_t i1, const size_t i2) {
        std::swap( bucket[i1], bucket[i2] );
    };

    /// return i:th photon data
    photon photonBucket::get( const size_t indx ) {
        std::array<double, 8> data = bucket[indx];
        photon ph(data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]);

        return ph;
    };

    /// return i:th photon raw data
    std::array<double, 8> photonBucket::get_data( const size_t indx ) {
        return bucket[indx];
    };

    /// Resize bucket
    void photonBucket::resize(const size_t N) {
        // TODO error check if N < N
        bucket.resize(N);
        nPhotons = N;

        return;
    };

    /// collect energy components of velocity
    std::vector<double> photonBucket::vecE() {
        std::vector<double> ve;
        ve.resize( size() );

        for (size_t i=0; i<size(); i++) {
            auto d = bucket[i];
            ve[i] = d[0];
        }
        return ve;
    };

    /// collect z components of velocity
    std::vector<double> photonBucket::vecZ() {
        std::vector<double> vz;
        vz.resize( size() );

        for (size_t i=0; i<size(); i++) {
            auto d = bucket[i];
            vz[i] = d[3];
        }
        return vz;
    };




    /*
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
                          pow(vx(), 2.0) 
                        + pow(vy(), 2.0) 
                        + pow(vz(), 2.0) 
                                ); };

            /// gamma factor
            const double gamma() { return 1.0 / std::sqrt( 1.0 - std::pow( v(), 2.0 ) ); };

            /// proper momentum
            const double p() { return gamma() * v(); };

    };

    electron::electron(double v0, double vx, double vy, double vz) {
        this->data = {{v0, vx, vy, vz}};
    };
    */


    /// Electron with all velocity transformations build-in
    // Uses Eigen vectors internally for performance
    class electron {

        public:

        /// spatial components of the four-velocity
        Vector3d fvel;

        /// load from vector
        void loadFvel(Vector3d u) { fvel = u; };

        /// Load from components
        void loadFvelComponents(double ux, double uy, double uz) { 
                fvel(0) = ux;
                fvel(1) = uy;
                fvel(2) = uz;
        };

        /// Load from velocity
        // u = v*gamma = v/sqrt(1-v^2)
        void loadVel(Vector3d v) {
            fvel = v / sqrt(1.0 - v.dot(v));
        };

        void loadVelComponents(double vx, double vy, double vz) {
            Vector3d v(vx, vy, vz);
            loadVel(v);    
        };

/*        // return 3d components of the electron 4-velocity
        const double vx(){ return fvel(0); };
        const double vy(){ return fvel(1); };
        const double vz(){ return fvel(2); };
*/        

        /// gamma factor
        const double gamma() { return sqrt(1.0 + fvel.dot(fvel) ); };
            
        /// beta = v/c = sqrt(1-1/gamma^2)
        const double beta() { return sqrt(1.0 - 1.0/(gamma()*gamma())); };

        /// p = sqrt(gamma^2-1)
        const double pel() { return sqrt(gamma()*gamma()-1.0); };

        /// coordinate velocity vector v = u/gamma
        Vector3d vel() { return fvel/gamma(); };

        /// vmod = sqrt(vx^2 + vy^2 + vz^2)
        const double vmod() { return sqrt(vel().dot(vel())); };
        
        // unit vector in the electron direction 
        Vector3d beta0() { return vel()/beta(); };


    };



    // --------------------------------------------------
    class electronBucket {

        /// size of the bucket (number of electrons)
        size_t nElectrons = 0;

        /// Electron container
//        std::vector<std::array<double, 3>> bucket;
        std::vector<Vector3d> bucket;

        public:
            void push_back( electron e );
            const size_t size( ) { return this->nElectrons; };

            void replace( size_t indx, electron e );

            electron get( size_t indx);

    };

    void electronBucket::push_back( electron e ) {
        bucket.push_back( e.vel() );
        nElectrons++;
    };

    void electronBucket::replace( size_t indx, electron e ) {
        bucket[indx] = e.vel();
    };

    electron electronBucket::get( size_t indx ) {
        //std::array<double, 3> data = bucket[indx];
        electron e;
        Vector3d vel = bucket[indx];
        e.loadVel(vel);

        return e;
    };

}






// --------------------------------------------------
PYBIND11_MODULE(el_ph, m) {


    py::class_<el_ph::photon>(m, "photon" )
        .def(py::init<double, double, double, double, double, double, double, double >())
        .def_readwrite("data",      &el_ph::photon::data)
        .def("hv",  &el_ph::photon::hv)
        .def("vx",  &el_ph::photon::vx)
        .def("vy",  &el_ph::photon::vy)
        .def("vz",  &el_ph::photon::vz)
        .def("x",  &el_ph::photon::x)
        .def("y",  &el_ph::photon::y)
        .def("z",  &el_ph::photon::z)
        .def("w",  &el_ph::photon::w);

    py::class_<el_ph::electron>(m, "electron" )
        .def(py::init<>())
        .def("vx",   [](el_ph::electron &e){return e.vel()(0);})
        .def("vy",   [](el_ph::electron &e){return e.vel()(1);})
        .def("vz",   [](el_ph::electron &e){return e.vel()(2);})
        .def("loadVelComponents",  &el_ph::electron::loadVelComponents)
        .def("loadFvelComponents", &el_ph::electron::loadFvelComponents)
        .def("beta",  &el_ph::electron::beta)
        .def("gamma", &el_ph::electron::gamma)
        .def("vmod", &el_ph::electron::vmod)
        .def("beta0", &el_ph::electron::beta0);

    py::class_<el_ph::electronBucket>(m, "electronBucket" )
        .def(py::init<>())
        .def("size",      &el_ph::electronBucket::size)
        .def("replace",   &el_ph::electronBucket::replace)
        .def("get",       &el_ph::electronBucket::get)
        .def("push_back", &el_ph::electronBucket::push_back);

    py::class_<el_ph::photonBucket>(m, "photonBucket" )
        .def(py::init<>())
        .def("size",      &el_ph::photonBucket::size)
        .def("replace",   &el_ph::photonBucket::replace)
        .def("get",       &el_ph::photonBucket::get)
        .def("vecE",      &el_ph::photonBucket::vecE)
        .def("vecZ",      &el_ph::photonBucket::vecZ)
        .def("push_back", &el_ph::photonBucket::push_back);



}
