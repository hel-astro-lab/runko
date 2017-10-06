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

        /// gamma factor
        const double gamma() { return sqrt(1.0 + fvel.dot(fvel) ); };
            
        /// beta = v/c = sqrt(1-1/gamma^2)
        const double beta() { return sqrt(1.0 - 1.0/(gamma()*gamma())); };

        /// coordinate velocity vector v = u/gamma
        Vector3d vel() { return fvel/gamma(); };

    };



    // --------------------------------------------------
    class photonBucket {

        /// size of the bucket (in terms of photons)
        size_t nPhotons = 0;

        /// Photon container
        std::vector<std::array<double, 4>> bucket;

        public:
            void push_back( photon ph );

            //const size_t size( ) { return this->nPhotons; };
            const size_t size( ) { return bucket.size(); };

            void replace( const size_t indx, photon ph );

            photon get( const size_t indx );
            
            std::array<double, 4> get_data( const size_t indx );

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
        std::array<double, 4> data = bucket[indx];
        photon ph(data[0], data[1], data[2], data[3]);

        return ph;
    };

    /// return i:th photon raw data
    std::array<double, 4> photonBucket::get_data( const size_t indx ) {
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




    //-------------------------------------------------- 
    class Slab {


        /// Box sizes
        double xmin, xmax, ymin, ymax, zmin, zmax;

        /// RNG seed
        uint32_t rngSeed = 1; 

        /// RNG engine (Mersenne twister)
        std::mt19937 rng;

        /// Ready-made distribution for flat variables in [0, 2pi[
        std::uniform_real_distribution<double> randPhi{0.0, 2.0*pi};
          
        /// Ready-made distribution for flat variables in [0,1[
        std::uniform_real_distribution<double> randmu{0.0, 1.0};

        /// Random spherical direction (r, theta, phi) with Lamberts law
        vec randSph() {
            return {{ 1.0, std::acos(std::pow(randmu(rng),0.5)), randPhi(rng) }};
        };

        /// Random direction in Cartesian (vx, vy, vx) coordinates
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



        /// Draw samples from Planck function using series sampling
        double Planck(double kT) {

            double z1 = randmu(rng);
            double z2 = randmu(rng);
            double z3 = randmu(rng);
            double x = -1.0*log(z1*z2*z3);

            double j = 1.0;
            double a = 1.0;
            z1 = randmu(rng);
            
            // 90/pi^4 (rounded down to get finite loop)
            while (1.08232*z1 > a) {
                j += 1.0;
                a += 1.0 / std::pow(j, 4.0);
            }
            return kT * x/j;
        };


        /// Sampling relativistic Maxwellian using the Sobol method
        // u = v*gamma
        // Ref: Sobol 1976, Pozdnyakov+ 1977
        double relMaxwellianVel(double Te) {
            double x4 = randmu(rng);
            double x5 = randmu(rng);
            double x6 = randmu(rng);
            double x7 = randmu(rng);

            double u = -Te*log(x4*x5*x6);
            double n = -Te*log(x4*x5*x6*x7);

            if (n*n - u*u < 1.0) { return relMaxwellianVel(Te); };

            return u;
        };

        /// Sampling from non-relativistic Maxwellian with rejection sampling
        // NOTE: these are coordinate velocities originally
        /*
        double MaxwellianVel(double Te) {
            double vmin = -5.0*Te; 
            double vmax =  5.0*Te; 
            double vf = vmin + (vmax-vmin)*randmu(rng);

            double f = vf*vf*std::exp(-(vf*vf)/(2.0*Te));
            double x = randmu(rng);

            if (x > f) { return MaxwellianVel(Te); };

            return vf;
        };
        */

        double MaxwellianVel(double Te) {
            double z1 = randmu(rng);
            double z2 = randmu(rng);
            double y = -2.0*log(z1);

            const double g = 1.647;

            if ( g*g*y < (z1/z2)*(z1/z2) ) {
                return MaxwellianVel(Te);
            }
            
            return sqrt(Te*y/2.0);
        };


        // Standard distribution sample using Box-Muller algorithm
        double BoxMuller(double Te) {
            return sqrt(-2.0*Te*log( randmu(rng) ));
        };



        /// Isotropic velocity components
        // using u = abs(u_i) to get cartesian (ux, uy, uz)
        Vector3d velXYZ(double u ) {
            double x1 = randmu(rng);
            double x2 = randmu(rng);

            Vector3d ret( 
                    u*(2.0*x1 - 1.0),
                    2.0*u*sqrt(x1*(1.0-x1))*cos(2.0*pi*x2),
                    2.0*u*sqrt(x1*(1.0-x1))*sin(2.0*pi*x2)
                        );
            return ret;
        };
        

        public:
        /// General boosted non-rel/rel Maxwellian
        // Te: electron temperature in m_e c^2
        // G:  Bulk Lorentz vector
        //
        // Ref: Zenitani 2015
        // NOTE: these are proper velocities and gamma = sqrt(1+u^2)
        Vector3d boostedMaxwellian(double Te, vec G) {

            double u;
            if (Te > 0.2) { // relativistic
                double u = relMaxwellianVel(Te);
                  
                // get isotropic velocity components
                Vector3d ui = velXYZ(u);

                return ui;

            } else { // non-relativistic
                double v1 = BoxMuller(Te);
                double v2 = BoxMuller(Te);

                double x1 = randPhi(rng);
                double x2 = randPhi(rng);

                Vector3d vi(
                        v1*sin(x1),
                        v1*cos(x1),
                        v2*sin(x2)
                           );

                // add bulk motion
                // vx += V0

                // double gamma = 1.0/sqrt(1.0 + vi.dot(vi));

                return vi;
            }
                

            /* TODO add general bulk transformation
            // check if bulk velocity
            if (G[0] == 0.0 && G[1] == 0.0 && G[2] == 0.0) {
                return ui;
            }

            // next boost in X dir; TODO generalize
            vec beta = 
            {{ 
                1.0/sqrt(1.0 + G[0]*G[0]),
                1.0/sqrt(1.0 + G[1]*G[1]),
                1.0/sqrt(1.0 + G[2]*G[2])
            }};


            double x8 = randmu(rng);
            if (-beta[0]*ui[0] > x8) { ui[0] = -ui[0]; };
            ui[0] = G[0]*(ui[0] + beta[0]*std::sqrt(1.0 + u*u));
            // u = std::sqrt(ui[0]*ui[0] + ui[1]*ui[1] + ui[2]*ui[2]);

            return ui;
            */
        };


        /// Draw electron sample from boosted Maxwellian distribution
        electron sampleBoostedMaxw(double Te, vec G) {
            Vector3d ui = boostedMaxwellian(Te, G);
            electron el;
            el.loadFvel(ui);
            return el;
        };


        /// Compton scattering using Sobol's algorithm
        std::pair<photon, electron> comptonScatter(photon ph, electron el) {

            // fmt::print("v: {}\n", el.v());
            // fmt::print("E: {}\n", ph.hv());

            // Vector3d ve( el.vx(), el.vy(), el.vz() );
            // Vector3d beta = ve.normalized();

            Vector3d omeg(ph.vx(), ph.vy(), ph.vz());

            double theta = acos( el.vel().dot(omeg) );

            // Create base vectors and matrix
            //-------------------------------------------------- 
            Vector3d kvec(0.0, -1.0, 0.0);                         // k
            Vector3d jvec = ( el.vel().cross(omeg) ).normalized(); // j
            Vector3d ivec = ( kvec.cross(jvec) ).normalized();     // i

            Matrix3d M;
            M << ivec, jvec, kvec;

            // --------------------------------------------------
            
            // unit vector of electron in scattering coordinates
            Vector3d v0 = M.transpose() * el.vel(); // rotating electron velocity 
                                                    // to scattering plane (i,k)
            double mu = v0(0)*sin(theta) + v0(2)*cos(theta);
            double rho = sqrt( pow(v0(0),2) + pow(v0(1),2) );

            // Compton parameter
            // here hv = hv/m_e c^2
            double y = ph.hv() * el.gamma() * (1.0 - mu*el.beta() );
            

            // Additional scattering angles (v0, w0, t0) that define a frame of reference
            Vector3d w0(v0(1)/rho,      -v0(0)/rho,        0.0);
            Vector3d t0(v0(0)*v0(2)/rho, v0(1)*v0(2)/rho, -rho);


            // --------------------------------------------------
            // scatter
            double OOp, z1, z2, z3, mup, phip, yp, Y; 
            while (true) {
                z1 = randmu(rng);
                z2 = randmu(rng);
                z3 = randmu(rng);

                mup  = (el.beta() + 2.0*z1 - 1.0)/(1.0 + el.beta()*(2.0*z1 - 1.0));
                phip = 2.0*pi*z2;

                OOp = mu*mup - sqrt(1.0-mup*mup) * (rho*sin(phip)*cos(theta) 
                      - (1.0/rho)*(v0(1)*cos(phip) + v0(0)*v0(2)*sin(phip))*sin(theta));

                yp = y/(1.0 + ph.hv()*(1.0 - OOp))/(el.gamma() * (1.0-mup*el.beta()));
                Y = yp/y + pow(yp/y,3) + pow(yp/y,2)*
                    ( pow(1.0/yp - 1.0/y, 2) - 2.0*( 1.0/yp - 1.0/y) );

                if (Y>2.0*z3) { break; };
            }

            // --------------------------------------------------
            // we have now scattered successfully
              
            // new energy
            double hvp = yp/( el.gamma()*(1.0 - mup*el.beta()) );

            // new direction from ijk base to xyz base
            Vector3d Op_ijk = mup*v0 
                          + sqrt(1.0-mup*mup)*( w0*cos(phip) + t0*sin(phip) );
            Vector3d Op = (M.transpose().inverse() * Op_ijk).normalized();


            // pack everything to classes and return
            photon phs(hvp, Op(0), Op(1), Op(2));

            return std::make_pair(phs, el);
        };


        public:
         
            /// photon container
            photonBucket bucket;

            /// Overflow bucket 
            photonBucket overflow;


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
                    double E = Planck(1.0);

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
                        // fmt::print("{} swapping {} to {} \n", i, yloc[i], yloc[size()-Nspills-1]);

                        std::swap( xloc[i], xloc[size() - Nspills - 1] );
                        std::swap( yloc[i], yloc[size() - Nspills - 1] );
                        std::swap( zloc[i], zloc[size() - Nspills - 1] );
                        bucket.swap( i, size()-Nspills - 1 );

                        Nspills++;
                    } else { 
                        i++;
                    }
                }

                // collecting spills
                // fmt::print("Scraping spills: {} // total size: {}\n", Nspills, size());
                for (size_t i=size()-Nspills; i<size(); i++) {
                    photon ph = bucket.get(i);
                    // fmt::print("  spill {}\n", ph.hv() );

                    overflow.push_back(ph);
                }

                // and finally remove 
                xloc.resize(size() - Nspills); 
                yloc.resize(size() - Nspills); 
                zloc.resize(size() - Nspills); 
                bucket.resize(size() - Nspills);

                //fmt::print("after size {}\n", size());
            };


            /// inject everything from point in bottom
            void floor() {
                for (size_t i=0; i<size(); i++) {
                    zloc[i] = zmin;

                    // center of floor
                    xloc[i] = 0.0;
                    yloc[i] = 0.0;
                }
            };


            /// Wrap into xy box bounded by [xmin, xmax] x [ymin, ymax]
            void wrap() {

                for (size_t i=0; i<size(); i++) {
                    if (xloc[i] < xmin) { xloc[i] += xmax; }
                    if (xloc[i] > xmax) { xloc[i] -= xmax; }

                    if (yloc[i] < ymin) { yloc[i] += ymax; }
                    if (yloc[i] > ymax) { yloc[i] -= ymax; }
                }

            };

            // Check the optical distance and then scatter
            void scatter(double Te) {

                double x, y, z, z0;

                double etau = std::exp(-dt/ne); // optical distance
                electron el; // target electron

                for (size_t i=0; i<size(); i++) {
                    z = zloc[i];

                    // scatter if e^-d\tau < rand()
                    z0 = randmu(rng);

                    if (etau < z0) {
                        photon ph = bucket.get(i);

                        // isotropic mono-energetic electron
                        // vec ve = velXYZ(0.8); // (beta)
                        // electron el(1.0, ve[0], ve[1], ve[2]);

                        // isotropic Maxwellian electrons
                        // fmt::print("sampling from Maxwellian...\n");
                        Vector3d ve = boostedMaxwellian(Te, {{0.0, 0.0, 0.0}});
                        el.loadFvel(ve);

                        // fmt::print("vx {} / vy {} / vz {}\n", el.vel()(0), el.vel()(1), el.vel()(2));
                        // fmt::print("target electron gamma: {} \n", el.gamma() );
                        // fmt::print("target electron beta: {} \n", el.beta() );

                        auto ret = comptonScatter(ph, el);
                        bucket.replace(i, ret.first );
                    }

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
        .def(py::init<>())
        .def("vx",   [](mcmc::electron &e){return e.vel()(0);})
        .def("vy",   [](mcmc::electron &e){return e.vel()(1);})
        .def("vz",   [](mcmc::electron &e){return e.vel()(2);})
        .def("loadVelComponents",  &mcmc::electron::loadVelComponents)
        .def("loadFvelComponents", &mcmc::electron::loadFvelComponents)
        .def("beta",  &mcmc::electron::beta)
        .def("gamma", &mcmc::electron::gamma);

    py::class_<mcmc::photonBucket>(m, "photonBucket" )
        .def(py::init<>())
        .def("size",      &mcmc::photonBucket::size)
        .def("replace",   &mcmc::photonBucket::replace)
        .def("get",       &mcmc::photonBucket::get)
        .def("vecE",      &mcmc::photonBucket::vecE)
        .def("vecZ",      &mcmc::photonBucket::vecZ)
        .def("push_back", &mcmc::photonBucket::push_back);


    py::class_<mcmc::Slab>(m, "Slab" )
        .def(py::init< mcmc::photonBucket >())
        .def_readwrite("xloc",    &mcmc::Slab::xloc)
        .def_readwrite("yloc",    &mcmc::Slab::yloc)
        .def_readwrite("zloc",    &mcmc::Slab::zloc)
        .def_readwrite("tau",     &mcmc::Slab::tau)
        .def_readwrite("height",  &mcmc::Slab::height)
        .def_readwrite("ne",      &mcmc::Slab::ne)
        .def_readonly("bucket",   &mcmc::Slab::bucket)
        .def_readonly("overflow", &mcmc::Slab::overflow)
        .def("size",              &mcmc::Slab::size)
        .def("push",              &mcmc::Slab::push)
        .def("inject",            &mcmc::Slab::inject)
        .def("set_dimensions",    &mcmc::Slab::set_dimensions)
        .def("set_numberDensity", &mcmc::Slab::set_numberDensity)
        .def("scrape",            &mcmc::Slab::scrape)
        .def("wrap",              &mcmc::Slab::wrap)
        .def("scatter",           &mcmc::Slab::scatter)
        .def("boostedMaxwellian", &mcmc::Slab::boostedMaxwellian)
        .def("sampleBoostedMaxw", &mcmc::Slab::sampleBoostedMaxw)
        .def("comptonScatter",    &mcmc::Slab::comptonScatter)
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
