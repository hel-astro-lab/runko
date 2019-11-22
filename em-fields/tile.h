#pragma once

#include <vector>
#include <mpi4cpp/mpi.h>

#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../tools/mesh.h"
#include "../tools/rotator.h"
#include "../definitions.h"


namespace fields {
  namespace mpi = mpi4cpp::mpi;


/// Yee lattice of plasma quantities
class YeeLattice 
//  public std::enable_shared_from_this<YeeLattice>
{

  public:

  int Nx;
  int Ny;
  int Nz;

  /// Electric field 
  toolbox::Mesh<float_t, 3> ex;
  toolbox::Mesh<float_t, 3> ey;
  toolbox::Mesh<float_t, 3> ez;
  
  /// Magnetic field 
  toolbox::Mesh<float_t, 3> bx;
  toolbox::Mesh<float_t, 3> by;
  toolbox::Mesh<float_t, 3> bz;
    
  /// Charge density
  toolbox::Mesh<float_t, 1> rho;

  /// Current vector 
  toolbox::Mesh<float_t, 3> jx;
  toolbox::Mesh<float_t, 3> jy;
  toolbox::Mesh<float_t, 3> jz;

  /// temporary current
  toolbox::Mesh<float_t, 3> jx1;
  toolbox::Mesh<float_t, 3> jy1;
  toolbox::Mesh<float_t, 3> jz1;


  // default empty constructor
  YeeLattice() 
  {
    std::cout << "YeeL ctor empty\n";
  };


  // real initializer constructor
  YeeLattice(int Nx, int Ny, int Nz) : 
    Nx(Nx), Ny(Ny), Nz(Nz),
    ex(Nx, Ny, Nz),
    ey(Nx, Ny, Nz),
    ez(Nx, Ny, Nz),
    bx(Nx, Ny, Nz),
    by(Nx, Ny, Nz),
    bz(Nx, Ny, Nz),
    rho(Nx, Ny, Nz),
    jx(Nx, Ny, Nz),
    jy(Nx, Ny, Nz),
    jz(Nx, Ny, Nz),
    jx1(Nx, Ny, Nz),
    jy1(Nx, Ny, Nz),
    jz1(Nx, Ny, Nz)
    { 
      std::cout << "YeeL ctor: (" << Nx << "," << Ny << "," << Nz << ")\n";
    }

  // copy ctor
  YeeLattice(YeeLattice& other) :
    Nx(other.Nx),
    Ny(other.Ny),
    Nz(other.Nz),
    ex(other.ex),
    ey(other.ey),
    ez(other.ez),
    bx(other.bx),
    by(other.by),
    bz(other.bz),
    rho(other.rho),
    jx(other.jx),
    jy(other.jy),
    jz(other.jz),
    jx1(other.jx1),
    jy1(other.jy1),
    jz1(other.jz1)
  {
    std::cout << "Yee copy ctor\n";
  }

  YeeLattice(const YeeLattice& other) :
    Nx(other.Nx),
    Ny(other.Ny),
    Nz(other.Nz),
    ex(other.ex),
    ey(other.ey),
    ez(other.ez),
    bx(other.bx),
    by(other.by),
    bz(other.bz),
    rho(other.rho),
    jx(other.jx),
    jy(other.jy),
    jz(other.jz),
    jx1(other.jx1),
    jy1(other.jy1),
    jz1(other.jz1)
  {
    std::cout << "const Yee copy ctor\n";
  }

  // move constructor
  YeeLattice(YeeLattice&& other)
      : YeeLattice() // initialize via default constructor, C++11 only
  {
    swap(*this, other);
  }


  // public swap for efficient memory management
  friend void swap(YeeLattice& first, YeeLattice& second)
  {
    using std::swap;
    swap(first.Nx,  second.Nx);
    swap(first.Ny,  second.Ny);
    swap(first.Nz,  second.Nz);
    swap(first.ex , second.ex);
    swap(first.ey , second.ey);
    swap(first.ez , second.ez);
    swap(first.bx , second.bx);
    swap(first.by , second.by);
    swap(first.bz , second.bz);
    swap(first.rho, second.rho);
    swap(first.jx , second.jx);
    swap(first.jy , second.jy);
    swap(first.jz , second.jz);
    swap(first.jx1, second.jx1);
    swap(first.jy1, second.jy1);
    swap(first.jz1, second.jz1);
  }

  // copy-and-swap algorithm
  YeeLattice& operator=(YeeLattice other) 
  {
    swap(*this, other);
    return *this;
  }

  //virtual ~YeeLattice() = default;
  ~YeeLattice() {
    std::cout << "~YeeLattice\n";
  }


  

};


/// Lattice to hold plasma moment values of separate species
class PlasmaMomentLattice {

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  // density
  toolbox::Mesh<float_t, 0> rho;

  /// energy density
  toolbox::Mesh<float_t, 0> edens;

  /// temperature
  toolbox::Mesh<float_t, 0> temp;

  /// (mean) flow speed
  toolbox::Mesh<float_t, 0> Vx;
  toolbox::Mesh<float_t, 0> Vy;
  toolbox::Mesh<float_t, 0> Vz;

  /// momentum density
  toolbox::Mesh<float_t, 0> momx;
  toolbox::Mesh<float_t, 0> momy;
  toolbox::Mesh<float_t, 0> momz;

  /// pressure
  toolbox::Mesh<float_t, 0> pressx;
  toolbox::Mesh<float_t, 0> pressy;
  toolbox::Mesh<float_t, 0> pressz;

  /// shear
  toolbox::Mesh<float_t, 0> shearxy;
  toolbox::Mesh<float_t, 0> shearxz;
  toolbox::Mesh<float_t, 0> shearyz;


  /// constructor; initializes Mesh objects at the same time
  PlasmaMomentLattice(size_t Nx_in, size_t Ny_in, size_t Nz_in) : 
    Nx(Nx_in), Ny(Ny_in), Nz(Nz_in),

    rho(Nx, Ny, Nz),
    edens(Nx, Ny, Nz),
    temp(Nx, Ny, Nz),
    Vx( Nx, Ny, Nz ),
    Vy( Nx, Ny, Nz ),
    Vz( Nx, Ny, Nz ),
    momx( Nx, Ny, Nz ),
    momy( Nx, Ny, Nz ),
    momz( Nx, Ny, Nz ),
    pressx( Nx, Ny, Nz ),
    pressy( Nx, Ny, Nz ),
    pressz( Nx, Ny, Nz ),
    shearxy( Nx, Ny, Nz ),
    shearxz( Nx, Ny, Nz ),
    shearyz( Nx, Ny, Nz )
  { }

  //virtual ~PlasmaMomentLattice() = default;

  // clear all internal storages
  void clear() 
  {
    rho    .clear();
    edens  .clear();
    temp   .clear();
    Vx     .clear();   
    Vy     .clear();  
    Vz     .clear();     
    momx   .clear();      
    momy   .clear();       
    momz   .clear();       
    pressx .clear();        
    pressy .clear();       
    pressz .clear();       
    shearxy.clear();
    shearxz.clear();
    shearyz.clear();  
  }
};




/*! \brief General Plasma tile for solving Maxwell's equations
 *
 * Internally everything for computations are stored in 
 * staggered Yee lattice.
 *
 * Analysis/reduced data is in PlasmaMomentLattice.
 */
template<std::size_t D>
class Tile : 
  virtual public corgi::Tile<D>,
  public std::enable_shared_from_this<Tile<D>>
{

  public:

  /// size of grid inside tile
  std::array<size_t, 3> mesh_lengths;

  /// Yee lattice of plasma quantities (with 1 timestep)
  //toolbox::Rotator<YeeLattice, 1> yee;
  YeeLattice yee;
  //std::vector<YeeLattice> yee;

  /// species-specific analysis results
  std::vector<PlasmaMomentLattice> analysis;

  /// explicitly show that we import tile limits from base class 
  using corgi::Tile<D>::mins;
  using corgi::Tile<D>::maxs;


  /// CFL number (corresponds to simulation light speed c)
  double_t cfl;

  /// grid size (assuming cubical cells)
  //double_t dx = 1.0;

  //--------------------------------------------------
  // constructor with internal mesh dimensions
  Tile(size_t nx, size_t ny, size_t nz) :
    mesh_lengths {{nx, ny, nz}},
    yee((int)nx, (int)ny, (int)nz)
  {
    if (D == 1) assert(ny == 1 && nz == 1);
    if (D == 2) assert(nz == 1);

    // initialize one Yee lattice into the grid 
    // TODO: into data rotator?
    //std::cout<<"len0 of yee vec:" << yee.size() << "\n";
    //add_yee_lattice();
    //std::cout<<"len1 of yee vec:" << yee.size() << "\n";
  }

  // avoid copies; TODO: is this needed?
  Tile(Tile& ) = delete;

  //~Tile() = default;
  ~Tile() {
    std::cout << "~Tile\n";
  }

  //--------------------------------------------------

  virtual void update_boundaries(  corgi::Grid<D>& grid);

  virtual void exchange_currents(  corgi::Grid<D>& grid);

  virtual void deposit_current();

  virtual YeeLattice& get_yee(size_t i=0);

  virtual YeeLattice& get_yee2();

  virtual void set_yee(YeeLattice& v);

  virtual std::shared_ptr<YeeLattice> get_yeeptr();

  virtual PlasmaMomentLattice& get_analysis(size_t i);

  //FIXME: naming is wrong because of pybind: it can not deduce correct function
  //       with only constness difference.
  virtual const YeeLattice& get_const_yee(size_t i=0) const; 

  virtual const PlasmaMomentLattice& get_const_analysis(size_t i) const;

  void cycle_yee();

  void cycle_current();

  void clear_current();

  void add_yee_lattice();

  void add_analysis_species();

  virtual std::vector<mpi::request> 
  send_data( mpi::communicator&, int orig, int mode, int tag) override;

  virtual std::vector<mpi::request> 
  recv_data( mpi::communicator&, int dest, int mode, int tag) override;
};



} // end of namespace fields





