#ifndef PRTCS_HPP
#define PRTCS_HPP


#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "cell.hpp"
#include "init.hpp"
#include <Eigen/Dense>


using namespace std;
using namespace dccrg;
using namespace Eigen;

typedef Array<double, 3, 2> Matrixd32;



class EB_Cube
{
public:

    // fixed size neighborhood array for the field tensors
    std::array<Matrixd32, 27> fields;

    // take grid and cell as pointers
    // build adjacent neighborhood and related functions around them
	template<
		class CellData,
		class Geometry
	> EB_Cube(
		dccrg::Dccrg<CellData, Geometry>& grid,
        uint64_t cell
	) {

        // first get cell itself (i.e. 0,0,0)
        //-------------------------------------------------- 
		auto* const cell_data = grid[cell];
    
        // TODO FIXME update to Yee staggeration
        fields[13](0,0) = cell_data->EY(0);
        fields[13](1,0) = cell_data->EY(1);
        fields[13](2,0) = cell_data->EY(2);

        fields[13](0,1) = cell_data->BY(0);
        fields[13](1,1) = cell_data->BY(1);
        fields[13](2,1) = cell_data->BY(2);


        // then neighboring cells
        //-------------------------------------------------- 
        const auto* const neighbors = grid.get_neighbors_of(cell);

        // TODO: start index depending on neighborhood;
        // now we just assume full Neumann neighbors and start from ijk=0
        // int ijk = 14;
        int ijk = 0;
        for (const auto neigh: *neighbors) {

            // cout << "hi from cube: " << neigh << " ijk:" << ijk << endl;
            // TODO: fixme 
            if(ijk == 13){ijk++;}; // skip (0,0,0)

            if (neigh == dccrg::error_cell){
                // TODO deal with boundary conditions

                fields[ijk] = Matrixd32::Zero();

            } else {
                auto* const cell_data = grid[neigh];

                fields[ijk](0,0) = cell_data->EY(0);
                fields[ijk](1,0) = cell_data->EY(1);
                fields[ijk](2,0) = cell_data->EY(2);

                fields[ijk](0,1) = cell_data->BY(0);
                fields[ijk](1,1) = cell_data->BY(1);
                fields[ijk](2,1) = cell_data->BY(2);
            }

            ijk++;
        } // end of neighborhood loop


        return;
    };

    Matrixd32& operator [](const std::size_t i){return this->fields[i];}
    const Matrixd32& operator [](const std::size_t i) const {return this->fields[i];}

    Matrixd32& F(int i, int j, int k)
    { 
        int ijk = (i+1) + (j+1)*3 + (k+1)*3*3;
        return this->fields[ijk];
    }
    
    double exY(int i, int j, int k)
    {
        return this->F(i,j,k)(0,0);
    }
    double eyY(int i, int j, int k)
    {
        return this->F(i,j,k)(1,0);
    }
    double ezY(int i, int j, int k)
    {
        return this->F(i,j,k)(2,0);
    }
    double bxY(int i, int j, int k)
    {
        return this->F(i,j,k)(0,1);
    }
    double byY(int i, int j, int k)
    {
        return this->F(i,j,k)(1,1);
    }
    double bzY(int i, int j, int k)
    {
        return this->F(i,j,k)(2,1);
    }


    /* basic trilinear interpolation (i.e. 3x linear interpolations) with the cube class
        x,y,z weights within the cells (0 means at i/j/k and 1 means i+1/j+1/k+1
    
        TODO: optimize this; now we call f000... elements each time whereas one time is enough
          only xd, yd, zd change
        TODO: staggered version is slooooow! And ugly.
    */

    /* XXX From TRISTAN-MP
        l=i+iy*(j-1)+iz*(k-1)

        f=ex(l,1,1)+ex(l-ix,1,1)+dx*(ex(l+ix,1,1)-ex(l-ix,1,1))

		f=f+dy*(ex(l+iy,1,1)+ex(l-ix+iy,1,1)+dx*(ex(l+ix+iy,1,1)-ex(l-ix &
		+iy,1,1))-f)

		g=ex(l+iz,1,1)+ex(l-ix+iz,1,1)+dx*(ex(l+ix+iz,1,1)-ex(l-ix+iz,1 &
		,1))

		g=g+dy* &
		(ex(l+iy+iz,1,1)+ex(l-ix+iy+iz,1,1)+dx*(ex(l+ix+iy+iz,1,1) &
		-ex(l-ix+iy+iz,1,1))-g)
		
		ex0=(f+dz*(g-f))*(.25*qm)

        -------------------------------------------------- 
        f=bx(l-iy,1,1)+bx(l-iy-iz,1,1)+dz*(bx(l-iy+iz,1,1)-bx(l-iy-iz,1 &
		,1))
		f=bx(l,1,1)+bx(l-iz,1,1)+dz*(bx(l+iz,1,1)-bx(l-iz,1,1))+f+dy &
		* (bx(l+iy,1,1)+bx(l+iy-iz,1,1)+dz*(bx(l+iy+iz,1,1)-bx(l+iy &
		-iz,1,1))-f)
		g=bx(l+ix-iy,1,1)+bx(l+ix-iy-iz,1,1)+dz*(bx(l+ix-iy+iz,1,1) &
		-bx(l+ix-iy-iz,1,1))
		g=bx(l+ix,1,1)+bx(l+ix-iz,1,1)+dz*(bx(l+ix+iz,1,1)-bx(l+ix-iz,1 &
		,1))+g+dy*(bx(l+ix+iy,1,1)+bx(l+ix+iy-iz,1,1)+dz*(bx(l+ix &
		+iy+iz,1,1)-bx(l+ix+iy-iz,1,1))-g)
		
		bx0=(f+dx*(g-f))*(.125*qm*cinv)
    */

    Matrixd32 trilinear_staggered( double xd, double yd, double zd) 
    {
        /*
        from nodal points
            f(i+dx) = f(i) + dx * (f(i+1)-f(i))
        to staggered grid stored at midpoints
            f at location i+dx  = half of f(i)+f(i-1) + dx*(f(i+1)-f(i-1))
		where now f(i) means f at location i+1/2.

        Then we apply the normal linear volume interpolation
        Note that E and B differ in staggered grid locations
        */

        // using cube class to play in +1/+1/+1 regime
        Matrixd32 c;


        // ex
        c(0,0) = (1 - zd)*((1 - yd)*((0.5 - 0.5*xd)*this->exY(-1,0,0) + 0.5*this->exY(0,0,0) + 0.5*xd*this->exY(1,0,0)) + 
      yd*((0.5 - 0.5*xd)*this->exY(-1,1,0) + 0.5*this->exY(0,1,0) + 0.5*xd*this->exY(1,1,0))) + 
   zd*((1 - yd)*((0.5 - 0.5*xd)*this->exY(-1,0,1) + 0.5*this->exY(0,0,1) + 0.5*xd*this->exY(1,0,1)) + 
      yd*((0.5 - 0.5*xd)*this->exY(-1,1,1) + 0.5*this->exY(0,1,1) + 0.5*xd*this->exY(1,1,1)));

        // ey
        c(1,0) =(1 - zd)*((1 - yd)*(-0.5*(-1. + xd)*(this->eyY(0,-1,0) + this->eyY(0,0,0)) + 0.5*xd*(this->eyY(1,-1,0) + this->eyY(1,0,0))) + 
      yd*(-0.5*(-1. + xd)*(this->eyY(0,0,0) + this->eyY(0,1,0)) + 0.5*xd*(this->eyY(1,0,0) + this->eyY(1,1,0)))) + 
   zd*((1 - yd)*(-0.5*(-1. + xd)*(this->eyY(0,-1,1) + this->eyY(0,0,1)) + 0.5*xd*(this->eyY(1,-1,1) + this->eyY(1,0,1))) + 
      yd*(-0.5*(-1. + xd)*(this->eyY(0,0,1) + this->eyY(0,1,1)) + 0.5*xd*(this->eyY(1,0,1) + this->eyY(1,1,1))));

        // ez
        c(2,0) = (1 - zd)*((1 - yd)*(-0.5*(-1. + xd)*(this->ezY(0,0,-1) + this->ezY(0,0,0)) + 0.5*xd*(this->ezY(1,0,-1) + this->ezY(1,0,0))) + 
      yd*(-0.5*(-1. + xd)*(this->ezY(0,1,-1) + this->ezY(0,1,0)) + 0.5*xd*(this->ezY(1,1,-1) + this->ezY(1,1,0)))) + 
   zd*((1 - yd)*(-0.5*(-1. + xd)*(this->ezY(0,0,0) + this->ezY(0,0,1)) + 0.5*xd*(this->ezY(1,0,0) + this->ezY(1,0,1))) + 
      yd*(-0.5*(-1. + xd)*(this->ezY(0,1,0) + this->ezY(0,1,1)) + 0.5*xd*(this->ezY(1,1,0) + this->ezY(1,1,1))));

        // bx
        c(0,1) =(1 - zd)*((1 - yd)*(-0.25*(-1. + xd)*(this->bxY(0,-1,-1) + this->bxY(0,-1,0) + this->bxY(0,0,-1) + this->bxY(0,0,0)) + 
         0.25*xd*(this->bxY(1,-1,-1) + this->bxY(1,-1,0) + this->bxY(1,0,-1) + this->bxY(1,0,0))) + 
      yd*(-0.25*(-1. + xd)*(this->bxY(0,0,-1) + this->bxY(0,0,0) + this->bxY(0,1,-1) + this->bxY(0,1,0)) + 
         0.25*xd*(this->bxY(1,0,-1) + this->bxY(1,0,0) + this->bxY(1,1,-1) + this->bxY(1,1,0)))) + 
   zd*((1 - yd)*(-0.25*(-1. + xd)*(this->bxY(0,-1,0) + this->bxY(0,-1,1) + this->bxY(0,0,0) + this->bxY(0,0,1)) + 
         0.25*xd*(this->bxY(1,-1,0) + this->bxY(1,-1,1) + this->bxY(1,0,0) + this->bxY(1,0,1))) + 
      yd*(-0.25*(-1. + xd)*(this->bxY(0,0,0) + this->bxY(0,0,1) + this->bxY(0,1,0) + this->bxY(0,1,1)) + 
         0.25*xd*(this->bxY(1,0,0) + this->bxY(1,0,1) + this->bxY(1,1,0) + this->bxY(1,1,1))));

        // by
        c(1,1) = (1 - zd)*((1 - yd)*(-0.25*(-1. + xd)*(this->byY(-1,0,-1) + this->byY(-1,0,0) + this->byY(0,0,-1) + this->byY(0,0,0)) + 
         0.25*xd*(this->byY(0,0,-1) + this->byY(0,0,0) + this->byY(1,0,-1) + this->byY(1,0,0))) + 
      yd*(-0.25*(-1. + xd)*(this->byY(-1,1,-1) + this->byY(-1,1,0) + this->byY(0,1,-1) + this->byY(0,1,0)) + 
         0.25*xd*(this->byY(0,1,-1) + this->byY(0,1,0) + this->byY(1,1,-1) + this->byY(1,1,0)))) + 
   zd*((1 - yd)*(-0.25*(-1. + xd)*(this->byY(-1,0,0) + this->byY(-1,0,1) + this->byY(0,0,0) + this->byY(0,0,1)) + 
         0.25*xd*(this->byY(0,0,0) + this->byY(0,0,1) + this->byY(1,0,0) + this->byY(1,0,1))) + 
      yd*(-0.25*(-1. + xd)*(this->byY(-1,1,0) + this->byY(-1,1,1) + this->byY(0,1,0) + this->byY(0,1,1)) + 
         0.25*xd*(this->byY(0,1,0) + this->byY(0,1,1) + this->byY(1,1,0) + this->byY(1,1,1))));

        // bz
        c(2,1) = (1 - zd)*((1 - yd)*(-0.25*(-1. + xd)*(this->bzY(-1,-1,0) + this->bzY(-1,0,0) + this->bzY(0,-1,0) + this->bzY(0,0,0)) + 
         0.25*xd*(this->bzY(0,-1,0) + this->bzY(0,0,0) + this->bzY(1,-1,0) + this->bzY(1,0,0))) + 
      yd*(-0.25*(-1. + xd)*(this->bzY(-1,0,0) + this->bzY(-1,1,0) + this->bzY(0,0,0) + this->bzY(0,1,0)) + 
         0.25*xd*(this->bzY(0,0,0) + this->bzY(0,1,0) + this->bzY(1,0,0) + this->bzY(1,1,0)))) + 
   zd*((1 - yd)*(-0.25*(-1. + xd)*(this->bzY(-1,-1,1) + this->bzY(-1,0,1) + this->bzY(0,-1,1) + this->bzY(0,0,1)) + 
         0.25*xd*(this->bzY(0,-1,1) + this->bzY(0,0,1) + this->bzY(1,-1,1) + this->bzY(1,0,1))) + 
      yd*(-0.25*(-1. + xd)*(this->bzY(-1,0,1) + this->bzY(-1,1,1) + this->bzY(0,0,1) + this->bzY(0,1,1)) + 
         0.25*xd*(this->bzY(0,0,1) + this->bzY(0,1,1) + this->bzY(1,0,1) + this->bzY(1,1,1))));


        return c;
    };


    Matrixd32 trilinear( double xd, double yd, double zd) 
    {
        // using cube class to play in +1/+1/+1 regime
        Matrixd32 c00, c01, c10, c11, c0, c1, c; 

        c00 = this->F(0,0,0)*(1-xd) + this->F(1,0,0)*xd;
        c01 = this->F(0,0,1)*(1-xd) + this->F(1,0,1)*xd;
        c10 = this->F(0,1,0)*(1-xd) + this->F(1,1,0)*xd;
        c11 = this->F(0,1,1)*(1-xd) + this->F(1,1,1)*xd;

        c0 = c00*(1-yd) + c10*yd;
        c1 = c01*(1-yd) + c11*yd;

        c = c0*(1-zd) + c1*zd;
    
    // cout << "tri lin elements" << endl;
    // cout << "000" << endl;
    // cout << n.F(0,0,0) << endl;
    // cout << "100" << endl;
    // cout << n.F(1,0,0)<< endl;
    // cout << "001" << endl;
    // cout << n.F(0,0,1)<< endl;
    // cout << "101" << endl;
    // cout << n.F(1,0,1)<< endl;
    // cout << "010" << endl;
    // cout << n.F(0,1,0) << endl;
    // cout << "110" << endl;
    // cout << n.F(1,1,0)<< endl;
    // cout << "011" << endl;
    // cout << n.F(0,1,1)<< endl;
    // cout << "111" << endl;
    // cout << n.F(1,1,1)<< endl;

        return c;
    };

};


class Particle_Mover
{
public:

    /*
        Boris/Vay pusher to update particle velocities
    */
	template<
		class CellData,
		class Geometry
	> void update_velocities(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {

        // get local cells
        auto cells = grid.get_cells();

        // loop over cells
        for (const uint64_t& cell: cells) {

            auto* const cell_data = grid[cell];
            if (cell_data->number_of_particles == 0){continue; }

            std::array<double, 3> start_coords = grid.geometry.get_min(cell);
            std::array<double, 3> ds = grid.geometry.get_length(cell);

            // create cell neighborhood; 
            // TODO improve
            EB_Cube n(grid, cell);

            for (auto& particle: cell_data->particles) {
                double 
                    x = particle[0],
                    y = particle[1],
                    z = particle[2],
                    ux = particle[3],
                    uy = particle[4],
                    uz = particle[5];

                // interpolate fields
                double 
                    xd = (x - start_coords[0])/ds[0],
                    yd = (y - start_coords[1])/ds[1],
                    zd = (z - start_coords[2])/ds[2];

                Matrixd32 F = n.trilinear_staggered(xd, yd, zd);

                double 
                    Exi = F(0,0),
                    Eyi = F(1,0),
                    Ezi = F(2,0),
                    Bxi = F(0,1),
                    Byi = F(1,1),
                    Bzi = F(2,1);

                double
                    uxm = ux + q*e*Exi*dt/(2.0*me*c),
                    uym = uy + q*e*Eyi*dt/(2.0*me*c),
                    uzm = uz + q*e*Ezi*dt/(2.0*me*c); 

                // Lorentz transform
                double
                    gamma = sqrt(1.0 + uxm*uxm + uym*uym + uzm*uzm);

                // calculate u'
                double
                   tx = q*e*Bxi*dt/(2.0*gamma*me*c),
                   ty = q*e*Byi*dt/(2.0*gamma*me*c),
                   tz = q*e*Bzi*dt/(2.0*gamma*me*c);

                double
                    ux0 = uxm + uym*tz - uzm*ty,
                    uy0 = uym + uzm*tx - uxm*tz,
                    uz0 = uzm + uxm*ty - uym*tx;

                // calculate u+
                double 
                    sx = 2.0*tx/(1.0 + tx*tx + ty*ty + tz*tz),
                    sy = 2.0*ty/(1.0 + tx*tx + ty*ty + tz*tz),
                    sz = 2.0*tz/(1.0 + tx*tx + ty*ty + tz*tz);

                double
                    uxp = uxm + uy0*sz - uz0*sy,
                    uyp = uym + uz0*sx - ux0*sz,
                    uzp = uzm + ux0*sy - uy0*sx;

                // t-dt/2 -> t + dt/2
                particle[3] = uxp + q*e*Exi*dt/(2.0*me*c);
                particle[4] = uyp + q*e*Eyi*dt/(2.0*me*c);
                particle[5] = uzp + q*e*Ezi*dt/(2.0*me*c);

                // TODO FIXME add particle propagate here
                /*
                double 
                    uxn = particle[3],
                    uyn = particle[4],
                    uzn = particle[5];

                    gamma = sqrt(1.0 + uxn*uxn + uyn*uyn + uzn*uzn);

                    particle[0] += (c*dt/gamma)*uxn;
                    particle[1] += (c*dt/gamma)*uyn;
                    particle[2] += (c*dt/gamma)*uzn;
                */

            }


        } // end of loop over cells


        return;
    };


    /*
        -------------------------------------------------------------------------------- 
        Particle propagator that pushes particles in space

         We basically just compute x_n+1 = x_n + v_x,n * dt/gamma
         and, hence, will advance the (relativistic) Newton-Lorentz 
         equation in time as experienced by the particle in her own frame

        TODO:
        Currently some non-necessary extra work is done because we move
        also particles in ghost cells, i.e., the neighboring cells at process
        boundary.

    */
	template<
		class CellData,
		class Geometry
	> void propagate(
		dccrg::Dccrg<CellData, Geometry>& grid
	) {


        // actually move particles in local cells and copies of remote neighbors
        auto cells = grid.get_cells();
        const std::vector<uint64_t>& remote_neighbors
            = grid.get_remote_cells_on_process_boundary();
        cells.insert(cells.begin(), 
                     remote_neighbors.begin(), 
                     remote_neighbors.end()
                     );

        double gamma = 1.0;

        for (const uint64_t& cell: cells) {

            auto* const cell_data = grid[cell];
            if (cell_data == NULL) {
                cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
                abort();
            }

            for (auto& particle: cell_data->particles) {

                double 
                    ux = particle[3],
                    uy = particle[4],
                    uz = particle[5];

                    gamma = sqrt(1.0 + ux*ux + uy*uy + uz*uz);

                    particle[0] += (c*dt/gamma)*ux;
                    particle[1] += (c*dt/gamma)*uy;
                    particle[2] += (c*dt/gamma)*uz;
            }
        }


        return;
    };






};

#endif
