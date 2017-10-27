#ifndef PRTCS_H
#define PRTCS_H


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wmissing-braces"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#pragma clang diagnostic pop

#include "cell.hpp"
#include "parameters.hpp"
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
	EB_Cube(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid, uint64_t cell);


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

    Matrixd32 trilinear_staggered( double xd, double yd, double zd);


    Matrixd32 trilinear( double xd, double yd, double zd);

};


class Particle_Mover
{
public:

    /*
        Boris/Vay pusher to update particle velocities
    */
	void update_velocities(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


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
	void propagate(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


};

#endif
