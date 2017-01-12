#ifndef INJECT_H
#define INJECT_H

#include "array"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wmissing-braces"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#pragma clang diagnostic pop

#include "cell.hpp"
#include "rng.hpp"


using namespace std;
using namespace dccrg;


/* 
    Auxiliary functions
*/

double maxwellian(double vb);

double double_maxwellian(double vb);


/* Main injector */
class Injector
{
public:
   
    // grid limits
    double grid_xmin,
           grid_xmax,
           grid_ymin,
           grid_ymax,
           grid_zmin,
           grid_zmax;


    Injector(double given_xmin,
             double given_xmax,
             double given_ymin,
             double given_ymax,
             double given_zmin,
             double given_zmax
             );


    // cylinder shape injection
	void cylinder( dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid, uint64_t Np, double vb);

};


#endif
