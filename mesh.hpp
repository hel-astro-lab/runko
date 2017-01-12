#ifndef MESH_H
#define MESH_H

#include "array"
#include "vector"

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

typedef Array<double, 4, 1> vec4;


// class to handle neighborhood cells
class Neighborhood_Cube
{
public:

    std::array<vec4, 27> fields;

    vec4& operator [](const std::size_t i){return this->fields[i];}
    const vec4& operator [](const std::size_t i) const {return this->fields[i];}

    void init_currents() 
    {
        for(int ijk=0; ijk<27; ijk++) {
            fields[ijk] = vec4::Zero();
        }
        return;
    };


    int index(int i, int j, int k) {
        return (i+1) + (j+1)*3 + (k+1)*3*3;
    };

    void add_J(int i, int j, int k, vec4 dJ)
    {
        int ijk = this->index(i,j,k);

        #ifdef DEBUG
        cout << "  addJ: " << i << j << k << " -> " << ijk 
        << " A" << dJ << " B " << fields[ijk] << endl;
        #endif

        fields[ijk] += dJ;
    };

    vec4 J(int i, int j, int k) {
        return this->fields[index(i,j,k)];
    };


	void distribute(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid, uint64_t cell);
};



class Mesh
{
public:

    // MPI parameters
    int 
        rank = 0;

    // Deposit currents into the mesh
	void sort_particles_into_cells(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);

    // Deposit currents into the mesh
	void deposit_currents(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);

    // Collect and compute staggered Yee currents
	void yee_currents(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);

    /* Shift fields from staggered Yee to nodal lattice
       in practice we interpolate linearly between neighboring Yee values
       TODO FIXME
    */
	void nodal_fields(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);

};



#endif
