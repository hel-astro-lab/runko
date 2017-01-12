#ifndef FIELD_H
#define FIELD_H


#include "vector"
#include "array"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wmissing-braces"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#pragma clang diagnostic pop

#include "cell.hpp"
#include "parameters.hpp"
#include "common.h"

#include <Eigen/Dense>

using namespace std;
using namespace dccrg;
using namespace Eigen;


typedef Array<double, 4, 1> vec4;
// typedef Array<double, 3, 1> vec3;
typedef Array<double, 4, 4> mat4;


// class to handle neighborhood cells
class Field_Cube
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


class Field_Solver
{
public:

    // Mesh() {};

    // MPI parameters
    int 
        rank = 0;

    // Advance half of the B field
    void push_half_B(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


    // Advance E field
	void push_E(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);

};


#endif
