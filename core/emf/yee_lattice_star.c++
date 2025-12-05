#include "core/emf/yee_lattice.h"
#include "tyvi/mdgrid.h"
#include "core/emf/tile.h"

#include <tuple>

namespace emf {

void
  YeeLattice::boundary_condition_b_star(const value_type boundary_value)
{
  tyvi::mdgrid_work w {};
  boundary_condition_b_star(w, boundary_value);
  w.wait();
}

void
  YeeLattice::boundary_condition_b_star(const tyvi::mdgrid_work& w, const value_type boundary_value)
//  YeeLattice::boundary_condition_b_star(const tyvi::mdgrid_work& w, Star& star)  //Should add also mins and maxs to the input arguments... Because this-mins[x] does not work in YeeLattice...
{
    const auto Bmds = B_.mds();
    //w.for_each_index(Bmds, [=] __device__ (auto idx) {
    w.for_each_index(Bmds, [=] (auto idx) {

        //const auto ig = idx[0]+this->mins[0];
        //const auto jg = idx[1]+this->mins[1];
        //const auto kg = idx[2]+this->mins[2];

        const auto ig = idx[0]; //+mins[0];
        const auto jg = idx[1]; //+mins[1];
        const auto kg = idx[2]; //+mins[2];


        if (ig < 10 && jg < 10 && kg < 10) {
        
            Bmds[idx][0] = boundary_value;
            Bmds[idx][1] = boundary_value;
            Bmds[idx][2] = boundary_value;
        }
    });
}

void
  YeeLattice::boundary_condition_e_star(const value_type dt)
{
  tyvi::mdgrid_work w {};
  boundary_condition_e_star(w, dt);
  w.wait();
}

void
  YeeLattice::boundary_condition_e_star(const tyvi::mdgrid_work& w, const value_type dt)
{
  return;
}

}  // namespace emf
