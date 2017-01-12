#include "inject.hpp"

double maxwellian(double vb)
{
    RNDM<double> rng;

    double fmax = 0.5*(1.0 + exp(-2.0*vb*vb));
    double vmin = 0.0;
    double vmax = 5.0*vb;
    double vf = rng.dunif(vmin, vmax);

    double f = exp(-(vb-vf)*(vb-vf)/2);
    double x = fmax*rng.dunif(0.0, 1.0);

    if(x > f){
        return maxwellian(vb);
    } 

    return vf;
}


double double_maxwellian(double vb)
{
    RNDM<double> rng;

    double fmax = 0.5*(1.0 + exp(-2.0*vb*vb));
    double vmin = -5.0*vb;
    double vmax =  5.0*vb;
    double vf = rng.dunif(vmin, vmax);

    double f = 0.5*(exp(-(vf-vb)*(vf-vb)/2) +
            exp(-(vf+vb)*(vf+vb)/2));
    double x = fmax*rng.dunif(0.0, 1.0);

    if(x > f){
        return double_maxwellian(vb);
    } 

    return vf;
}

        
Injector::Injector(double given_xmin,
        double given_xmax,
        double given_ymin,
        double given_ymax,
        double given_zmin,
        double given_zmax
) {

        this->grid_xmin = given_xmin;
        this->grid_xmax = given_xmax;

        this->grid_ymin = given_ymin;
        this->grid_ymax = given_ymax;

        this->grid_zmin = given_zmin;
        this->grid_zmax = given_zmax;

}

    // cylinder shape injection
void Injector::cylinder(
        dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid, 
        uint64_t Np, 
        double vb
) {

    RNDM<double> rng;

    // XXX: not needed in cylinder
    // double grid_xmid = this->grid_xmin + 0.5*(this->grid_xmax - this->grid_xmin);
    double grid_ymid = this->grid_ymin + 0.5*(this->grid_ymax - this->grid_ymin);
    double grid_zmid = this->grid_zmin + 0.5*(this->grid_zmax - this->grid_zmin);


    uint64_t p=0;
    while(p < Np){
        // cout << rank << ": inserting particle " << p << endl;

        // cylinder inside inj_r radius from xmin to xmax
        double 
            inj_x = rng.dunif(this->grid_xmin, this->grid_xmax),
            inj_r = rng.dunif(0.0, 0.05),
            inj_a = rng.dunif(0.0, 6.28318530718); // [0,2pi[

        double 
            inj_y = grid_ymid + inj_r*cos(inj_a),
                  inj_z = grid_zmid + inj_r*sin(inj_a);

        std::array<double, 3> coords = {{inj_x, inj_y, inj_z}};

        uint64_t cell = grid.geometry.get_cell(0, coords);

        if(grid.is_local(cell)) {

            // cout << rank << ": is local!" << endl;
            auto* const cell_data = grid[cell];

            std::array<double, 6> loc_vel = {{
                inj_x,
                inj_y,
                inj_z,
                maxwellian(vb),
                0.0,
                0.0
            }};

            // insert particle to the cell
            cell_data->particles.push_back(loc_vel);
            cell_data->number_of_particles = cell_data->particles.size();

        } 

        p++;
    };


    return;
}


