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

/// Draw y location from cosh distribution using rejection sampling
double cosh_distribution(
        double delta,
        double ymin,
        double ymax,
        double Ly
){
    RNDM<double> rng;

    double fmax = 1.0;
    double y = rng.dunif(ymin, ymax);

    double f = 0.0;
    if(y < Ly/2)
    {
        f = pow(cosh((y - Ly/4.0)/delta),-2);
    } else {
        f = pow(cosh((y - 3.0*Ly/4.0)/delta),-2);
    }

    double x = rng.dunif(0.0, fmax);

    if(x > f){
        return cosh_distribution(delta, ymin, ymax, Ly);
    } 

    return y;
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
            cell_data->electrons.push_back(loc_vel);
            cell_data->number_of(Population::ELECTRONS) = cell_data->electrons.size();



        } 

        p++;
    };


    return;
}


void Injector::two_sheets(
        dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid, 
        uint64_t Np, 
        double vb,
        double delta
) {

    RNDM<double> rng;

    double ymin = this->grid_ymin;
    double ymax = this->grid_ymax;
    double Ly = ymax - ymin;


    uint64_t p=0;
    while(p < Np){
        // cout << rank << ": inserting particle " << p << endl;

        double 
            inj_x = rng.dunif(this->grid_xmin, this->grid_xmax),
            inj_y = cosh_distribution(delta, ymin, ymax, Ly),
            inj_z = rng.dunif(this->grid_zmin, this->grid_zmax);


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
                maxwellian(vb)
            }};

            // insert particle to the cell
            cell_data->electrons.push_back(loc_vel);
            cell_data->number_of(Population::ELECTRONS) = cell_data->electrons.size();

        } 

        p++;
    };


    return;
}



void Injector::uniform(
        dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid, 
        uint64_t Np, 
        double vb
) {

    RNDM<double> rng;

    // Uniform electron-positron plasma injection, i.e. we inject equal amount of both
    for (int ipops=Population::ELECTRONS; ipops <= Population::POSITRONS; ipops++)
    {
        Population ptype = (Population) ipops; ///< type cast into Population type for later usage

             
        uint64_t p=0;
        while(p < Np){
            // cout << rank << ": inserting particle " << p << endl;

            double 
                inj_x = rng.dunif(this->grid_xmin, this->grid_xmax),
                inj_y = rng.dunif(this->grid_ymin, this->grid_ymax),
                inj_z = rng.dunif(this->grid_zmin, this->grid_zmax);


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
                        0.5*maxwellian(vb),
                        maxwellian(vb)
                }};

                // insert particle to the cell
                cell_data->particles(ptype).push_back(loc_vel);
                cell_data->number_of(ptype) = cell_data->particles(ptype).size();
            } 

            p++;
        }
    };


    return;
}
