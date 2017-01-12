#ifndef IO_H
#define IO_H


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


class Save
{
public:

    // MPI parameters
    int 
        rank = 0,
        comm_size = 1;

    // current simulation file step
    unsigned int
        step = 0;

    string 
        save_dir,
        filename;

    // checker if master file is created
    bool init_flag = false;


    std::ofstream 
        visit_particles, 
        visit_grid, 
        visit_fields;
    
    // default constructor
    Save(){
        this->rank = 0;
        this->save_dir = "out/";
        this->filename = "pic";

        return;
    }

    // Initialize files and streams for writing
    void init();

    void finalize();


    /*
        Write grid to file using internal dccrg functions
    */
	void save_grid(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


    /*
        Write *ALL* particles to file:
        NOTE this is very IO heavy
    */
	void save_particles(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


    /*
        Write field/mesh quantities
        NOTE: this should be the choice for visualization 
    */
	void save_fields(dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid);


    /*
        Add current file names to a master list
    */
    void update_master_list();


};

#endif
