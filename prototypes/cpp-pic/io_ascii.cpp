#include "io.hpp"

void Save::init()
{

    if(this->rank == 0) { 
        this->visit_particles.open(this->save_dir + this->filename +"_prtcls.visit");
        this->visit_grid.open(this->save_dir + this->filename +"_grid.visit");
        this->visit_fields.open(this->save_dir + this->filename +"_fields.visit");


        this->visit_particles << "!NBLOCKS " << this->comm_size << "\n";
        this->visit_grid      << "!NBLOCKS " << this->comm_size << "\n";
        this->visit_fields    << "!NBLOCKS " << this->comm_size << "\n";
    }

    this->init_flag = true;

    return;
}

void Save::finalize()
{

    // close streams
    if (this->rank == 0) {
        this->visit_particles.close();
        this->visit_grid.close();
        this->visit_fields.close();
    }

    return;
}

void Save::save_grid(
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
    ) {

    // write the grid
    const string grid_file_name(
            this->save_dir + this->filename + "_"
            + std::to_string(this->rank) + "_" 
            + std::to_string(this->step) + "_grid.vtk"
            ); 
    grid.write_vtk_file(grid_file_name.c_str());

    return;
}

void Save::save_particles(
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    // write the grid
    const string outname(
            this->save_dir + this->filename + "_"
            + std::to_string(this->rank) + "_" 
            + std::to_string(this->step) + "_prtcls.vtk"
            ); 

    std::ofstream outfile;
    outfile.open(outname.c_str());

    const vector<uint64_t> cells = grid.get_cells();


    // write header
    outfile << "# vtk DataFile Version 2.0\n"
        "Particle test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    // calculate the total number of particles
    uint64_t total_particles = 0;
    for (const auto& cell: cells) {
        auto* const cell_data = grid[cell];
        if (cell_data == NULL) {
            cerr << __FILE__ << ":" << __LINE__ << ": no data for cell " << cell << endl;
            abort();
        }

        total_particles += cell_data->number_of(Population::ELECTRONS);
    }
    outfile << "POINTS " << total_particles << " float\n";

    // write out coordinates
    uint64_t written_particles = 0;
    for (const auto& cell: cells) {

        auto* const cell_data = grid[cell];
        for (vector<std::array<double, 6> >::const_iterator
                coordinates = cell_data->particles(Population::ELECTRONS).begin();
                coordinates != cell_data->particles(Population::ELECTRONS).end();
                coordinates++, written_particles++
            ) {
            outfile << (*coordinates)[0] << "\t"
                << (*coordinates)[1] << "\t"
                << (*coordinates)[2] << "\n";
            // << (*coordinates)[3] << "\t"
            // << (*coordinates)[4] << "\t"
            // << (*coordinates)[5] << "\n";
        }
    }
    if (written_particles != total_particles) {
        cerr << __FILE__ << ":" << __LINE__
            << ": Incorrect number of particles written: " << written_particles
            << ", should be " << total_particles
            << endl;
        abort();
    }

    outfile << "CELLS " << total_particles << " " << 2 * total_particles << "\n";
    for (uint64_t i = 0; i < total_particles; i++) {
        outfile << "1 " << i << "\n";
    }

    outfile << "CELL_TYPES " << total_particles << "\n";
    for (uint64_t i = 0; i < total_particles; i++) {
        outfile << "1 ";
    }

    // process numbers of particles
    outfile << "\nCELL_DATA " << total_particles
        << "\nSCALARS process int 1\nLOOKUP_TABLE default\n";
    for (uint64_t i = 0; i < total_particles; i++) {
        outfile << rank << " ";
    }
    outfile << "\n";

    // cell numbers of particles
    outfile << "SCALARS cell int 1\nLOOKUP_TABLE default\n";
    for (const auto& cell: cells) {
        auto* const cell_data = grid[cell];
        for (vector<std::array<double, 6> >::const_iterator
                coordinates = cell_data->particles(Population::ELECTRONS).begin();
                coordinates != cell_data->particles(Population::ELECTRONS).end();
                coordinates++, written_particles++
            ) {
            outfile << cell << " ";
        }
    }
    outfile << "\n";

    outfile.close();

    return;
}

void Save::save_fields(
	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>& grid
) {

    // write the grid
    const string outname(
            this->save_dir + this->filename + "_"
            + std::to_string(this->rank) + "_" 
            + std::to_string(this->step) + "_fields.vtk"
            ); 

    std::ofstream outfile;
    outfile.open(outname.c_str());
    if (!outfile.is_open()) {
        std::cerr << "Couldn't open file " << outname << std::endl;
        exit(EXIT_FAILURE);
    }


    outfile << "# vtk DataFile Version 2.0" << std::endl;
    outfile << "GUMICS data" << std::endl;
    outfile << "ASCII" << std::endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // write separate points for every cells' corners
    const vector<uint64_t> cells = grid.get_cells();

    outfile << "POINTS " << cells.size() * 8 << " float" << std::endl;
    for (const auto& cell: cells) {
        // auto* const cell_data = grid[cell];

        // get x,y,z coordinates of cell
        const std::array<double, 3>
            cell_min = grid.geometry.get_min(cell),
                     cell_max = grid.geometry.get_max(cell);

        outfile
            << cell_min[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
            << cell_max[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
            << cell_min[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
            << cell_max[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
            << cell_min[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
            << cell_max[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
            << cell_min[0] << " " << cell_max[1] << " " << cell_max[2] << "\n"
            << cell_max[0] << " " << cell_max[1] << " " << cell_max[2] << "\n";
    }

    // map cells to written points
    outfile << "CELLS " << cells.size() << " " << cells.size() * 9 << std::endl;
    for (unsigned int j = 0; j < cells.size(); j++) {
        outfile << "8 ";
        for (int i = 0; i < 8; i++) {
            outfile << j * 8 + i << " ";
        }
        outfile << std::endl;
    }

    // cell types
    outfile << "CELL_TYPES " << cells.size() << std::endl;
    for (unsigned int i = 0; i < cells.size(); i++) {
        outfile << 11 << std::endl;
    }

    /*
       Write simulation data
       */
    outfile << "CELL_DATA " << cells.size() << endl;

    // cells' process
    outfile << "SCALARS process int 1" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
        // outfile << grid.cell_process.at(*cell) << "\n";
        outfile << rank << "\n";
    }

    // cells' refinement level
    /*
       outfile << "SCALARS refinement_level int 1" << endl;
       outfile << "LOOKUP_TABLE default" << endl;
       for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
       outfile << mapping.get_refinement_level(*cell) << "\n";
       }
       */


    // density rho
    outfile << "SCALARS rho float 1" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (const auto& cell: cells) {
        auto* const cell_data = grid[cell];
        outfile << cell_data->J(0) << "\n";
    }


    // current vector
    outfile << "VECTORS current float" << endl;
    for (const auto& cell: cells) {
        auto* const cell_data = grid[cell];
        outfile << cell_data->J(1) << " ";
        outfile << cell_data->J(2) << " ";
        outfile << cell_data->J(3) << " ";

        /* TODO J or JY
           outfile << cell_data->JY(0) << " ";
           outfile << cell_data->JY(1) << " ";
           outfile << cell_data->JY(2) << " ";
           */
        outfile << "\n";
    }


    outfile << "VECTORS E float" << endl;
    for (const auto& cell: cells) {
        auto* const cell_data = grid[cell];
        outfile << cell_data->EY(0) << " ";
        outfile << cell_data->EY(1) << " ";
        outfile << cell_data->EY(2) << " ";
        outfile << "\n";
    }

    outfile << "VECTORS B float" << endl;
    for (const auto& cell: cells) {
        auto* const cell_data = grid[cell];
        outfile << cell_data->BY(0) << " ";
        outfile << cell_data->BY(1) << " ";
        outfile << cell_data->BY(2) << " ";
        outfile << "\n";
    }


    outfile << "\n";
    outfile.close();

    return;
}

void Save::update_master_list()
{
    // if it does not exists create the list
    if (!this->init_flag) {
        this->init();
    }

    // now actually append lists (if rank0)
    if(this->rank == 0) {
        for (int i = 0; i < this->comm_size; i++) {
            this->visit_particles << this->filename << "_"
                << i << "_"
                << this->step << "_prtcls.vtk\n";
            this->visit_grid << this->filename << "_"
                << i << "_"
                << this->step << "_grid.vtk\n";
            this->visit_fields << this->filename << "_"
                << i << "_"
                << this->step << "_fields.vtk\n";
        }
    }

    return;
}

