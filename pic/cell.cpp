
#include "cell.hpp"


std::tuple<void*, int, MPI_Datatype> Cell::get_mpi_datatype()
{

    void* address = NULL;
    int count = -1;
    MPI_Datatype datatype = MPI_DATATYPE_NULL;

    switch(Cell::transfer_mode) {
        case Cell::INIT:
            break;
        case Cell::ELECTRONS:
            if (this->electrons.size() > 0) {
                address = &(this->electrons[0]);
            } else {
                // return a sane address just in case
                address = &(this->number_of(Population::ELECTRONS));
            }

            count = this->electrons.size() * 6;
            datatype = MPI_DOUBLE;
            break;
        case Cell::POSITRONS:
            if (this->positrons.size() > 0) {
                address = &(this->positrons[0]);
            } else {
                // return a sane address just in case
                address = &(this->number_of(Population::POSITRONS));
            }

            count = this->positrons.size() * 6;
            datatype = MPI_DOUBLE;
            break;
        case Cell::NUMBER_OF_ELECTRONS:
            address = &(this->number_of(Population::ELECTRONS));
            count = 1;
            datatype = MPI_UNSIGNED;
            break;
        case Cell::NUMBER_OF_POSITRONS:
            address = &(this->number_of(Population::POSITRONS));
            count = 1;
            datatype = MPI_UNSIGNED;
            break;
        case Cell::REMOTE_NEIGHBOR_LIST:
            address = &(this->remote_neighbor_list[0]);
            count = 27;
            datatype = MPI_UINT64_T;
            break;
        case Cell::INCOMING_CURRENTS:
            // address = this->incoming_currents.data();
            address = &(this->incoming_currents[0]);
            count = 108;
            datatype = MPI_DOUBLE;
            break;
        case Cell::TYPE:
            address = &(this->cell_type);
            count = 1;
            datatype = MPI_DOUBLE;
            break;
        case Cell::CURRENT:
            address = this->J.data();
            count = 4;
            datatype = MPI_DOUBLE;
            break;
        case Cell::YEE_CURRENT:
            address = this->JY.data();
            // FIXME count = 4;
            count = 3;
            datatype = MPI_DOUBLE;
            break;
        case Cell::YEE_B:
            address = this->BY.data();
            count = 3;
            datatype = MPI_DOUBLE;
            break;
        case Cell::YEE_E:
            address = this->EY.data();
            count = 3;
            datatype = MPI_DOUBLE;
            break;
        default:
            std::cerr << __FILE__ << ":" << __LINE__
                << " Invalid transfer switch: " << Cell::transfer_mode
                << std::endl;
            abort();
            break;
    }

    return std::make_tuple(address, count, datatype);
}
