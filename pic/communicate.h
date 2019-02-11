#pragma once

#include <mpi4cpp/mpi.h>

// introduce class to mpi
namespace mpi4cpp { namespace mpi {



template <>
struct is_mpi_datatype<::pic::Particle> 
  : mpl::true_ { };


template <>
inline MPI_Datatype get_mpi_datatype<::pic::Particle>(
    const ::pic::Particle& obj) 
{
    
  // get start of the class
  MPI_Aint base;
  MPI_Get_address( &obj, &base ); 

  // how many elements per each type
  std::array<int, 3> block_lengths{
    { 7, 1, 1} 
  };

  // and then the actual members
  std::array<MPI_Aint, 3> member_offsets; // relative offsets
  MPI_Get_address( &obj.data[0], &member_offsets[0]);
  MPI_Get_address( &obj._id,     &member_offsets[1]);
  MPI_Get_address( &obj._proc,   &member_offsets[2]);

  // create real (absolute) offsets (=rel - base)
  std::array<MPI_Aint, 3> offsets {
    member_offsets[0] - base,
    member_offsets[1] - base,
    member_offsets[2] - base,
  };

  // introduce datatypes
  std::array<MPI_Datatype, 3> datatypes{
    {MPI_DOUBLE, MPI_INT, MPI_INT}
  };

  //--------------------------------------------------
  // create datatype; this is standard format and should not be changed
  MPI_Datatype obj_type;
  MPI_Type_create_struct(
      block_lengths.size(),
      block_lengths.data(),
      offsets.data(),
      datatypes.data(),
      &obj_type);

  MPI_Type_commit(&obj_type);
  return obj_type;
}



} } // ns mpi4cpp::mpi

