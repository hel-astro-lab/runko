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


  /* check extent (not really necessary on most platforms) */
  MPI_Aint lb, extent, dist[2];
  MPI_Type_get_extent(obj_type, &lb, &extent);
  std::array<pic::Particle,10> particles;

  MPI_Get_address(&particles[0], &dist[0]);
  MPI_Get_address(&particles[1], &dist[1]);

  MPI_Datatype temptype;
  if (extent != (dist[1] - dist[0])) {
    std::cout << "WARNING! Extra padding between arrays in MPI datatype\n";
    temptype = obj_type;
    lb = 0;
    extent = dist[1] - dist[0];

    MPI_Type_create_resized(temptype, lb, extent, &obj_type);
    MPI_Type_commit(&obj_type);
    MPI_Type_free(&temptype);
  }

  return obj_type;
}



} } // ns mpi4cpp::mpi

