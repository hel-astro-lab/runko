#pragma once

#include <type_traits>
#include <typeinfo>
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
  //--------------------------------------------------
  // check if we have already created the object once
  MPI_Datatype obj_type;
  //obj_type = mpi4cpp::mpi::detail::mpi_datatype_cache().datatype<::pic::Particle>();
  std::type_info const* t = &typeid(::pic::Particle);
  obj_type = mpi4cpp::mpi::detail::mpi_datatype_cache().get(t);

  if(obj_type != MPI_DATATYPE_NULL) return obj_type;


  //--------------------------------------------------
  // if not; we'll create it now
    
  // get start of the class
  MPI_Aint base;
  MPI_Get_address( &obj, &base ); 

  // how many elements per each type
  std::array<int, 9> block_lengths{ {1,1,1,1,1,1,1,1,1} };

  // and then the actual members
  std::array<MPI_Aint, 9> member_offsets; // relative offsets
  MPI_Get_address( &obj.x,    &member_offsets[0]);
  MPI_Get_address( &obj.y,    &member_offsets[1]);
  MPI_Get_address( &obj.z,    &member_offsets[2]);
  MPI_Get_address( &obj.ux,   &member_offsets[3]);
  MPI_Get_address( &obj.uy,   &member_offsets[4]);
  MPI_Get_address( &obj.uz,   &member_offsets[5]);
  MPI_Get_address( &obj.w,    &member_offsets[6]);
  MPI_Get_address( &obj.id,   &member_offsets[7]);
  MPI_Get_address( &obj.proc, &member_offsets[8]);

  // create real (absolute) offsets (=rel - base)
  std::array<MPI_Aint, 9> offsets {
    member_offsets[0] - base,
    member_offsets[1] - base,
    member_offsets[2] - base,
    member_offsets[3] - base,
    member_offsets[4] - base,
    member_offsets[5] - base,
    member_offsets[6] - base,
    member_offsets[7] - base,
    member_offsets[8] - base,
  };


  // introduce datatypes; 
  std::array<MPI_Datatype, 9> datatypes{
    {
      MPI_FLOAT_TP, MPI_FLOAT_TP, MPI_FLOAT_TP, 
      MPI_FLOAT_TP, MPI_FLOAT_TP, MPI_FLOAT_TP, 
      MPI_FLOAT_TP, 
      MPI_INT, MPI_INT}
  };

  //--------------------------------------------------
  // create datatype; this is standard format and should not be changed
  MPI_Type_create_struct(
      block_lengths.size(),
      block_lengths.data(),
      offsets.data(),
      datatypes.data(),
      &obj_type);

  MPI_Type_commit(&obj_type);


  // check and correct for extent in arrays
  MPI_Aint lb, extent, dist[2];
  MPI_Type_get_extent(obj_type, &lb, &extent);
  std::array<pic::Particle,10> particles;

  MPI_Get_address(&particles[0], &dist[0]);
  MPI_Get_address(&particles[1], &dist[1]);

  MPI_Datatype temptype;
  if (extent != (dist[1] - dist[0])) {
    //std::cout << "WARNING! Extra padding between arrays in MPI datatype\n";
    temptype = obj_type;
    lb = 0;
    extent = dist[1] - dist[0];

    MPI_Type_create_resized(temptype, lb, extent, &obj_type);
    MPI_Type_commit(&obj_type);
    MPI_Type_free(&temptype);
  }
  
  // save datatype to internal cache
  mpi4cpp::mpi::detail::mpi_datatype_cache().set(t, obj_type);

  return obj_type;
}



} } // ns mpi4cpp::mpi

