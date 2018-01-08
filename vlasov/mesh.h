#pragma once


#include "algorithm"
#include "array"
#include <vector>
#include <unordered_map>




namespace vlasov {



template<typename T, int D>
class Mesh {

  // TODO what should be implemented

  typedef std::array<uint64_t, D> indices_t;
  typedef std::array<T, D> value_t;


  static const uint64_t error_cid = 0;
  static const uint64_t error_index = 0xFFFFFFFFFFFFFFFF;
  int max_refinement_level = 10;
	uint64_t last_cell;
  size_t number_of_blocks = 0;
  
  std::unordered_map<uint64_t, T> data;


  /*
  template<>
  uint64_t index(const std::array<uint64_t, 1>& indices) const {}

  template<>
  uint64_t index(const std::array<uint64_t, 2>& indices) const {}

  template<>
  uint64_t index(const std::array<uint64_t, 3>& indices) const {
		uint64_t cid = 1;

    return cid;
  }

  uint64_t index(const std::array<uint64_t, D>& indices) const {}
  */


  
  indices_t get_indices(uint64_t cid) const 
  {

    if (cid == error_cid || cid > this->last_cid) {
      const indices_t error_indices = {{error_index, error_index, error_index}};
      return error_indices;
    }

    // subtract ids of larger cells
    const int refinement_level = this->get_refinement_level(cid);
    for (int i = 0; i < refinement_level; i++) {
      cid -=
        this->length.get()[0]
        * this->length.get()[1]
        * this->length.get()[2]
        * (uint64_t(1) << (i * 3));
    }

    cid -= 1;	// cell numbering starts from 1
    const indices_t indices = {{

      (cid % (this->length.get()[0] * (uint64_t(1) << refinement_level)))
        * (uint64_t(1) << (max_refinement_level - refinement_level)),

        ((cid / (this->length.get()[0] * (uint64_t(1) << refinement_level)))
         % (this->length.get()[1] * (uint64_t(1) << refinement_level)))
          * (uint64_t(1) << (max_refinement_level - refinement_level)),

        (cid / (
                 this->length.get()[0]
                 * this->length.get()[1]
                 * (uint64_t(1) << (2 * refinement_level))
                ))
          * (uint64_t(1) << (max_refinement_level - refinement_level))
    }};

    return indices;


  }


  uint64_t get_cell(
      const indices_t& indices,
      const int refinement_level
  ) const 
  { 

    if (indices[0] >= this->length.get()[0] * (uint64_t(1) << this->max_refinement_level)) {
      return error_cid;
    }

    if (indices[1] >= this->length.get()[1] * (uint64_t(1) << this->max_refinement_level)) {
      return error_cid;
    }

    if (indices[2] >= this->length.get()[2] * (uint64_t(1) << this->max_refinement_level)) {
      return error_cid;
    }

    if (refinement_level < 0) {
      return error_cid;
    }

    if (refinement_level > this->max_refinement_level) {
      return error_cid;
    }

    // cell numbering starts at 1
    uint64_t cid = 1;

    // add ids of larger cells
    for (int i = 0; i < refinement_level; i++) {
      cid +=
        this->length.get()[0]
        * this->length.get()[1]
        * this->length.get()[2]
        * (uint64_t(1) << (i * 3));
    }

    // convert to indices of this cell's refinement level
    const indices_t this_level_indices = {{
      indices[0] / (uint64_t(1) << (this->max_refinement_level - refinement_level)),
        indices[1] / (uint64_t(1) << (this->max_refinement_level - refinement_level)),
        indices[2] / (uint64_t(1) << (this->max_refinement_level - refinement_level))
    }};

    // get the length of the grid in terms of cells of this refinement level
    const std::array<uint64_t, 2> this_level_length = {{
      this->length.get()[0] * (uint64_t(1) << refinement_level),
        this->length.get()[1] * (uint64_t(1) << refinement_level)
    }};

    cid
      += this_level_indices[0]
      + this_level_indices[1] * this_level_length[0]
      + this_level_indices[2] * this_level_length[0] * this_level_length[1];

    return cid;

  
  
  }
    
  int get_refinement_level(const uint64_t cid) const {}

  uint64_t get_parent(const uint64_t cid) const {}

  uint64_t get_level_0_parent(const uint64_t cid) const {}

  value_t get_level_0_cell_lenght() const {}

	value_t get_length(const uint64_t cid) const {}

	value_t get_center(const uint64_t cid) const {}

	value_t get_min(const uint64_t cid) const {}

	value_t get_max(const uint64_t cid) const {}

	value_t get_center(
		const int refinement_level,
		const indices_t index
	) const {}






};
  












} // end of namespace vlasov
