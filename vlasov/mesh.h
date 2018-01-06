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


  static const uint64_t error_block = 0;
  static const uint64_t error_index = 0xFFFFFFFFFFFFFFFF;
  int max_refinement_level;
	uint64_t last_cell;
  size_t number_of_blocks = 0;
  
  std::unordered_map<uint64_t, T> data;


  template<>
  uint64_t index(const std::array<uint64_t, 1>& indices) const {}

  template<>
  uint64_t index(const std::array<uint64_t, 2>& indices) const {}

  template<>
  uint64_t index(const std::array<uint64_t, 3>& indices) const {
		uint64_t cid = 1;

    return cid
  }

  uint64_t index(const std::array<uint64_t, D>& indices) const {}



  
  indices_t get_indices(uint64_t cid) const {}

  uint64_t get_cell(
      const indices_t& indices,
      const int refinement_level
  ) const { }
    
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
