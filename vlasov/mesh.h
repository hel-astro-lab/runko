#pragma once


#include <algorithm>
#include <array>
#include <vector>
#include <unordered_map>




namespace toolbox {



template<typename T, int D>
class AdaptiveMesh {

  public:

  typedef std::array<uint64_t, D> indices_t;
  typedef std::array<T, D> value_array_t;
  using iterator       = typename std::unordered_map<uint64_t, T>::iterator;
  using const_iterator = typename std::unordered_map<uint64_t, T>::const_iterator;


  static const uint64_t error_cid = 0;
  static const uint64_t error_index = 0xFFFFFFFFFFFFFFFF;
  int max_refinement_level = 10;
	uint64_t last_cid;
  size_t number_of_blocks = 0;
  
  std::unordered_map<uint64_t, T> data;

  std::array<uint64_t, D> length;


  /*
  AdaptiveMesh(indices_t length) : length(length) 
  {
    static_assert(sizeof(length) == D, "number of given array lengths must match dimension");

  }
  */


  AdaptiveMesh() {}

  void resize(indices_t& given_length)
  {
    length = given_length;
  }

  void set(uint64_t key, T val)
  {
    data[key] = val;
  }

  T get(uint64_t key) const 
  {
    const_iterator it = data.find(key);
    return it == data.end() ? T(0) : it->second;
  }


  // set item
  T operator() (indices_t indx, int refinement_level) 
  {
    uint64_t cid = get_cell(indx, refinement_level);
    return data[cid];
  }


  // get item
  const T operator() (indices_t indx, int refinement_level) const
  {
    uint64_t cid = get_cell(indx, refinement_level);
    const_iterator it = data.find(cid);
    return it == data.end() ? T(0) : it->second;
  }



	uint64_t get_last_cid() const
	{
		return last_cid;
	}


	void update_last_cid()
	{
		const uint64_t grid_length
			= length[0] * length[1] * length[2];

		last_cid = 0;
		for (int i = 0; i <= max_refinement_level; i++) {
			last_cid += grid_length * (uint64_t(1) << (i * 3));
		}
	}
  

  indices_t get_indices(uint64_t cid) const 
  {

    if (cid == error_cid || cid > last_cid) {
      const indices_t error_indices = {{error_index, error_index, error_index}};
      return error_indices;
    }

    // subtract ids of larger cells
    const int refinement_level = get_refinement_level(cid);
    for (int i = 0; i < refinement_level; i++) {
      cid -=
        length[0]
        * length[1]
        * length[2]
        * (uint64_t(1) << (i * 3));
    }

    cid -= 1;	// cell numbering starts from 1
    const indices_t indices = {{

      (cid % (length[0] * (uint64_t(1) << refinement_level)))
        * (uint64_t(1) << (max_refinement_level - refinement_level)),

        ((cid / (length[0] * (uint64_t(1) << refinement_level)))
         % (length[1] * (uint64_t(1) << refinement_level)))
          * (uint64_t(1) << (max_refinement_level - refinement_level)),

        (cid / (
                 length[0]
                 * length[1]
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

    if (indices[0] >= length[0] * (uint64_t(1) << max_refinement_level)) {
      return error_cid;
    }

    if (indices[1] >= length[1] * (uint64_t(1) << max_refinement_level)) {
      return error_cid;
    }

    if (indices[2] >= length[2] * (uint64_t(1) << max_refinement_level)) {
      return error_cid;
    }

    if (refinement_level < 0) {
      return error_cid;
    }

    if (refinement_level > max_refinement_level) {
      return error_cid;
    }

    // cell numbering starts at 1
    uint64_t cid = 1;

    // add ids of larger cells
    for (int i = 0; i < refinement_level; i++) {
      cid +=
          length[0]
        * length[1]
        * length[2]
        * (uint64_t(1) << (i * 3));
    }

    // convert to indices of this cell's refinement level
    const indices_t this_level_indices = 
    {{
        indices[0] / (uint64_t(1) << (max_refinement_level - refinement_level)),
        indices[1] / (uint64_t(1) << (max_refinement_level - refinement_level)),
        indices[2] / (uint64_t(1) << (max_refinement_level - refinement_level))
    }};

    // get the length of the grid in terms of cells of this refinement level
    const std::array<uint64_t, 2> this_level_length = 
    {{
        length[0] * (uint64_t(1) << refinement_level),
        length[1] * (uint64_t(1) << refinement_level)
    }};

    cid
      += this_level_indices[0]
      + this_level_indices[1] * this_level_length[0]
      + this_level_indices[2] * this_level_length[0] * this_level_length[1];

    return cid;

  
  
  }
    

  int get_refinement_level(const uint64_t cid) const 
  {
		if (cid == error_cid || cid > last_cid) {
			return -1;
		}

		int refinement_level = 0;
		uint64_t current_last = 0;

		while (refinement_level <= max_refinement_level) {
			current_last +=
				  length[0]
				* length[1]
				* length[2]
				* (uint64_t(1) << 3 * refinement_level);

			if (cid <= current_last) {
				break;
			}

			refinement_level++;
		}

		if (refinement_level > max_refinement_level) {
			return -1;
		}

		return refinement_level;
	}


  int get_maximum_refinement_level() 
	{
		const uint64_t grid_length
			= length[0] * length[1] * length[2];
		int refinement_level = 0;
		double current_last = 0;
		while (current_last <= double(~uint64_t(0))) {
			// TODO: don't assume 3 dimensions
        
			current_last += double(grid_length) * pow(double(8), double(refinement_level));
			refinement_level++;
		}

		return refinement_level - 2;
	}


  bool set_maximum_refinement_level(const int given_refinement_level)
  {
    if(given_refinement_level < 0) return false;

    max_refinement_level = given_refinement_level;

    return true;
  }



  uint64_t get_parent(const uint64_t cid) const
	{
		const int refinement_level = get_refinement_level(cid);

		if (refinement_level < 0
		|| refinement_level > max_refinement_level) {
			return error_cid;
		}

		if (refinement_level == 0) {
			return cid;
		}

		return get_cell_from_indices(get_indices(cid), refinement_level - 1);
	}


  uint64_t get_level_0_parent(const uint64_t cid) const
	{
		const int refinement_level = get_refinement_level(cid);

		if (refinement_level < 0
		|| refinement_level > max_refinement_level) {
			return error_cid;
		}

		if (refinement_level == 0) {
			return cid;
		}

		return get_cell_from_indices(get_indices(cid), 0);
	}


  std::vector<uint64_t> et_all_children(const uint64_t cid) const 
	{
		std::vector<uint64_t> children;

		if (cid == error_cid) {
			return children;
		}

		if (this->cell_process.count(cid) == 0) {
			return children;
		}

		// given cell cannot have children
		int refinement_level = get_refinement_level(cid);
		if (refinement_level >= get_maximum_refinement_level()) {
			return children;
		}

		children.reserve(8);

		indices_t indices = get_indices(cid);

		// get indices of next refinement level within this cell
		refinement_level++;
		const uint64_t index_offset
			= (uint64_t(1) << (get_maximum_refinement_level() - refinement_level));

		for (uint64_t
			z_index_offset = 0;
			z_index_offset < 2 * index_offset;
			z_index_offset += index_offset
		)
		for (uint64_t
			y_index_offset = 0;
			y_index_offset < 2 * index_offset;
			y_index_offset += index_offset
		)
		for (uint64_t
			x_index_offset = 0;
			x_index_offset < 2 * index_offset;
			x_index_offset += index_offset
		) {
			const indices_t index = {{
				indices[0] + x_index_offset,
				indices[1] + y_index_offset,
				indices[2] + z_index_offset,
			}};

			children.push_back(
				get_cell_from_indices(index, refinement_level)
			);
		}

		return children;
	}




  /*
  value_array_t get_level_0_cell_lenght() const {}

	value_array_t get_length(const uint64_t cid) const {}

	value_array_t get_center(const uint64_t cid) const {}

	value_array_t get_min(const uint64_t cid) const {}

	value_array_t get_max(const uint64_t cid) const {}

	value_array_t get_center(
		const int refinement_level,
		const indices_t index
	) const {}

  */



};






} // end of namespace vlasov
