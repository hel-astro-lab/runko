#pragma once

#include <vector>
#include <iostream>


namespace toolbox {

/// \brief Container for storing multiple time steps of the simulation
template <class T, size_t L>
class Rotator {
  std::vector<T> container;

  public:

  size_t current_step = 0;


  /// method to add data into the container
  void push_back(T vm) {container.push_back(vm); };

  /// general index
  inline size_t index(size_t i) const {return (i + current_step) % L ; };

  /// get reference to the current time step
  inline T& get(size_t i=0) { return container[ index(i) ]; };

  /// get reference to the current time step
  inline const T& get(size_t i=0) const { return container[ index(i) ]; };

  // FIXME raw cycling for time step index
  void cycle() {
    current_step++;

    // check bounds and cycle back
    if (current_step >= L) current_step = 0;

    //std::cout << "ROTATOR STEP IS " << current_step << '\n';
  }


};


} // end of namespace toolbox


