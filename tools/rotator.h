#pragma once

#include <vector>
#include <iostream>


namespace toolbox {

/// \brief Container for storing multiple time steps of the simulation
template <class T, size_t L>
class Rotator {
  std::vector<T> container;

  public:

  size_t currentStep = 0;


  /// method to add data into the container
  void push_back(T vm) {container.push_back(vm); };

  /// general index
  inline size_t index(size_t i) {return (i + currentStep) % L ; };

  /// get reference to the current time step
  inline T& get(size_t i=0) { return container[ index(i) ]; };


  // FIXME raw cycling for time step index
  void cycle() {
    currentStep++;

    // check bounds and cycle back
    if (currentStep >= L) currentStep = 0;

    //std::cout << "ROTATOR STEP IS " << currentStep << '\n';
  }


};


} // end of namespace toolbox


