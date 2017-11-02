#pragma once

#include <vector>

namespace datarotators {

/// \brief Container for storing multiple time steps of the simulation
template <class T>
class DataContainer {
  std::vector<T> container;

  public:

  size_t currentStep = 0;

  void push_back(T vm);

  T* get();

  T* getNew();

  T* getAll(size_t cs);

  // FIXME raw cycling for time step index
  void cycle();

};


} // end of namespace datarotators
