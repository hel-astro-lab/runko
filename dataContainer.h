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
  T& getRef();


  T* getNew();
  T& getNewRef();

  T* getAll(size_t cs);

  // FIXME raw cycling for time step index
  void cycle();

};


} // end of namespace datarotators


/// method to add data into the container
template <class T>
void datarotators::DataContainer<T>::push_back(T vm) {
  container.push_back(vm);
}


/// Get current element
template <class T>
T* datarotators::DataContainer<T>::get() {
  // fmt::print("getting from DataContainer with {}\n", currentStep);
  return (T*) &(container[ currentStep ]);
}

/// Return reference to current element
template <class T>
T& datarotators::DataContainer<T>::getRef() {
  return &(container[ currentStep ]);
}


/// get a fresh container that we can update into
template <class T>
T* datarotators::DataContainer<T>::getNew() {
  if (currentStep == 0) return (T*) &(container[1]);
  if (currentStep == 1) return (T*) &(container[0]);

  return NULL;
}

/// Get reference to a new element
template <class T>
T& datarotators::DataContainer<T>::getNewRef() {
  if (currentStep == 0) return &(container[1]);
  if (currentStep == 1) return &(container[0]);
}



/// Get any arbitrary snapshot from the container
template <class T>
T* datarotators::DataContainer<T>::getAll(size_t cs) {
  // fmt::print("pulling from DataContainer with {}\n", cs);
  return (T*) &(container[cs]);
}


/// raw cycling for time step index
// NOTE: could be done better
template <class T>
void datarotators::DataContainer<T>::cycle() {
  // fmt::print(" calling cycle (originally {})\n", currentStep);
  currentStep++;

  // check bounds and cycle back
  if (currentStep > 1) currentStep = 0;
}



