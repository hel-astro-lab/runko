#pragma once

#include <vector>

namespace toolbox {

/// \brief Container for storing multiple time steps of the simulation
template <class T>
class Rotator {
  std::vector<T> container;

  public:

  size_t currentStep = 0;

  void push_back(T vm);

  T* getPtr();
  T& getRef();


  T* getNewPtr();
  T& getNewRef();

  T* getAll(size_t cs);

  // FIXME raw cycling for time step index
  void cycle();

};


} // end of namespace toolbox


/// method to add data into the container
template <class T>
void toolbox::Rotator<T>::push_back(T vm) {
  container.push_back(vm);
}


/// Get current element
template <class T>
T* toolbox::Rotator<T>::getPtr() {
  // fmt::print("getting from Rotator with {}\n", currentStep);
  return (T*) &(container[ currentStep ]);
}

/// Return reference to current element
template <class T>
T& toolbox::Rotator<T>::getRef() {
  return container[ currentStep ];
}


/// get a fresh container that we can update into
template <class T>
T* toolbox::Rotator<T>::getNewPtr() {
  if (currentStep == 0) return (T*) &(container[1]);
  if (currentStep == 1) return (T*) &(container[0]);

  return NULL;
}

/// Get reference to a new element
template <class T>
T& toolbox::Rotator<T>::getNewRef() {
  if (currentStep == 0) return container[1];
  if (currentStep == 1) return container[0];
}



/// Get any arbitrary snapshot from the container
template <class T>
T* toolbox::Rotator<T>::getAll(size_t cs) {
  // fmt::print("pulling from Rotator with {}\n", cs);
  return (T*) &(container[cs]);
}


/// raw cycling for time step index
// NOTE: could be done better
template <class T>
void toolbox::Rotator<T>::cycle() {
  // fmt::print(" calling cycle (originally {})\n", currentStep);
  currentStep++;

  // check bounds and cycle back
  if (currentStep > 1) currentStep = 0;
}



