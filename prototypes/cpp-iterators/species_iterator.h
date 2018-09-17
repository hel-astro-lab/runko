// iterator to cycle through species

#include <vector>



template<typename M, typename T>
class PairPlasmaIterator {

  private:

    /// First species 
    size_t _beginning = 0;

    /// beyond last species 
    size_t _ending    = 2;

    /// internal pointer to the parent object/container (Mother object)
    M* ptr;

    /// Iterators own internal counting system
    size_t spcs = 0; 


  public:

    PairPlasmaIterator(M& rhs) : ptr(&rhs) {}
    PairPlasmaIterator(const M& rhs) : ptr(&rhs) {}

    PairPlasmaIterator(const PairPlasmaIterator& rhs) : ptr(rhs.ptr) {}


    /// Assignment
    PairPlasmaIterator& operator= (const PairPlasmaIterator& rhs) = default;

    /// iterate
    PairPlasmaIterator& operator++ () {
      ++this->spcs;
      return *this;
    }

    /// Referencing tile interiors
    T& operator *() {
      if(spcs == 0) return (T&) (ptr->electrons);
      else if(spcs == 1) return (T&) (ptr->positrons);
      else throw std::range_error("iterator goes beyond electrons (0) or positrons (1)");
    }


    /// return charge to mass ratio of the current species under consideration
    T getChargeToMass() {
      return ptr->getQ(spcs);
    }
      

    /// iterate with steps
    // PairPlasmaIterator operator++ (int) {
    //   PairPlasmaIterator temp(*ptr);
    //   ++*this;
    //   return (temp);
    // }

    /// equal comparison done by comparing internal spcs value
    bool operator== (PairPlasmaIterator& rhs) const {
      return (ptr == rhs.ptr) && (spcs == rhs.spcs);
    }

    /// unequal comparison done by comparing internal spcs value
    bool operator!= (PairPlasmaIterator& rhs) const {
      return (ptr != rhs.ptr) || (spcs != rhs.spcs);
    }

    /// Returns an iterator pointing to the first element in the sequence
    PairPlasmaIterator begin() {
      PairPlasmaIterator temp(*ptr);
      temp.spcs = _beginning;
      return temp;
    }

    /// Returns an iterator pointing to the past-the-end element in the sequence
    PairPlasmaIterator end() {
      PairPlasmaIterator temp(*ptr);
      temp.spcs = _ending;
      return temp ;
    }


};


