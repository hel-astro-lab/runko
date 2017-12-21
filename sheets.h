#pragma once

#include <vector>
#include <array>


#include "definitions.h"



namespace sheets {
  // -------------------------------------------------- 
  /// Sheet of velocities from vmesh
  // These are always slices of the full mesh along some dimension at some 
  // location.
  class Sheet {

    public:

      /// guiding grid along horizontal dimensions of the sheet (i.e., x)
      std::vector<Realf> iGrid;

      // sheet size in horizontal dim
      size_t Ni = 0;

      /// guiding grid along vertical dimensions of the sheet (i.e., y)
      std::vector<Realf> jGrid;

      // sheet size in horizontal dim
      size_t Nj = 0;

      /// coordinate of the sheet in the slicing dimension
      double sliceVal;

      /// Value storage of the sheet
      std::vector<Realf> values;


      void resize(size_t Ni_, size_t Nj_);

      size_t getIndex(size_t i, size_t j);

      void loadValue(size_t i, size_t j, Realf val);

      void loadZeroBlock(size_t i, size_t j);

      void loadBlock(size_t i, size_t j, vblock_t block);

      vblock_t getBlock(size_t i, size_t j);

      bool isNonZero(size_t i, size_t j);

      void checkSizes(const Sheet& s);

      Realf sum();

      // Sheet diff();


      // Sheet arithmetics
      Sheet& operator+=(const Sheet& rhs);

      Sheet& operator-=(const Sheet& rhs);

      Sheet& operator*=(const Sheet& rhs);

      Sheet& operator*=(const Realf rhs);
  };

  inline Sheet operator+(Sheet lhs, const Sheet& rhs) {
    lhs += rhs;
    return lhs;
  };

  inline Sheet operator-(Sheet lhs, const Sheet& rhs) {
    lhs -= rhs;
    return lhs;
  };

  inline Sheet operator*(Sheet lhs, const Sheet& rhs) {
    lhs *= rhs;
    return lhs;
  };

  inline Sheet operator*(Sheet lhs, const Realf rhs) {
    lhs *= rhs;
    return lhs;
  };

}
