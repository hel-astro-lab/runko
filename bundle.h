
#include <array>
#include <vector>

#include "definitions.h"


namespace bundles {
  /* \brief Bundle of pencils
   *
   *
   */
  class Bundle {

    /// values along the pencil
    std::vector<Realf> pencil;

    /// guiding grid for the pencil
    std::vector<Realf> grid;

    public:

    void resize( size_t N);

    size_t size();

    void loadZeroBlock(size_t q);

    void loadBlock(size_t q, vblock_t block);

    void loadGrid(size_t q, Realf val);

    std::vector<Realf> getGrid();

    std::vector<Realf> getPencil();

    bool isNonZero(size_t q);

    vblock_t getSlice(size_t q);

    Realf getDx(size_t q);

  }; // end of Bundle class

} // end of bundles namespace
