#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <unordered_map>


#include "definitions.h"
#include "bundles.h"
#include "sheets.h"



namespace vmesh {


  static const uint64_t error_block = 0;
  static const uint64_t error_index = 0xFFFFFFFFFFFFFFFF;

  /** \brief Sparse block-based velocity mesh implementation
   *
   *
   *
   */
  class VeloMesh {

    public:
      size_t number_of_blocks = 0;

      /// main hashmap block container
      std::unordered_map<uint64_t, vblock_t> blockContainer;

      /// returns a pointer to the data of given block
      vblock_t* operator [] (const uint64_t cellID) const
      {
        if (this->blockContainer.count(cellID) > 0) {
          return (vblock_t*) &(this->blockContainer.at(cellID));
        } else {
          return NULL;
        }
      }


      std::array<double, 3> mins, maxs, lens; // Geometry parameters
      indices_t Nblocks = {{ NBLOCKS, NBLOCKS, NBLOCKS }};

      /// Clipping threshold
      Real threshold = 1.0e-4;


      void zFill( std::array<double, 3> mins_,
          std::array<double, 3> maxs_);

      vblock_t getBlock( const uint64_t cellID ) const;

      uint64_t getBlockID( const indices_t index ) const;

      indices_t getIndices( uint64_t cellID );

      std::array<double, 3> getSize( const uint64_t cellID );

      std::array<double, 3> getCenter( const uint64_t cellID );
      std::array<double, 3> getCenterIndx( const indices_t indx );

      std::vector<double> getXGrid();
      std::vector<double> getYGrid();
      std::vector<double> getZGrid();

      std::vector<uint64_t> allBlocks( bool sorted = false);

      bundles::Bundle getBundle(size_t, size_t, size_t);

      void addBundle(size_t, size_t, size_t, bundles::Bundle);

      sheets::Sheet getSheet(size_t, size_t);
      void addSheet(size_t, size_t, sheets::Sheet);
      void setSheet(size_t, size_t, sheets::Sheet);

      bool clip( );

      size_t sizeInBytes() const;
      size_t capacityInBytes() const;





      // --------------------------------------------------
      /* For python interface we need to define __getitem__ and __setitem__
       * TODO: implement with lambda operator instead in pyplasmatools.c++
      */
      vblock_t __getitem__(const uint64_t cellID) const {
        return this->blockContainer.at( cellID );
      };

      vblock_t __getitem2__(const size_t i, 
                            const size_t j, 
                            const size_t k) const {
        uint64_t cellID = this->getBlockID( {{i,j,k}} );
        // fmt::print("({},{},{}) = {}\n",i,j,k,cellID);
        return this->__getitem__(cellID);
      };

      void __setitem__(const uint64_t cellID, const vblock_t vals) {
        blockContainer[cellID] = vals;
      }

      void __setitem2__(const size_t i, 
                        const size_t j, 
                        const size_t k, 
                        vblock_t vals) {
        uint64_t cellID = this->getBlockID( {{i,j,k}} );
        blockContainer[cellID] = vals;
      };

  }; // end of VeloMesh class header

} // end of vmesh namespace





