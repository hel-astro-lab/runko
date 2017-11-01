#pragma once


#include "corgi/corgi.h"


namespace vlasov {


/*! \brief Vlasov cell 
 *
 * Grid infrastructure methods are inherited from corgi::Cell
 */
class VCell : public corgi::Cell {
  public:
    VCell(size_t i, size_t j, int o) : corgi::Cell(i, j, o) { }

    // void bark();

    void bark() { fmt::print("Woof!\n"); };
      

};







} // end of namespace vlasov
