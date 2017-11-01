#pragma once

#include "corgi/corgi.h"



namespace vlasov {

/*! \brief Plasma grid
 *
 * Grid infrastructure methods are inherited from corgi::Node
 */
class Grid : public corgi::Node {
  public:
    Grid() : corgi::Node() { }
  
    void howl() { fmt::print("Auuuuuu!\n"); };





};






} // end of namespace vlasov

