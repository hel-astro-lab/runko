#pragma once

#include "corgi/corgi.h"
#include "cell.h"
#include "dataContainer.h"


namespace vlasov {

/*! \brief Plasma grid
 *
 * Grid infrastructure methods are inherited from corgi::Node
 */
class Grid : public corgi::Node {
  public:
    std::unordered_map< uint64_t, std::shared_ptr<vlasov::VCell> > cells;

    Grid() : corgi::Node() { }
  
    void howl() { fmt::print("Auuuuuu!\n"); };


    void cycle() {

      // std::unordered_map<uint64_t, std::shared_ptr<corgi::Cell>>::iterator it = cells.begin();
      // while (it != cells.end()) {
      //   vlasov::VCell c = it->second;
      //     
      //   c->data.cycle();
      //   it++;
      // }

      for (auto it: cells) {
        it.second->data.cycle();
      }

    }



};






} // end of namespace vlasov

