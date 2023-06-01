#pragma once

#include <string>
#include <tuple>
#include <random>
#include "../../definitions.h"
#include "../../pic/tile.h"


namespace qed {
  using std::string;
  using std::tuple;



// Monte-Carlo pairing of particles
class Pairing
{
private:
  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float_p> uni_dis;


public:

  // constructor with incident/target types
  Pairing() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)       
  { }


  template<size_t D>
  void solve(pic::Tile<D>& tile)
  {
    for(auto&& container : tile.containers)
    {
      std::cout << "container type:" << container.type << std::endl;
    }
  }



};





} // end of namespace qed
