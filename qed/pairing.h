#pragma once

#include <string>
#include <tuple>
#include <random>
#include <memory>

#include "interactions/interaction.h"
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

  using InteractionPtr = std::shared_ptr<Interaction>;

public:

  // constructor with incident/target types
  Pairing() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)       
  { }

  //using Tile_map = std::unordered_map<TileID_t, Tileptr>;
  //Tile_map tiles; /// Map with tile_id & tile data

  std::vector<InteractionPtr> interactions;


  // add interactions to internal memory of the class; 
  // done via pointers to handle pybind interface w/ python
  void add_interaction(InteractionPtr iptr)
  {
    assert(iptr); // check that we are not appending nullptr

    auto name = iptr->name;
    auto t1 = iptr->t1;
    auto t2 = iptr->t2;

    std::cout << " adding: " << name << " of t1/t2 " << t1 << " " << t2 << std::endl;
    interactions.push_back(iptr);
  }

  // check if interaction list is empty for type t1
  bool is_empty(string& t1) 
  {
    int i=0;
    for(auto iptr : interactions){
      if(t1 == iptr->t1) i += 1;
    }

    return (i > 0) ? false : true; 
  }

  template<size_t D>
  void solve(pic::Tile<D>& tile)
  {
    for(auto&& con1 : tile.containers)
    {
      auto t1 = con1.type;
      if(is_empty(t1)) continue; // no interactions with incident type t1

      std::cout << "container type:" << t1 << std::endl;
    }
  }



};





} // end of namespace qed
