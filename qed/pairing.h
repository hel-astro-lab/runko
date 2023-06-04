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
  using std::min;
  using std::max;



// duplicate particle info into fresh variables
inline auto duplicate_prtcl(
    string t1, float_p ux1, float_p uy1, float_p uz1, float_p w1
    ) -> std::tuple<string, float_p, float_p, float_p, float_p>
{
  return {t1, ux1, uy1, uz1, w1};
}



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

  // normalization factor for probabilities
  float_p prob_norm = 1.0f;

  //--------------------------------------------------
  //--------------------------------------------------
  //--------------------------------------------------
    
  // random numbers between [0, 1[
  float_p rand() { return uni_dis(gen); };
  
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

  //--------------------------------------------------
  // check if interaction list is empty for type t1
  bool is_empty(string& t1) 
  {
    int i=0;
    for(auto iptr : interactions){
      if(t1 == iptr->t1) i += 1;
    }
    return (i > 0) ? false : true; 
  }


  //--------------------------------------------------
  template<size_t D>
  void solve(pic::Tile<D>& tile)
  {

    // TODO
    for(auto&& con1 : tile.containers)
    {
      con1.sort_in_rev_energy();
      //con1.update_cumulative_arrays();
    }

    // variables inside loop
    float_p ux1, uy1, uz1, w1, ux2, uy2, uz2, w2;
    float_p ux3, uy3, uz3, w3, ux4, uy4, uz4, w4;
    float_p cross_s, wmin, wmax, prob;
    float_p p_ini, p_tar;


    //--------------------------------------------------
    // loop over incident types
    for(auto&& con1 : tile.containers)
    {
      auto t1 = con1.type;
      if(is_empty(t1)) continue; // no interactions with incident type t1

      std::cout << "container type:" << t1 << std::endl;

      // loop over target types
      for(auto&& con2 : tile.containers)
      {
        auto t2 = con2.type;

        // loop over interactions
        for(auto iptr : interactions){

          // do type matching of containers and interactions
          if( (t1 == iptr->t1) && (t2 == iptr->t2) ){


            // loop over incident particles
            //UniIter::iterate([=] DEVCALLABLE (
            //          size_t n, 
            //          pic::ParticleContainer<D>& con
            //          ){
            for(size_t n1=0; n1<con1.size(); n1++) {

              //unpack incident 
              ux1 = con1.vel(0,n1);
              uy1 = con1.vel(1,n1);
              uz1 = con1.vel(2,n1);
              w1  = con1.wgt(n1);

              // loop over targets
              for(size_t n2=0; n2<con2.size(); n2++) {

                // unpack target
                ux2 = con2.vel(0,n2);
                uy2 = con2.vel(1,n2);
                uz2 = con2.vel(2,n2);
                w2  = con2.wgt(  n2);

                // interaction cross section
                cross_s = iptr->comp_cross_section(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2 );
                
                // interaction probability
                wmin = min(w1, w2);
                wmax = max(w1, w2);
                prob = cross_s*w1*w2/prob_norm;
                // NOTE: difference of all2all scheme is here where prob depends on w1*w2

                //-------------------------------------------------- 
                if(rand() < prob){

                  // particle values after interaction
                  auto [t3, ux3, uy3, uz3, w3] = duplicate_prtcl(t1, ux1, uy1, uz1, w1);
                  auto [t4, ux4, uy4, uz4, w4] = duplicate_prtcl(t2, ux2, uy2, uz2, w2);

                  // interact and udpate variables in-place
                  iptr->interact( t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4);

                  p_ini = w2/wmax;
                  p_tar = w1/wmax;

                  //std::cout << " interacting:" << prob << " pini/tar" << p_ini << " " << p_tar << std::endl;

                  //-------------------------------------------------- 
                  if(rand() < p_ini){
                    if(t1 == t3){ // if type is conserved only update the prtcl info
                      con1.vel(0,n1) = ux3;
                      con1.vel(1,n1) = uy3;
                      con1.vel(2,n1) = uz3;
                    } else { // else destroy previous and add new 
                      // TODO 

                    }
                  }
                  //-------------------------------------------------- 

                  //-------------------------------------------------- 
                  if(rand() < p_tar){
                    if(t2 == t4){ // if type is conserved only update the prtcl info
                      con2.vel(0,n2) = ux4;
                      con2.vel(1,n2) = uy4;
                      con2.vel(2,n2) = uz4;
                    } else { // else destroy previous and add new 
                      // TODO 

                    }
                  }
                  //-------------------------------------------------- 

                } // end of if prob
                //-------------------------------------------------- 
              } // end of loop over con2 particles
            } // end of loop over con1 particles
            //}, con.size(), con);
            //UniIter::sync();
          }
        }
      }
    }

    return;
  }



};





} // end of namespace qed
