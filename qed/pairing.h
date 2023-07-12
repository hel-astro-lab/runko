#pragma once

#include <string>
#include <tuple>
#include <random>
#include <memory>
#include <map>
#include <functional>

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



// Binary all2all pairing of particles
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

    // using raw pointer instead of smart ptrs; it does not take ownership of the object
    // so the container is not deleted when the temporary storage goes out of scope.
    using ConPtr = pic::ParticleContainer<D>* ;
    //using ConPtr = std::weak_ptr< pic::ParticleContainer<D> >; // this could also work
    //using ConPtr = std::reference_wrapper< pic::ParticleContainer<D> >; // this maybe as well

      
    // build pointer map of types to containers; used as a helper to access particle tyeps
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    // test cons storage calling
    //std::cout << "calling cons with e- and returns:" << cons["e-"]->type << " " << cons["e-"]->size() << std::endl;
    //std::cout << "calling cons with e+ and returns:" << cons["e+"]->type << " " << cons["e+"]->size() << std::endl;
    //std::cout << "calling cons with ph and returns:" << cons["ph"]->type << " " << cons["ph"]->size() << std::endl;


    //--------------------------------------------------
    // call pre-iteration functions to update internal arrays 
    // TODO
    for(auto&& con : tile.containers)
    {
      con.sort_in_rev_energy();
      //con.update_cumulative_arrays();
      con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
    }


    //--------------------------------------------------
    // variables inside loop
    float_p lx1, ly1, lz1,     lx2, ly2, lz2;
    float_p ux1, uy1, uz1, w1, ux2, uy2, uz2, w2;
    float_p ux3, uy3, uz3, w3, ux4, uy4, uz4, w4;
    float_p cross_s, wmin, wmax, prob;
    float_p p_ini, p_tar;

    // loop over interactions
    for(auto iptr : interactions){

      //--------------------------------------------------
      // loop over incident types
      for(auto&& con1 : tile.containers)
      {
        auto t1 = con1.type;
        if(is_empty(t1)) continue; // no interactions with incident type t1

        //std::cout << "container type:" << t1 << std::endl;

        // loop over target types
        for(auto&& con2 : tile.containers)
        {
          auto t2 = con2.type;

          // do type matching of containers and interactions
          if( (t1 == iptr->t1) && (t2 == iptr->t2) ){

            // loop over incident particles
            //UniIter::iterate([=] DEVCALLABLE (
            //          size_t n, 
            //          pic::ParticleContainer<D>& con
            //          ){
            for(size_t n1=0; n1<con1.size(); n1++) {

              // loop over targets
              for(size_t n2=0; n2<con2.size(); n2++) {

                // NOTE: incident needs to be unpacked in the innermost loop, since 
                // some interactions modify its value during the iteration
                  
                //unpack incident 
                lx1 = con1.loc(0,n1);
                ly1 = con1.loc(1,n1);
                lz1 = con1.loc(2,n1);
                  
                ux1 = con1.vel(0,n1);
                uy1 = con1.vel(1,n1);
                uz1 = con1.vel(2,n1);
                w1  = con1.wgt(n1);

                if(w1 < EPS) continue; // omit zero-w incidents

                // unpack target
                lx2 = con2.loc(0,n2);
                ly2 = con2.loc(1,n2);
                lz2 = con2.loc(2,n2);

                ux2 = con2.vel(0,n2);
                uy2 = con2.vel(1,n2);
                uz2 = con2.vel(2,n2);
                w2  = con2.wgt(  n2);

                if(w2 < EPS) continue; // omit zero-w targets

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
                                    
                      // NOTE: we keep location the same
                      con1.vel(0,n1) = ux3;
                      con1.vel(1,n1) = uy3;
                      con1.vel(2,n1) = uz3;
                    } else { // else destroy previous and add new 

                      // destroy current
                      con1.to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU version
                      con1.wgt(n1) = 0.0f; // make zero wgt so its omitted from loop

                      // add new
                      // TODO add_particle
                      cons[t3]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w1);

                      //std::cout << "killing t1" << t1 << std::endl;
                      //std::cout << "adding t3" << t3 << std::endl;
                    }
                  }
                  //-------------------------------------------------- 

                  //-------------------------------------------------- 
                  if(rand() < p_tar){
                    if(t2 == t4){ // if type is conserved only update the prtcl info

                      // NOTE: we keep location the same
                      con2.vel(0,n2) = ux4;
                      con2.vel(1,n2) = uy4;
                      con2.vel(2,n2) = uz4;
                    } else { // else destroy previous and add new 

                      // destroy current
                      con2.to_other_tiles.push_back( {0,0,0,n2} ); // NOTE: CPU version
                      con2.wgt(n2) = 0.0f; // make zero wgt so its omitted from loop

                      // TODO add_particle
                      cons[t4]->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w2);

                      //std::cout << "killing t2" << t2 << std::endl;
                      //std::cout << "adding t4" << t4 << std::endl;
                    }
                  }
                  //-------------------------------------------------- 

                } // end of if prob
                //-------------------------------------------------- 
              } // end of loop over con2 particles
            } // end of loop over con1 particles
            //}, con.size(), con);
            //UniIter::sync();
          } // con types match interaction
        }//con2
      } // con1
    }//end of loop over types


    //--------------------------------------------------
    // final book keeping routines

    for(auto&& con : tile.containers)
    {
      con.delete_transferred_particles(); // remove annihilated prtcls; this transfer storage 
                                          // is used as a tmp container for storing the indices
    }


    return;
  }



  //--------------------------------------------------
  template<size_t D>
  void solve_mc(pic::Tile<D>& tile)
  {

    // using raw pointer instead of smart ptrs; it does not take ownership of the object
    // so the container is not deleted when the temporary storage goes out of scope.
    using ConPtr = pic::ParticleContainer<D>* ;

    // build pointer map of types to containers; used as a helper to access particle tyeps
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //--------------------------------------------------
    // call pre-iteration functions to update internal arrays 
    // TODO
    for(auto&& con : tile.containers)
    {
      con.sort_in_rev_energy();
      //con.update_cumulative_arrays();
      con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
    }





};





} // end of namespace qed
