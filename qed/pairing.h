#pragma once

#include <algorithm>
#include <string>
#include <tuple>
#include <random>
#include <memory>
#include <map>
#include <functional>

#include "interactions/interaction.h"
#include "../../definitions.h"
#include "../../pic/tile.h"
#include "../../tools/sample_arrays.h"



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

  using InteractionPtr = std::shared_ptr<qed::Interaction>;
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
  // auxiliary containers for Monte Carlo sampling 
  std::vector<double> 
      probs,    // maximum probabillity
      wsums,    // sum over possible target weights
      cmaxs;    // maximum cross section
                 
  std::vector<int> 
      jmins,    // minimum indices of prtcls that can particiapte in initeraction
      jmaxs;    // maximum -||-
                 
  std::vector<size_t> 
      ids;     // internal id of the interaction in the storage

  std::map<std::string, double> info_max_int_cs;
                 
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

    info_max_int_cs[name] = 0.0;
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


  // compute maximum partial interaction rates for each process 
  // that LP of type t1 and energy of e1 can experience.
  template<size_t D>
  void comp_pmax(string t1, float_p e1, 
                 std::map<std::string, pic::ParticleContainer<D>*>& cons)
  {

    //size_t n_ints = interactions.size(); // get number of interactions

    probs.clear(); // maximum probabillity
    wsums.clear(); // sum over possible target weights
    cmaxs.clear(); // maximum cross section
    jmins.clear(); // minimum indices of prtcls that can particiapte in initeraction
    jmaxs.clear(); // maximum -||-
    ids.clear(); // internal id of the interaction in the storage
                 
    size_t id = 0;
    for(auto iptr : interactions){
      if(t1 == iptr->t1)
      {
        auto t2 = iptr->t2;      // get target
        auto con_tar = cons[t2]; // target container
        size_t N2 = con_tar->size(); // total number of particles

        double cross_max = iptr->cross_section; // maximum cross section (including x2 for head-on collisions)

        // NOTE: assumes that target distribution remains static for the duration of the time step.
        // In that case, no LP changes energy and the limits finding is ok.
        // This assumption is valid in the limit of small time step dt << mean time between interactions

        //#find the min/max interaction energies and the corresponding array indices jmin/jmax
        auto [emin, emax] = iptr->get_minmax_ene(t1, t2, e1);

        //--------------------------------------------------
        // double counting prevention since we only consider targets with energy more than incident particle

        if(false){

          // require e1 < e2 = e_tar
          // e_tar = [emin, emax]
          // therefore, if e1 > emin we need to set emin' = max(e1, emin)
          // also, if e1 > emax -> interaction is not possible
            
          // NOTE: not needed 
          //if(e1 > emax) {
          //  id++; // remembering to increase counter nevertheless
          //  continue; // no possible targets to interact with
          //}
            
          //if(emin < e1) emin = e1; 
          emin = std::max(e1, emin);
          
          // other way around; require e1 > e2 = etarget
          // e_tar = [emin, emax]
          //
          // therefore emax' = min(e1, emax)
          // and if emin > e1 --> not possible
          //if(emin > e1){
          //  id++; // remembering to increase counter nevertheless
          //  continue; // no possible targets to interact with
          //}
          //emax = std::min(e1, emax);
        }
        //--------------------------------------------------


        // NOTE: assuming reverse order here since ene is sorted in decreasing order
        // NOTE: assumes that sort_in_rev_energy is called in the container; this creates eneArr
        int jmin = toolbox::find_rev_sorted_nearest( con_tar->eneArr, emax );
        int jmax = toolbox::find_rev_sorted_nearest( con_tar->eneArr, emin );
        int jme  = toolbox::find_rev_sorted_nearest( con_tar->eneArr, e1 );

        // TEST remove energy check
        //jmin = 0;
        //jmax = N2;
                
        // total weight of prtcls between jmin-jmax
        // ver1: calculates total weight for every particle; value evolves dynamically
        double wsum = 0.0;
        for(size_t j=jmin; j<jmax; j++) wsum += con_tar->wgt( j ); // sum( w[jmin:jmax] )

        // FIXME: switch to this version as it is faster
        // ver2: assumes static targets; value calculated in the beginnig of loop
        double wsum2 = 0.0;
        if(jmin <= jmax) { // in the opposite case arrays dont span a range and so wsum2 = 0
          double wsum_min = jmin == 0 ? 0.0 : con_tar->wgtCumArr[jmin];
          wsum2 = con_tar->wgtCumArr[jmax-1] - wsum_min;
        }

        // FIXME: using ver2 now; is it correct?

        //std::cout << "wsum " <<  wsum << " " << wsum2 << " jminmax " << jmin << " " << jmax << std::endl;
        //assert( std::abs( wsum - wsum2 ) < EPS );


        //--------------------------------------------------
        //total sum of target weights
        double wtot = 0.0;
        for(size_t j=0; j<N2; j++) wtot += con_tar->wgt( j ); // sum( w )
        if(wtot <= 0.0) wtot = 1.0; // guard for NaNs

        //--------------------------------------------------
        // maximum partial interaction rate
        double par_int_rate = 2.0*cross_max * wsum2; // factor 2 comes from v_rel = 2c

        // update/accept interaction only if it has non-zero prob to occur
        if(par_int_rate > 0.0) {

          // add stuff 
          probs.push_back(par_int_rate);
          cmaxs.push_back(cross_max);
          wsums.push_back(wsum2);

          jmins.push_back(jmin);
          jmaxs.push_back(jmax);

          ids.push_back(id); // update interaction numbering
                           
          //std::cout << "comp_pmax: " << id << " " << iptr->name << " " << t1 << "/" << t2;
          //std::cout << " cmaxs:" << cross_max;
          //std::cout << " probs:" << par_int_rate << " prob: " << par_int_rate/prob_norm;
          //std::cout << " wsum:" << wsum << " wsum/wtot: " << wsum/wtot;
          //std::cout << " wsum2:" << wsum2;
          //std::cout << " js:" << jmin << " " << jme << " " << jmax;
          //std::cout << " emin:" << emin << " e1" << e1 << " emax:" << emax;
          //std::cout << std::endl;
        }
                           
      }

      id++; // increase counter
    }

    return;
  }

  int draw_rand_proc()
  {
    size_t N = probs.size();
    //double ptot = toolbox::sum(probs); // TODO
    float_p ptot = 0.0;
    for(size_t i=0; i<N; i++) ptot += probs[i];

    // calculate cumulative sum of all interactions
    ManVec<float_p> cumsum_int_probs;
    cumsum_int_probs.resize(N);

    cumsum_int_probs[0] = probs[0];
    for(size_t i=1; i<N; i++) cumsum_int_probs[i] = cumsum_int_probs[i-1] + probs[i];

    // normalize
    for(size_t i=0; i<N; i++) cumsum_int_probs[i] = cumsum_int_probs[i]/ptot;

    // draw random interaction from cumulative distribution
    int reti = toolbox::sample_prob(cumsum_int_probs, rand() );

    //std::cout << "draw_rand_proc: ";
    //std::cout << " ptot:" << ptot;
    //std::cout << " reti:" << reti;

    //std::cout << " probs:";
    //for(size_t i=0; i<N; i++) std::cout << probs[i] << " , ";

    //std::cout << " cumsum:";
    //for(size_t i=0; i<N; i++) std::cout << cumsum_int_probs[i] << " , ";

    //std::cout << std::endl;

    // get corresponding id
    return reti;
  }


  // energy-dependent weight adapation functions
  float_p ene_weight_funs(std::string t, float_p x) 
  {
    if(       t == "ph") { return 1.0/x; //std::pow(x, -0.5); 
    } else if(t == "e-") { return 1.0; //std::pow(x, +0.2);
    } else if(t == "e+") { return 1.0; //std::pow(x, +0.2);
    }

    assert(false);
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
    float_p e1, e2;
    float_p wmin, wmax, prob;
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
            //for(int n1=con1.size()-1; n1>=0; n1--) {
            for(size_t n1=0; n1<con1.size(); n1++) { // FIXME

              // loop over targets
              //for(int n2=con2.size()-1; n2>=0; n2--) {
              for(size_t n2=0; n2<con2.size(); n2++) { // FIXME

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

                e1  = con1.eneArr[n1]; 

                if(w1 < EPS) continue; // omit zero-w incidents

                // unpack target
                lx2 = con2.loc(0,n2);
                ly2 = con2.loc(1,n2);
                lz2 = con2.loc(2,n2);

                ux2 = con2.vel(0,n2);
                uy2 = con2.vel(1,n2);
                uz2 = con2.vel(2,n2);
                w2  = con2.wgt(  n2);

                e2  = con2.eneArr[n2]; 

                if(w2 < EPS) continue; // omit zero-w targets

                //--------------------------------------------------
                //avoid double counting by considering only e1 < e2 cases
                //if(e1 > e2) continue;

                // interaction cross section
                auto [cm, vrel] = iptr->comp_cross_section(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2 );
                
                // interaction probability
                wmin = min(w1, w2);
                wmax = max(w1, w2);
                prob = cm*vrel*w1*w2/prob_norm;
                // NOTE: difference of all2all scheme is here where prob depends on w1*w2

                // exponential waiting time between interactions
                double t_free = -log( rand() )/prob;

                //-------------------------------------------------- 
                if(t_free < 1.0){

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
    for(auto&& con : tile.containers) {
      con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
    }

    // keep this ordering; initialization of arrays assumes this way of calling the functions
    // NOTE: cannot move this inside the loop because particle removal assumes that indices remain static
    for(auto&& con : tile.containers) {
      con.sort_in_rev_energy();
      con.update_cumulative_arrays();
    }

    // FIXME: is it better to keep the cumulative wgt array in Particles container or create/destroy it here?


    //--------------------------------------------------
    // collect statistics for bookkeeping
    std::map<std::string, int> info_prtcl_num;
    for(auto&& con : tile.containers) {
      auto t1 = con.type;
      info_prtcl_num[t1] = con.size();
    }


    //--------------------------------------------------

    // initialize temp variable storages
    float_p lx1, ly1, lz1,     lx2, ly2, lz2;
    float_p ux1, uy1, uz1, w1, ux2, uy2, uz2, w2;
    float_p ux3, uy3, uz3, w3, ux4, uy4, uz4, w4;
    float_p e1, e2, e3, e4;
    float_p m3, m4;

    float_p wmin, wmax, prob;

    //--------------------------------------------------
    // loop over incident types
    for(auto&& con1 : tile.containers)
    {
      auto t1 = con1.type;
      if(is_empty(t1)) continue; // no interactions with incident type t1

      size_t Ntot1 = con1.size();

      // create randomized order for particle access
      //std::vector<size_t> inds(Ntot1);
      //for(size_t n1=0; n1<Ntot1; n1++) inds[n1] = n1;
      //std::shuffle(std::begin(inds), std::end(inds), gen);

      // loop over incident particles
      //UniIter::iterate([=] DEVCALLABLE (
      //          size_t n, 
      //          pic::ParticleContainer<D>& con
      //          ){
      for(size_t n1=0; n1<Ntot1; n1++) {
      //for(int n1=con1.size()-1; n1>=0; n1--) {

        //unpack incident 
        lx1 = con1.loc(0,n1);
        ly1 = con1.loc(1,n1);
        lz1 = con1.loc(2,n1);
          
        ux1 = con1.vel(0,n1);
        uy1 = con1.vel(1,n1);
        uz1 = con1.vel(2,n1);
        w1  = con1.wgt(n1);

        e1  = con1.eneArr[n1]; // TODO or compute here? or call function in container?

        if(w1 < EPS) continue; // omit zero-w incidents

        //pre-calculate maximum partial interaction rates
        comp_pmax(t1, e1, cons); 

        // maximum interaction rate
        //= sigma_max * w2_sum * w1/prob_norm
        double prob_vir_max = 0.0;
        for(size_t i=0; i<ids.size(); i++) prob_vir_max += 2.0*cmaxs[i]*wsums[i]/prob_norm;
        // TODO: no w1 here ???

        if(prob_vir_max>=1.0){
          std::cout << " prob_vir_max:" << prob_vir_max << std::endl;
          std::cout << " norm" << prob_norm << std::endl;
          for(size_t i=0; i<ids.size(); i++) std::cout << cmaxs[i] << " " << wsums[i] << std::endl;

          assert(false);
        }


        // exponential waiting time between interactions
        double t_free = -log( rand() )/prob_vir_max;

        //if np.random.rand() < prob_max: # virtual maximum channel
        if(t_free < 1.0)  // virtual maximum channel
        { 
          // get random interaction
          // NOTE: incident type t1 must be the same as what comp_pmax was called with
          int i = draw_rand_proc(); // i:th interaction in probs array
          
          //--------------------------------------------------
          // unpack interaction
          auto jmin = jmins[i];
          auto jmax = jmaxs[i];
          auto wsum = wsums[i];
          auto cmax = cmaxs[i];

          auto int_id = ids[i];             // id in global array
          auto iptr = interactions[int_id]; // pointer to interaction 
          auto t2   = iptr->t2;             // target type
          auto con2 = cons[t2];             // target container

          //--------------------------------------------------
          // get random target with energy between jmin/jmax
          // propability is proptional to weight of LPs

          size_t n2 = toolbox::sample_prob_between(con2->wgtCumArr, rand(), jmin, jmax);

          if( (t1 == t2) && (n1 == n2) ) continue; // do not interact with self

          // unpack target
          lx2 = con2->loc(0,n2);
          ly2 = con2->loc(1,n2);
          lz2 = con2->loc(2,n2);

          ux2 = con2->vel(0,n2);
          uy2 = con2->vel(1,n2);
          uz2 = con2->vel(2,n2);
          w2  = con2->wgt(  n2);
          e2  = con2->eneArr[n2]; 

          if(w2 < EPS) continue; // omit zero-w targets
                                   
          //--------------------------------------------------
          // double counting prevention by considering only e1 < e2 cases
          //if(e1 > e2) continue;
                  
          //avoid double counting by considering only e1 > e2 cases
          //if(e1 < e2) continue;
          //could also put double counting check to comp_pmax
          //--------------------------------------------------
          
          wmin = min(w1, w2);
          wmax = max(w1, w2);

          //--------------------------------------------------
          
          // real probablity of interaction
          auto [cm, vrel] = iptr->comp_cross_section(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2 );

          // collect max cross section
          float_p cm_cur = info_max_int_cs[iptr->name];
          info_max_int_cs[iptr->name] = std::max( cm_cur, cm);


          // FIXME which max should this prob be compared against?  I think ver3 or 4

          // ver1
          // maximum partial interaction rate of current interaction type;
          // NOTE: should be sum of all interactions
          //double prob_int_max = 2.0*cmax*wsum*w1/prob_norm; // fac 2 from v_rel = 2c
          //prob = cm*wsum*w1/prob_norm; 
          //double prob_vir = prob/prob_int_max

          // ver2
          //prob = cm*wsum*w1/prob_norm; 
          //double prob_vir = prob/prob_vir_max;

          // ver3
          //double cm_hat_max = 0.0;
          //for(size_t i=0; i<ids.size(); i++) cm_hat_max += 2.0*cmaxs[i];
          //double prob_vir = cm/cm_hat_max;

          // ver4
          // comparison of interaction to max interaction
          double prob_vir = cm*vrel/(2.0*cmax);

          // FIXME remove check if sure this works
          if(prob_vir >= 1.0){
            std::cout << " prob_vir > 1: " << prob_vir << std::endl;
            std::cout << " int  " << iptr->name << std::endl;
            std::cout << " cm   " << cm << std::endl;
            std::cout << " vrel " << vrel << std::endl;
            std::cout << " cmax " << cmax << std::endl;

            assert(false);
          }

          if(rand() < prob_vir)  // check if this interaction is chosen among the samples ones
          {
            // particle values after interaction
            auto [t3, ux3, uy3, uz3, w3] = duplicate_prtcl(t1, ux1, uy1, uz1, w1);
            auto [t4, ux4, uy4, uz4, w4] = duplicate_prtcl(t2, ux2, uy2, uz2, w2);

            // interact and udpate variables in-place
            iptr->interact( t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4);

            // new energies
            m3 = (t3 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
            m4 = (t4 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
            e3 = std::sqrt( m3*m3 + ux3*ux3 + uy3*uy3 + uz3*uz3 );
            e4 = std::sqrt( m4*m4 + ux4*ux4 + uy4*uy4 + uz4*uz4 );

            // weight adaptation
            float_p fw3 = ene_weight_funs(t1, e1)/ene_weight_funs(t3, e3); 
            float_p fw4 = ene_weight_funs(t2, e2)/ene_weight_funs(t4, e4); 

            // limit particle creation
            float_p n3 = std::min( fw3, 32.0f );
            float_p n4 = std::min( fw4, 32.0f );


            //# NOTE these two expressions are equal: w_i = w_j/wmax == wmin/w_i
            //# this is where the algorithm differs from original LP MC method by Stern95;
            //# here, instead of killing the LP we do not update its energy.
            float_p prob_upd3 = wmin/w1; //w2/wmax;
            float_p prob_upd4 = wmin/w2; //w1/wmax;

            // redistribute weights among the new copies
            w3 = w1/n3;
            w4 = w2/n4;


            //-------------------------------------------------- 
            int n_added = 0, ncop = 0;

            if(t1 == t3) { // same type before/after interactions; update energy with prob_upd

              double z1 = rand();
              while(n3 > z1 + ncop) {

                // optimized routine that does not to leave holes in arrays 
                // it first replaces the original value and only then adds if necessary
                if(ncop == 0) {
                  if( rand() < prob_upd3 ) {
                    //cons[t1].replace(iold, enew, wnew) //# replace with new energy
                                                         //
                    // NOTE: we keep location the same
                    con1.vel(0, n1) = ux3;
                    con1.vel(1, n1) = uy3;
                    con1.vel(2, n1) = uz3;
                    con1.wgt(   n1) = w3;
                    con1.eneArr[n1] = e3;
                  } else {
                    //cons[told].replace(iold, eold, wnew) # replace with old energy
                    con1.wgt(n1) = w3;
                  }
                } else {
                  if( rand() < prob_upd3 ) {
                    //cons[tnew].add(enew, wnew) # add new energy
                    con1.add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w3); // new ene & w
                  } else {
                    //con1->add(eold, wnew) # add old energy
                    con1.add_particle( {{lx1, ly1, lz1}}, {{ux1, uy1, uz1}}, w3); // new w
                  }
                }

                ncop += 1;
              } // end of while

              // remove parent prtcl if nothing was added
              if( ncop == 0 ) {
                //con1.delete(iold)
                con1.to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU version
                con1.wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              }

            //--------------------------------------------------
            } else { //# different before/after type; kill parent with a prob_upd

              double z1 = rand();
              while( n3 > z1 + ncop ){
                //cons[t3].add(enew, wnew)
                //#cons[tnew].buffer_add(enew, wnew)
                cons[t3]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w3); // new ene & w

                ncop += 1;
              }

              // kill parent
              if( prob_upd3 > rand() ) {
                //cons[told].delete(iold)
                //#cons[told].buffer_del(iold)

                con1.to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU version
                con1.wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              }

            } // end of prtcl t1/t3 addition


            //-------------------------------------------------- 
            // add prtcl t2/t4
            n_added = 0; 
            ncop = 0;

            if(t2 == t4) { // same type before/after interactions; update energy with prob_upd

              double z1 = rand();
              while(n4 > z1 + ncop) {

                // optimized routine that does not to leave holes in arrays 
                // it first replaces the original value and only then adds if necessary
                if(ncop == 0) {
                  if( rand() < prob_upd4 ) {
                    //cons[t1].replace(iold, enew, wnew) //# replace with new energy
                                                         //
                    // NOTE: we keep location the same
                    con2->vel(0, n2) = ux4;
                    con2->vel(1, n2) = uy4;
                    con2->vel(2, n2) = uz4;
                    con2->wgt(   n2) = w4;
                    con2->eneArr[n2] = e4;
                  } else {
                    //cons[told].replace(iold, eold, wnew) # replace with old energy
                    con2->wgt(n2) = w4;
                  }
                } else {
                  if( rand() < prob_upd4 ) {
                    //cons[tnew].add(enew, wnew) # add new energy
                    con2->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w4); // new ene & w
                  } else {
                    //con1->add(eold, wnew) # add old energy
                    con2->add_particle( {{lx2, ly2, lz2}}, {{ux2, uy2, uz2}}, w4); // new w
                  }
                }

                ncop += 1;
              } // end of while

              // remove parent prtcl if nothing was added
              if( ncop == 0 ) {
                //con1.delete(iold)
                con2->to_other_tiles.push_back( {0,0,0,n2} ); // NOTE: CPU version
                con2->wgt(n2) = 0.0f; // make zero wgt so its omitted from loop
              }

            //--------------------------------------------------
            } else { //# different before/after type; kill parent with a prob_upd

              double z1 = rand();
              while( n4 > z1 + ncop ){
                //cons[t3].add(enew, wnew)
                //#cons[tnew].buffer_add(enew, wnew)
                cons[t4]->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w4); // new ene & w

                ncop += 1;
              }

              // kill parent
              if( prob_upd4 > rand() ) {
                //cons[told].delete(iold)
                //#cons[told].buffer_del(iold)

                con2->to_other_tiles.push_back( {0,0,0,n2} ); // NOTE: CPU version
                con2->wgt(n2) = 0.0f; // make zero wgt so its omitted from loop
              }

            } // end of prtcl t1/t3 addition





            //--------------------------------------------------
            // old prtcl add routine
            //if(false) {

            //  if(rand() < prob_upd3){
            //    if(t1 == t3){ // if type is conserved only update the prtcl info
            //                    
            //      // NOTE: we keep location the same
            //      con1.vel(0,n1) = ux3;
            //      con1.vel(1,n1) = uy3;
            //      con1.vel(2,n1) = uz3;
            //    } else { // else destroy previous and add new 

            //      // destroy current
            //      con1.to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU version
            //      con1.wgt(n1) = 0.0f; // make zero wgt so its omitted from loop

            //      // add new
            //      cons[t3]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w1);

            //      //std::cout << "killing t1" << t1 << std::endl;
            //      //std::cout << "adding t3" << t3 << std::endl;
            //    }
            //  }
            //  //-------------------------------------------------- 

            //  //-------------------------------------------------- 
            //  if(rand() < prob_upd4){
            //    if(t2 == t4){ // if type is conserved only update the prtcl info

            //      // NOTE: we keep location the same
            //      con2->vel(0,n2) = ux4;
            //      con2->vel(1,n2) = uy4;
            //      con2->vel(2,n2) = uz4;
            //    } else { // else destroy previous and add new 

            //      // destroy current
            //      con2->to_other_tiles.push_back( {0,0,0,n2} ); // NOTE: CPU version
            //      con2->wgt(n2) = 0.0f; // make zero wgt so its omitted from loop

            //      cons[t4]->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w2);

            //      //std::cout << "killing t2" << t2 << std::endl;
            //      //std::cout << "adding t4" << t4 << std::endl;
            //    }
            //  }
            //}
            //--------------------------------------------------


            //-------------------------------------------------- 
          }
        }
      } // end of loop over con1 particles
      //}, con.size(), con);
      //UniIter::sync();
    }// end con1


    //--------------------------------------------------
    // info for bookkeeping
    for(auto&& con : tile.containers) {
      auto t1 = con.type;
      info_prtcl_num[t1] = con.size() - info_prtcl_num[t1]; // change in prtcl num
    }


    // calculate how many will be killed
    std::map<std::string, int> info_prtcl_kill;
    for(auto&& con : tile.containers) {
      auto t1 = con.type;
      info_prtcl_kill[t1] = cons[t1]->to_other_tiles.size();
    }

    // print
    //std::cout << "-------prtcl statistics----------\n";
    //for(auto const& [key, val] : info_prtcl_num){
    //  int n_kill = info_prtcl_kill[key];
    //  std::cout << "    " << key << " Delta N_p: " << val << " N_kill: " << n_kill << std::endl;
    //}

    //int n_tot = info_prtcl_num["ph"] - (info_prtcl_kill["e-"] + info_prtcl_kill["e+"]);
    //assert(n_tot == 0);

    //--------------------------------------------------
    //std::cout << " max cross sections:" << std::endl;
    //for(auto const& [key, val] : info_max_int_cs) std::cout << "   " << key << " : " << val << std::endl;

    //--------------------------------------------------
    // final book keeping routines
    for(auto&& con : tile.containers)
    {
      con.delete_transferred_particles(); // remove annihilated prtcls; this transfer storage 
                                          // is used as a tmp container for storing the indices
    }

    return;
  }

};



} // end of namespace qed
