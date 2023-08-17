#pragma once

#include <algorithm>
#include <string>
#include <tuple>
#include <random>
#include <memory>
#include <map>
#include <functional>
#include <cmath>

#include "interactions/interaction.h"
#include "../../definitions.h"
#include "../../pic/tile.h"
#include "../../tools/sample_arrays.h"
#include "../../tools/linlogspace.h"



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
template<size_t D>
class Pairing
{
private:
  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float_p> uni_dis;

  using InteractionPtr = std::shared_ptr<qed::Interaction>;

  // using raw pointer instead of smart ptrs; it does not take ownership of the object
  // so the container is not deleted when the temporary storage goes out of scope.
  using ConPtr = pic::ParticleContainer<D>* ;
  //using ConPtr = std::weak_ptr< pic::ParticleContainer<D> >; // this could also work
  //using ConPtr = std::reference_wrapper< pic::ParticleContainer<D> >; // this maybe as well

public:

  // constructor with incident/target types
  Pairing() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)       
  { 
    update_hist_lims(hist_emin, hist_emax, hist_nbin);
  }

  //using Tile_map = std::unordered_map<TileID_t, Tileptr>;
  //Tile_map tiles; /// Map with tile_id & tile data

  std::vector<InteractionPtr> interactions;

  // normalization factor for probabilities
  float_p prob_norm = 1.0f;

  // force e- e+ to have unit weights irrespective of weighting functions
  bool force_ep_uni_w = true; 

  //--------------------------------------------------
  //histogram for the leaking/escaping photons
  double hist_emin = -4.0; // log10(emin)
  double hist_emax =  2.0; // log10(emax)
  int    hist_nbin = 200;

  double leaked_ene  = 0.0;
  double leaked_wsum = 0.0;
  int    leaked_pnum = 0;

  double inj_ene_ph  = 0.0;
  double inj_ene_ep  = 0.0;

  double tau_measured = 0.0;
  
  ManVec<double> hist;
  ManVec<double> hist_ene_edges;

  void update_hist_lims(double emin, double emax, int nbin)
  {
    hist_emin = emin;
    hist_emax = emax;
    hist_nbin = nbin;

    // resize histogram storate
    hist.resize(nbin);
    hist_ene_edges.resize(nbin);

    // get logspace edges
    toolbox::logspace(hist_emin, hist_emax, hist_nbin, hist_ene_edges);   

    //std::cout << "hist init:" << std::endl;
    //for(size_t i=0; i<hist_nbin; i++) std::cout << "   " << i << " " << hist_ene_edges[i] << std::endl;
  }

  // update histogram
  void add_to_histogram(double x, double w)
  {
    int i = toolbox::find_sorted_nearest(hist_ene_edges, x);
    hist[i] += w;

    //std::cout << "adding to hist " << i << " x" << x << " w:" << w << std::endl;
  }

  // refill histogram with zeros 
  void clear_hist()
  {
    for(size_t i=0; i<hist_nbin; i++) hist[i] = 0.0;

    leaked_ene  = 0.0;
    leaked_wsum = 0.0;
    leaked_pnum = 0;

    inj_ene_ph  = 0.0;
    inj_ene_ep  = 0.0;
  }


  //--------------------------------------------------
  // auxiliary containers for Monte Carlo sampling 
  std::vector<float_p> 
      probs,    // maximum probabillity
      wsums,    // sum over possible target weights
      cmaxs;    // maximum cross section
      //faccs,    // w-weighted average target's accumulation factor
                 
  std::vector<int> 
      jmins,    // minimum indices of prtcls that can particiapte in initeraction
      jmaxs;    // maximum -||-
                 
  std::vector<size_t> 
      ids;     // internal id of the interaction in the storage

  std::map<std::string, double> info_max_int_cs; // maximum cross s measured
                 
  std::map<std::string, double> info_int_nums; // number of interactions performed

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
    auto long_name = name + "_" + t1 + "_" + t2;

    std::cout << " adding: " << name << " of t1/t2 " << t1 << " " << t2 << std::endl;
    interactions.push_back(iptr);

    info_max_int_cs[name] = 0.0;
    info_int_nums[long_name] = 0.0;
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
  void comp_pmax(string t1, float_p e1, 
                 std::map<std::string, pic::ParticleContainer<D>*>& cons)
  {

    //size_t n_ints = interactions.size(); // get number of interactions

    probs.clear(); // maximum probabillity
    wsums.clear(); // sum over possible target weights
    cmaxs.clear(); // maximum cross section
    jmins.clear(); // minimum indices of prtcls that can particiapte in initeraction
    jmaxs.clear(); // maximum -||-
    //faccs.clear(); // w-weigthed average target's accumulation factor // NOTE: not needed; now in wsums
    ids.clear(); // internal id of the interaction in the storage
                 
    size_t id = 0;
    for(auto iptr : interactions){
      if(t1 == iptr->t1)
      {
        auto t2 = iptr->t2;      // get target
        auto con_tar = cons[t2]; // target container
        size_t N2 = con_tar->size(); // total number of particles

        float_p cross_max = iptr->cross_section; // maximum cross section (including x2 for head-on collisions)

        // NOTE: assumes that target distribution remains static for the duration of the time step.
        // In that case, no LP changes energy and the limits finding is ok.
        // This assumption is valid in the limit of small time step dt << mean time between interactions

        //#find the min/max interaction energies and the corresponding array indices jmin/jmax
        auto [emin, emax] = iptr->get_minmax_ene(t1, t2, e1);

        //--------------------------------------------------
        // double counting prevention since we only consider targets with energy more than incident particle

        if(true){

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
          emin = std::max(e1, emin); // double counting prevention; only consider LPs with energy larger than incident
          
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
        float_p wsum = 0.0f;
        for(size_t j=jmin; j<jmax; j++) wsum += con_tar->wgt( j ); // sum( w[jmin:jmax] )

        // FIXME: switch to this version as it is faster
        // ver2: assumes static targets; value calculated in the beginnig of loop

        float_p wsum2 = 0.0f;
          
        //--------------------------------------------------
        if(! iptr->do_accumulate){ // normal mode; no accumulation

          if(jmin < jmax) { // in the opposite case arrays dont span a range and so wsum2 = 0
            float_p wsum_min = jmin == 0 ? 0.0f : con_tar->wgtCumArr[jmin];
            wsum2 = con_tar->wgtCumArr[jmax-1] - wsum_min;
          }

        //--------------------------------------------------
        } else { // accumulate interactions; effectively reduces weight

          float_p w, wprev, e2, f;

          wprev = jmin == 0 ? 0.0 : con_tar->wgtCumArr[jmin-1];
          for(size_t j=jmin; j<jmax; j++) {
            
            // weight between [j-1, j]
            //w = con_tar->wgt(j); // real weight
            w = con_tar->wgtCumArr[j] - wprev; // calc via cum array 
            wprev = con_tar->wgtCumArr[j];

            e2 = con_tar->eneArr[j];
            auto [f1,f2] = iptr->accumulate(t1, e1, t2, e2);
            f = f1*f2; //std::max(f1,f2);  // FIXME is max ok here? or product?

            //facc += f; // w-weighted average
            //facc += w*f/wsum2; // w-weighted average // ver1
            wsum2 += w/f; 
          }
        }

        if(wsum2 > 1.0e6){
          std::cout << "ERROR: wsum" << std::endl;
          std::cout << " wsum:" << wsum << std::endl;
          std::cout << " wsum2:" << wsum2 << std::endl;
          std::cout << " jmin:" << jmin << std::endl;
          std::cout << " jmax:" << jmax << std::endl;
          //std::cout << " wsum_min:" << wsum_min << std::endl;
          assert(false);
        }


        //--------------------------------------------------
        // average accumulation factor

        //double facc = 0.0;
        //if(iptr->do_accumulate){
        //  double w, wprev, e2, f;
        //  wprev = jmin == 0 ? 0.0 : con_tar->wgtCumArr[jmin-1];
        //  for(size_t j=jmin; j<jmax; j++) {
        //    
        //    // weight between [j-1, j]
        //    //w = con_tar->wgt(j);
        //    w = con_tar->wgtCumArr[j] - wprev; // calc via cum array 
        //    wprev = con_tar->wgtCumArr[j];

        //    e2 = con_tar->eneArr[j];
        //    auto [f1,f2] = iptr->accumulate(t1, e1, t2, e2);
        //    facc = std::max(f1,f2); // FIXME is max ok here? or product?

        //    //facc += f; // w-weighted average
        //    facc += w*f/wsum2; // w-weighted average // ver1
        //    
        //    // or max
        //    //facc = std::max(facc, f);                        // ver2

        //    //std::cout << " facc: " << facc;
        //    //std::cout << " f: " << f;
        //    //std::cout << " w: " << w;
        //    //std::cout << " wp: " << wprev;
        //    //std::cout << " e1: " << e1;
        //    //std::cout << " e2: " << e2;
        //    //std::cout << std::endl;

        //  }

        //} else {
        //  facc = 1.0;
        //}
         

        ////--------------------------------------------------
        //// debug
        ////std::cout << "wsum " <<  wsum << " " << wsum2 << " jminmax " << jmin << " " << jmax << std::endl;
        ////assert( std::abs( wsum - wsum2 ) < EPS );
        //--------------------------------------------------
        //total sum of target weights
        float_p wtot = 0.0f;
        for(size_t j=0; j<N2; j++) wtot += con_tar->wgt( j ); // sum( w )
        if(wtot <= 0.0f) wtot = 1.0f; // guard for NaNs

        //--------------------------------------------------
        // maximum partial interaction rate
        float_p par_int_rate = 2.0f*cross_max * wsum2; // factor 2 comes from v_rel = 2c

        // update/accept interaction only if it has non-zero prob to occur
        if(par_int_rate > 0.0f) {

          // add stuff 
          probs.push_back(par_int_rate);
          cmaxs.push_back(cross_max);
          wsums.push_back(wsum2);
          //faccs.push_back(facc);

          jmins.push_back(jmin);
          jmaxs.push_back(jmax);

          ids.push_back(id); // update interaction numbering
                           
          //std::cout << "comp_pmax: " << id << " " << iptr->name << " " << t1 << "/" << t2;
          //std::cout << " cmaxs:" << cross_max;
          //std::cout << " probs:" << par_int_rate << " prob: " << par_int_rate/prob_norm;
          //std::cout << " wsum:" << wsum << " wsum/wtot: " << wsum/wtot;
          //std::cout << " wsum2:" << wsum2;
          //std::cout << " js:" << jmin << " " << jme << " " << jmax;
          //std::cout << " emin:" << emin << " e1:" << e1 << " emax:" << emax;
          //std::cout << " facc:" << facc;
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
  //
  // value corresponds to number of particles (of type t) produced with given energy x
  // e.g., constant value means that every interactions splits the particle into that many pieces
  float_p ene_weight_funs(std::string t, float_p x) 
  {

    // standard unit weights
    if(       t == "ph") { return 1.0; 
    } else if(t == "e-") { return 1.0; 
    } else if(t == "e+") { return 1.0; 
    }
      
    // photon emphasis
    //if(       t == "ph") { return std::pow(x/0.01f, 0.5f); 
    //} else if(t == "e-") { return 1.0; 
    //} else if(t == "e+") { return 1.0; 
    //}
      
    // FIXME not useful
    //if(       t == "ph") { return x > 1.0 ? 2.0 : 1.0; //std::pow(x, 0.2); //std::pow(x, -0.5); 
    //} else if(t == "e-") { return 1.0; //std::pow(x, +0.2);
    //} else if(t == "e+") { return 1.0; //std::pow(x, +0.2);
    //}
    //if(       t == "ph") { return std::pow(x, +1.0); 
    //} else if(t == "e-") { return std::pow(x, -1.0);
    //} else if(t == "e+") { return std::pow(x, -1.0);
    //}

    assert(false);
  }


  //--------------------------------------------------
  void solve(pic::Tile<D>& tile)
  {
      
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

                //e1  = con1.eneArr[n1]; 
                e1  = con1.get_prtcl_ene(n1);

                if(w1 < EPS) continue; // omit zero-w incidents

                // unpack target
                lx2 = con2.loc(0,n2);
                ly2 = con2.loc(1,n2);
                lz2 = con2.loc(2,n2);

                ux2 = con2.vel(0,n2);
                uy2 = con2.vel(1,n2);
                uz2 = con2.vel(2,n2);
                w2  = con2.wgt(  n2);

                //e2  = con2.eneArr[n2]; 
                e1 = con2.get_prtcl_ene(n2);

                if(w2 < EPS) continue; // omit zero-w targets

                //--------------------------------------------------
                //avoid double counting by considering only e1 < e2 cases
                //if(e1 > e2) continue;

                // interaction cross section
                auto [cm, vrel] = iptr->comp_cross_section(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2 );
                
                // interaction probability
                wmin = min(w1, w2);
                wmax = max(w1, w2);
                prob = cm*vrel*w1*w2;
                // NOTE: difference of all2all scheme is here where prob depends on w1*w2

                // exponential waiting time between interactions
                double t_free = -log( rand() )*prob_norm/prob;

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
  void solve_mc(pic::Tile<D>& tile)
  {
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

    // ver2

    // random shuffled indices of containers; makes iteration order random in every step
    //std::vector<std::string> t_inis = { "e-", "e+", "ph" };
    //std::shuffle(std::begin(t_inis), std::end(t_inis), gen);
    //for(auto t1 : t_inis)
    //{
    //  if(is_empty(t1)) continue; // no interactions with incident type t1
    //  auto con1 = cons[t1];

    for(auto&& con1 : tile.containers) // ver1
    {
      auto t1 = con1.type;
      if(is_empty(t1)) continue; // no interactions with incident type t1

      //size_t Ntot1 = con1.size();
      size_t Ntot1 = info_prtcl_num[t1]; // read particle number from here; 
                                         // it changes adaptively and routines assume non-evolving arrays

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

        e1 = con1.get_prtcl_ene(n1);

        if(w1 < EPS) continue; // omit zero-w incidents

        //pre-calculate maximum partial interaction rates
        comp_pmax(t1, e1, cons); 

        if(ids.size() == 0) continue; // no targets to interact with 

        // maximum interaction rate
        //= sigma_max * w2_sum * w1/prob_norm

        // TODO add 1/facc here
        double prob_vir_max = 0.0;
        //for(size_t i=0; i<ids.size(); i++) prob_vir_max += 2.0*cmaxs[i]*wsums[i]/faccs[i]; //*prob_norm;
        //for(size_t i=0; i<ids.size(); i++) prob_vir_max += 2.0*cmaxs[i]*std::max(wsums[i], static_cast<double>(w1))/faccs[i]; // FIXME max(w1, w2) version
        //for(size_t i=0; i<ids.size(); i++) prob_vir_max += 2.0*cmaxs[i]*wsums[i]*w1/faccs[i]; // FIXME additional w1 here
        for(size_t i=0; i<ids.size(); i++) prob_vir_max += 2.0f*cmaxs[i]*wsums[i]*w1; // FIXME facc is in wsums
                                                                                             
                                                                                             
        // NOTE: no w1 here. putting w1 gives different result from uni-weight sim. Hence, its wrong.
        // instead, it is taken into account later on via prob_update \propto 1/w1

        //if(prob_vir_max>=1.01){ // some tolerance here
        //  std::cout << " prob_vir_max:" << prob_vir_max << std::endl;
        //  std::cout << " norm:" << prob_norm << std::endl;
        //  for(size_t i=0; i<ids.size(); i++) std::cout << "c: " << cmaxs[i] << " w: " << wsums[i] << " f: " << faccs[i] << std::endl;
        //  //assert(false);
        //}

        // exponential waiting time between interactions
        double t_free = -log( rand() )*prob_norm/prob_vir_max; ///w1; // FIXME added /w1


        //if( t_free > 1.0) {
        //if( true ) {
        //  std::cout<< "t_free: " << t_free << " N_Q/p_int: " << prob_norm/prob_vir_max << std::endl;
        //  std::cout<< "prob_vir_max: " << prob_vir_max << std::endl;
        //  std::cout<< " t1 " << t1 << std::endl;
        //  std::cout<< " ids " << ids.size() << std::endl;
        //  for(size_t i=0; i<ids.size(); i++) {
        //  std::cout<< " ints " << i << " " << cmaxs[i] << " " << wsums[i] << " " << faccs[i] << std::endl;
        //  }
        //  //assert(false);
        //}

        //if np.random.rand() < prob_max: # virtual maximum channel
        if(t_free < 1.0) // virtual maximum channel
        { 
          // get random interaction
          // NOTE: incident type t1 must be the same as what comp_pmax was called with
          int i = draw_rand_proc(); // i:th interaction in probs array
          
          //--------------------------------------------------
          // unpack interaction
          auto jmin = jmins[i];
          auto jmax = jmaxs[i];
          auto cmax = cmaxs[i];
          //auto wsum = wsums[i];
          //auto facc_max = faccs[i];

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
          e2  = con2->get_prtcl_ene(n2);

          if(w2 < EPS) continue; // omit zero-w targets
                                   
          //--------------------------------------------------
          // double counting prevention by considering only e1 < e2 cases
          //if(e1 > e2) continue;
                  
          //avoid double counting by considering only e1 > e2 cases
          //if(e1 < e2) continue;
          //could also put double counting check to comp_pmax
          //--------------------------------------------------
          
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

          // correct average accumulation factor with the real value
          // TODO add 1/facc here?
          // or facc_real/facc_vir ??
          //auto [facc3, facc4] = iptr->do_accumulate ? iptr->accumulate(t1, e1, t2, e2) : 1.0;
          auto [facc3, facc4] = iptr->accumulate(t1, e1, t2, e2);
          facc3 = iptr->do_accumulate ? facc3 : 1.0f;
          facc4 = iptr->do_accumulate ? facc4 : 1.0f;


          // FIXME remove check if sure this works
          //if(true){
          if(prob_vir >= 1.0){
            std::cout << " prob_vir > 1: " << prob_vir << std::endl;
            std::cout << " int  " << iptr->name << std::endl;
            std::cout << " cm   " << cm << std::endl;
            std::cout << " vrel " << vrel << std::endl;
            std::cout << " cmax " << cmax << std::endl;
            std::cout << " facc3 " << facc3 << std::endl;
            std::cout << " facc4 " << facc4 << std::endl;
            //assert(false);
          }

          // correct for real/max facc factor
          //prob_vir *= facc_max/std::max(facc3, facc4); // NOTE: prob \propto 1/facc 
          // NOTE: not needed anymore since facc is in wsum -> \sum w/facc

          if(rand() < prob_vir)  // check if this interaction is chosen among the sampled ones
          {
            auto long_name = iptr->name + "_" + t1 + "_" + t2;
            info_int_nums[long_name] += 1;

            // particle values after interaction
            auto [t3, ux3, uy3, uz3, w3] = duplicate_prtcl(t1, ux1, uy1, uz1, w1);
            auto [t4, ux4, uy4, uz4, w4] = duplicate_prtcl(t2, ux2, uy2, uz2, w2);

            // interact and udpate variables in-place
            iptr->interact( t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4 );

            // new energies; NOTE: could use container.m to get the mass
            m3 = (t3 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
            m4 = (t4 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
            e3 = std::sqrt( m3*m3 + ux3*ux3 + uy3*uy3 + uz3*uz3 );
            e4 = std::sqrt( m4*m4 + ux4*ux4 + uy4*uy4 + uz4*uz4 );

            // weight adaptation; original Stern95 version
            //float_p fw3 = (e3/e1)*ene_weight_funs(t1, e1)/ene_weight_funs(t3, e3); 
            //float_p fw4 = (e4/e2)*ene_weight_funs(t2, e2)/ene_weight_funs(t4, e4); 

            // more intuitive version (flipped)
            float_p fw3 = ene_weight_funs(t3, e3)/ene_weight_funs(t1, e1);
            float_p fw4 = ene_weight_funs(t4, e4)/ene_weight_funs(t2, e2);

            // limit explosive particle creation
            float_p n3 = std::min( fw3, 32.0f );
            float_p n4 = std::min( fw4, 32.0f );

            //# NOTE these two expressions are equal: w_i = w_j/wmax == wmin/w_i
            //# this is where the algorithm differs from original LP MC method by Stern95;
            //# here, instead of killing the LP we do not update its energy.

            // we can think here that only wmin = min(w1, w2) fraction of the LPs weight is 
            // used in the interaction. If everything is used then (i.e., wmin=w1 for LP1 
            // or wmin=w2 for LP2)) then prob_upd = 1

            // NOTE: update probability is reduced as if particle weights are w -> w*f_acc

            wmin = min(w1, w2);
            wmax = max(w1, w2);

            // probability for the particle energy to be updated
            float_p prob_upd3 = wmin/w1; //w2/wmax;
            float_p prob_upd4 = wmin/w2; //w1/wmax;

            // probability for the particle to die
            float_p prob_kill3 = INF;
            float_p prob_kill4 = INF;

            //--------------------------------------------------
            // accumulation means that we treat the particles weight as w -> w*f_cc
            // but consider the interaction less often as prob -> prob/f_acc
            // to compensate, this means that we need to also genereate n_cop -> ncop/f_acc 
            // less of particles.
            //
            // NOTE: this increases the weight of the new LP3 automatically, i.e., we do not need 
            // to multiply w -> w*f_acc anymore.

            // TODO
            //n3 *= 1.0/facc3; // FIXME
            //n4 *= 1.0/facc4;
            //prob_upd3 *= facc3;

            // alternatively do not kill accumulated particles
            //prob_upd3 *= 1.0f/facc3;

            if( w3*facc3 > 1.0e8f || w4*facc4 > 1.0e8f) {
              std::cout << " ERR: pairing" << std::endl;
              std::cout << " w3: " << w3 << std::endl;
              std::cout << " w4: " << w4 << std::endl;
              std::cout << " e3: " << e3 << std::endl;
              std::cout << " e4: " << e4 << std::endl;
              std::cout << " n3: " << n3 << std::endl;
              std::cout << " n4: " << n4 << std::endl;
              std::cout << " f3: " << facc3 << std::endl;
              std::cout << " f4: " << facc4 << std::endl;

              assert(false);
            }



            //-------------------------------------------------- 
            int n_added = 0; //ncop = 0;
            double ncop = 0.0;

            if(t1 == t3) { // same type before/after interactions; update energy with prob_upd
                             
              // -------------scattering interactions go here----------------
              w3 = w1/n3; // redistribute weights among the new copies

              if(facc3 > 1.0f) {
                w3 *= facc3; // increase weight for accumulated interactions
                //prob_kill3 = facc3*w1/w2;
                //n3 *= 1.0f/facc3; // less incidents produced 
                //n3 *= w2/(w1*facc3);

                //n3 *= 1.0f/facc3; // FIXME why decrease numbers? to conserve energy, we should increase w
                //prob_upd3 = 1.0f; // FIXME always update, why?
              }

              if(force_ep_uni_w && (t3 == "e-" || t3 == "e+") ){
                w3 = 1.0f;
                n3 = facc3*w1/w3; // remembering to increase prtcl num w/ facc
              }

              // split too big LPs
              if(w3 > 1.0e3f) {
                w3 *= 0.5f;
                n3 *= 2.0f;
              }

            
              double z1 = rand();
              while(n3 > z1 + ncop) {

                // optimized routine that does not to leave holes in arrays 
                // it first replaces the original value and only then adds if necessary
                if(ncop < EPS) {
                  if( rand() < prob_upd3 ) {
                    // NOTE: we keep location the same
                    cons[t1]->vel(0, n1) = ux3;
                    cons[t1]->vel(1, n1) = uy3;
                    cons[t1]->vel(2, n1) = uz3;
                    cons[t1]->wgt(   n1) = w3;
                  } else {
                    cons[t1]->wgt(n1) = w3;
                  }
                } else {
                  if( rand() < prob_upd3 ) {
                    cons[t1]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w3); // new ene & w
                  } else {
                    cons[t1]->add_particle( {{lx1, ly1, lz1}}, {{ux1, uy1, uz1}}, w3); // new w
                  }
                }

                // ver2; same as above but different order of if statements
                //if( rand() < prob_upd3 ) {
                //  if(ncop == 0) {
                //    // NOTE: we keep location the same
                //    cons[t1]->vel(0, n1) = ux3;
                //    cons[t1]->vel(1, n1) = uy3;
                //    cons[t1]->vel(2, n1) = uz3;
                //    cons[t1]->wgt(   n1) = w3;
                //  } else { 
                //    cons[t1]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w3); // new ene & w
                //  }

                //} else { 
                //  if(ncop == 0) {
                //    cons[t1]->wgt(n1) = w3;
                //  } else {
                //    cons[t1]->add_particle( {{lx1, ly1, lz1}}, {{ux1, uy1, uz1}}, w3); // new w
                //  }
                //}

                ncop += 1.0;
              } // end of while

              // remove parent prtcl if nothing was added
              if( ncop < EPS ) {
                cons[t1]->to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU version
                cons[t1]->wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              }

              // FIXME remove?
              //if( rand() < 1.0f/prob_kill3) {
              //  cons[t1]->to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU version
              //  cons[t1]->wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              //}
              //if( ncop > 0 ) {
              //  float_p n5 = facc3*w1/w2;
              //}

            //--------------------------------------------------
            } else { //# different before/after type; kill parent with a prob_upd

              // annihilation interactions go here

              w3 = wmin/n3;

              if(force_ep_uni_w && (t3 == "e-" || t3 == "e+") ){
                w3 = 1.0;
                n3 = wmin/w3;
              }

              double z1 = rand();
              while( n3 > z1 + ncop ){
                cons[t3]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w3); // new ene & w
                ncop += 1.0;
              }

              // kill parent
              if( prob_upd3 > rand() ) {
                cons[t1]->to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU version
                cons[t1]->wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              }

            } // end of prtcl t1/t3 addition


            //-------------------------------------------------- 
            // add prtcl t2/t4
            n_added = 0; 
            ncop = 0.0;

            if(t2 == t4) { // same type before/after interactions; update energy with prob_upd

              // scattering interactions go here

              w4 = w2/n4; // redistribute weights among the new copies

              if(facc4 > 1.0f) {
                w4 *= facc4; // increase weight for accumulated interactions
              }

              if(force_ep_uni_w && (t4 == "e-" || t4 == "e+") ){
                w4 = 1.0f;
                n4 = facc4*w2/w4; // remembering to increase prtcl num w/ facc
              }
                
              // split too big LPs
              if(w4 > 1.0e3f) {
                //n4 = w4/1.0e3;
                //w4 = 1.0e3;
                w4 *= 0.5f;
                n4 *= 2.0f;
              }

              double z1 = rand();
              while(n4 > z1 + ncop) {

                // optimized routine that does not to leave holes in arrays 
                // it first replaces the original value and only then adds if necessary
                if(ncop < EPS) {
                  if( rand() < prob_upd4 ) {
                    // NOTE: we keep location the same
                    cons[t2]->vel(0, n2) = ux4;
                    cons[t2]->vel(1, n2) = uy4;
                    cons[t2]->vel(2, n2) = uz4;
                    cons[t2]->wgt(   n2) = w4;
                  } else {
                    cons[t2]->wgt(n2) = w4;
                  }
                } else {
                  if( rand() < prob_upd4 ) {
                    cons[t2]->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w4); // new ene & w
                  } else {
                    cons[t2]->add_particle( {{lx2, ly2, lz2}}, {{ux2, uy2, uz2}}, w4); // new w
                  }
                }

                ncop += 1.0;
              } // end of while

              // remove parent prtcl if nothing was added
              if( ncop < EPS ) {
                cons[t2]->to_other_tiles.push_back( {0,0,0,n2} ); // NOTE: CPU version
                cons[t2]->wgt(n2) = 0.0f; // make zero wgt so its omitted from loop
              }

            //--------------------------------------------------
            } else { //# different before/after type; kill parent with a prob_upd
                       
              // annihilation interactions go here

              w4 = wmin/n4;

              if( force_ep_uni_w && ( t4 == "e-" || t4 == "e+") ){
                w4 = 1.0;
                n4 = wmin/w4;
              }
            
              double z1 = rand();
              while( n4 > z1 + ncop ){
                cons[t4]->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w4); // new ene & w
                ncop += 1.0;
              }

              // kill parent
              if( prob_upd4 > rand() ) {
                cons[t2]->to_other_tiles.push_back( {0,0,0,n2} ); // NOTE: CPU version
                cons[t2]->wgt(n2) = 0.0f; // make zero wgt so its omitted from loop
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
    //std::cout << " interaction evaluations:" << std::endl;
    //for(auto const& [key, val] : info_int_nums) std::cout << "   " << key << " : " << val << std::endl;

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
  // normalize container of type t1
  void rescale(pic::Tile<D>& tile, string& t1, double f_kill)
  {

    // build pointer map of types to containers; used as a helper to access particle tyeps
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    cons[t1]->to_other_tiles.clear(); // clear book keeping array

    size_t N1 = cons[t1]->size(); // read particle number from here; 

    // total weight = sum(ws)
    float_p wtot = 0.0;
    for(size_t n1=0; n1<N1; n1++) wtot += cons[t1]->wgt(n1);

    // TODO it is not energy conserving to select particles equally w/o w-weighting

    // loop over particles
    float_p w1, zeta, prob_kill;
    for(size_t n1=0; n1<N1; n1++) {
      w1  = cons[t1]->wgt(n1);

      prob_kill = 1.0f - 1.0f/f_kill;

      zeta = rand();
      if( zeta < prob_kill) {
        cons[t1]->to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU deletion version
        cons[t1]->wgt(n1) = 0.0f; 
      } else {
        cons[t1]->wgt(n1) = w1*f_kill; // compensate lost particles by increasing w
      }
    }

    // remove annihilated prtcls; this transfer storage 
    cons[t1]->delete_transferred_particles(); 

    return;
  }

  // inject soft photons
  void inject_photons(pic::Tile<D>& tile, 
      float_p temp_inj, 
      float_p wph_inj,
      float_p Nph_inj) 
  {
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    auto mins = tile.mins;
    auto maxs = tile.maxs;

    float_p x_min = D >= 1 ? mins[0] : 0.0f;
    float_p y_min = D >= 2 ? mins[1] : 0.0f;
    float_p z_min = D >= 3 ? mins[2] : 0.0f;

    float_p lenx = D >= 1 ? (maxs[0] - mins[0]) : 0.0f;
    float_p leny = D >= 2 ? (maxs[1] - mins[1]) : 0.0f;
    float_p lenz = D >= 3 ? (maxs[2] - mins[2]) : 0.0f;

    // inject Nph_inj photons
    float_p vx, vy, vz, xi;
    float_p ux, uy, uz, xinj;
    float_p xloc, yloc, zloc;

    float_p ncop = 0.0f;
    float_p z1 = rand();

    while(Nph_inj > z1 + ncop)
    {
      //--------------------------------------------------
      // draw random isotropic 3d vector
      vz = 2.0f*rand() - 1.0f;
      xi = 2.0f*PI*rand();

      vx = std::sqrt(1.0f-vz*vz)*std::cos(xi);
      vy = std::sqrt(1.0f-vz*vz)*std::sin(xi);

      xloc = rand()*lenx + x_min;
      yloc = rand()*leny + y_min;
      zloc = rand()*lenz + z_min;

      //--------------------------------------------------
      // draw energy from a black body distribution
      float_p xi1 = rand();
      float_p xi2 = rand();
      float_p xi3 = rand();
      float_p xi4 = rand();
    
      float_p xi, jj, fsum;
      if( 1.202f*xi1 < 1.0f ){
          xi = 1.0f;
      } else {
          jj = 1.0f;
          fsum = std::pow(jj, -3);
          while( 1.202*xi1 > fsum + std::pow(jj + 1.0f, -3) )
          {
              jj   += 1.0f;
              fsum += std::pow(jj, -3);
          }
          xi = jj + 1.0f;
      }
      xinj = -temp_inj*std::log( xi2*xi3*xi4 )/xi;

      //--------------------------------------------------
      ux = xinj*vx;
      uy = xinj*vy;
      uz = xinj*vz;

      cons["ph"]->add_particle( {{xloc, yloc, zloc}}, {{ux, uy, uz}}, wph_inj );
      ncop += 1.0f;

      inj_ene_ph += wph_inj*xinj; // bookkeeping of injected photon energy
    }

    return;
  }


  // inject soft photons
  void inject_plaw_pairs(pic::Tile<D>& tile, 
      float_p slope, 
      float_p pmin,
      float_p pmax,
      float_p w_inj,
      float_p N_inj) 
  {

    //const float_p pmin = 10.0f;
    //const float_p pmax = 100.0f;

    assert(pmin > 1.0f); // pmin interpreted as gamma 
    assert(pmax > 1.0f); // pmax interpreted as gamma 

    float_p pminp = std::pow(pmin, 1.0f+slope);  // pmin^(1-g)
    float_p pmaxp = std::pow(pmax, 1.0f+slope);  // pmax^(1-g)

    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    auto mins = tile.mins;
    auto maxs = tile.maxs;

    float_p x_min = D >= 1 ? mins[0] : 0.0f;
    float_p y_min = D >= 2 ? mins[1] : 0.0f;
    float_p z_min = D >= 3 ? mins[2] : 0.0f;

    float_p lenx = D >= 1 ? (maxs[0] - mins[0]) : 0.0f;
    float_p leny = D >= 2 ? (maxs[1] - mins[1]) : 0.0f;
    float_p lenz = D >= 3 ? (maxs[2] - mins[2]) : 0.0f;

    // inject N_inj pairs
    float_p vx, vy, vz, xi;
    float_p ux, uy, uz, ginj, pinj;
    float_p xloc, yloc, zloc;

    float_p ncop = 0.0f;
    float_p z1 = rand();

    while(N_inj > z1 + ncop)
    {
      xloc = rand()*lenx + x_min;
      yloc = rand()*leny + y_min;
      zloc = rand()*lenz + z_min;

      // loop over both species
      for(size_t t=0; t<2; t++){

        //--------------------------------------------------
        // draw random isotropic 3d vector
        vz = 2.0f*rand() - 1.0f;
        xi = 2.0f*PI*rand();

        vx = std::sqrt(1.0f-vz*vz)*std::cos(xi);
        vy = std::sqrt(1.0f-vz*vz)*std::sin(xi);

        //--------------------------------------------------
        // get random powerlaw energy
        ginj = std::pow( pminp + rand()*( pmaxp - pminp ), 1.0f/(1.0f+slope));

        // u = sqrt( gamma^2 - 1)
        pinj = std::sqrt(ginj*ginj - 1.0f);

        //--------------------------------------------------
        ux = pinj*vx;
        uy = pinj*vy;
        uz = pinj*vz;

        if(t==0) cons["e-"]->add_particle( {{xloc, yloc, zloc}}, {{ux, uy, uz}}, w_inj );
        if(t==1) cons["e+"]->add_particle( {{xloc, yloc, zloc}}, {{ux, uy, uz}}, w_inj );

        inj_ene_ep += w_inj*ginj; // bookkeeping of injected photon energy
      }

      ncop += 1.0f;
    }

    return;
  }


  //--------------------------------------------------
  // photon escape from the box w/ escape probability formalism
  void leak_photons(
      pic::Tile<D>& tile, 
      double w2tau_units,
      double tc_per_dt,
      double tau_ext
      )
  {

    // build pointer map of types to containers; used as a helper to access particle tyeps
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //--------------------------------------------------
    // pair number density
    float_p wsum_ep = 0.0f;
    float_p wvrel   = 0.0f;
    float_p w1, beta, gam;

    size_t Ne = cons["e-"]->size(); 
    for(size_t n1=0; n1<Ne; n1++) {
      w1  = cons["e-"]->wgt(n1);
      gam = cons["e-"]->get_prtcl_ene(n1);

      beta = std::sqrt(1.0f - 1.0f/(gam*gam) );
      wvrel += beta*w1;
      wsum_ep += w1;
    }

    size_t Np = cons["e+"]->size(); 
    for(size_t n1=0; n1<Np; n1++) {
      w1  = cons["e+"]->wgt(n1);
      gam = cons["e+"]->get_prtcl_ene(n1);

      beta = std::sqrt(1.0f - 1.0f/(gam*gam) );
      wvrel += beta*w1;
      wsum_ep += w1;
    }
    //--------------------------------------------------
    // TODO which one?
    float_p tauT  = wvrel*w2tau_units; // \tau_T = \sigma_T <v_rel> wsum
    float_p tauT2 = wsum_ep*w2tau_units; // \tau_T = \sigma_T wsum

    tau_measured = tauT; // store for book keeping

    //std::cout << "escape" << std::endl;
    //std::cout << "   tauT: " << tauT << std::endl;
    //std::cout << "   tauT2:" << tauT2 << std::endl;
    //std::cout << "   wsum: " << wsum_ep << std::endl;
    //std::cout << "   N_-:  " << cons["e-"]->size() << std::endl;
    //std::cout << "   N_+:  " << cons["e+"]->size() << std::endl;
    //std::cout << "   N_w:  " << w2tau_units << std::endl;

    if(tau_ext > 0.0) tauT = tau_ext; // use external tau if given

    //--------------------------------------------------
    std::string t1 = "ph";
    size_t Nx = cons[t1]->size(); // read particle number from here; 
                                  //
    cons[t1]->to_other_tiles.clear(); // clear book keeping array

    float_p w, x, f, sKN, P_esc;
    for(size_t n1=0; n1<Nx; n1++) {
      w = cons[t1]->wgt(n1);
      x = cons[t1]->get_prtcl_ene(n1);

      // Klein-Nishina cross section
      //sKN = 0.75f*( 
      //    ( (1.0f + x)/x*x*x)*( 2.0f*x*(1.0f + x)/(1.0f + 2.0f*x) - std::log(1.0f + 2.0f*x)) 
      //      + std::log(1.0f + 2.0f*x)/(2.0f*x)  - (1.0f + 3.0f*x)/pow(1.0f + 2.0f*x, 2) );
      sKN = 1.0f; // Thomson cross-section
                  //
      // TODO sometimes gives sKN < 0.0 values

      // empirical escape probability function to account for forward scattering pile-up
      // see Lightman \& Zdarskiaki 1987
      if(        x <= 0.1) { f = 1.0f;
      } else if (x >  1.0) { f = 0.0f;
      } else               { f = (1.0f - x)/0.9f; }

      //(c/R)*dt = dt/t_c
      //P_esc = dt_per_tc/( 0.75f + 0.188f*tauT*f*sKN ); // sphere
      float_p t_esc = ( 1.0f + tauT*f*sKN ); // slab
      //float_p t_esc = ( 1.0f + tauT*f*sKN + (1.0f-f)*2.886f ); // slab w/ asymptic scaling to 5/sqrt(3)
                                                               // this mimics pair-production opacity
                                                               // asymptotic solution to slab geometry
                                                               // with rad. transf. when tau -> infty


      //P_esc = t_esc/dt_per_tc; // P = R/c / dt for tau -> 0 

      // tc_per_dt = 1/x = 20

      // photon has an escape probability rate of P_esc = 1.0/(t_c*t_esc) to escape
      // probability is p_esc = P_esc*dt = dt/(t_c*t_esc) 

      // FIXME
      //tc_per_dt *= w; // compensate by weight; heavier prtcl has less prob of escaping

      //std::cout << "esc:" << 1.0f/tc_per_dt/t_esc << " " << t_esc << " " << tc_per_dt << std::endl;

      if( w > 1.0e15 or x > 1.0e4) {
        std::cout << "ERR: leak ph" << std::endl;
        std::cout << "  x:" << x << std::endl;
        std::cout << "  w:" << w << std::endl;
        std::cout << "  f:" << f << std::endl;
        std::cout << "  t_esc:" << t_esc << std::endl;
        std::cout << "  sKN   :" << sKN << std::endl;
        std::cout << "  1/tc_per_dt :" << 1.0f/tc_per_dt << std::endl;
        assert(false);
      }


      if( 1.0f/tc_per_dt/t_esc > rand() ) {
      //if( 1.0f/(tc_per_dt*t_esc) > rand() ) {
      //if( 1.0f/rand() > tc_per_dt*t_esc) {

        // escape
        //if( t_esc/tc_per_dt > rand() ){
        //if( 1.0f/(t_esc*tc_per_dt) > rand() ){
        //if( rand() > t_esc*tc_per_dt ){
        //if( rand() < dt_per_tc/t_esc ){
        //if( rand() < 1.0f/P_esc ){
        //if( P_esc  < 1.0f/rand() ){
        //if( 1.0/rand() > P_esc ){

        leaked_ene  += x*w;
        leaked_wsum += w;
        leaked_pnum += 1;

        // book keeping of escaped flux
        add_to_histogram(x, w);

        cons[t1]->to_other_tiles.push_back( {0,0,0,n1} ); // NOTE: CPU deletion version
        cons[t1]->wgt(n1) = 0.0f; 
      }

    }

    //// remove annihilated prtcls; this transfer storage 
    cons[t1]->delete_transferred_particles(); 

    return;
  }


};



} // end of namespace qed
