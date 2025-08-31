#pragma once

#include <algorithm>
#include <string>
#include <tuple>
#include <random>
#include <memory>
#include <map>
#include <functional>
#include <cmath>

#include "definitions.h"
#include "core/pic/tile.h"
#include "tools/sample_arrays.h"
#include "tools/linlogspace.h"
#include "tools/staggered_grid.h"

#ifdef DEBUG
#define USE_INTERNAL_TIMER // comment this out to remove the profiler
#endif

#include "external/timer/timer.h"

#include "core/qed/interactions/interaction.h"


namespace qed {
  using std::string;
  using std::tuple;
  using std::min;
  using std::max;

  using toolbox::shape; // tanh function


// duplicate particle info into fresh variables
inline auto duplicate_prtcl(
    string t1, float ux1, float uy1, float uz1, float w1
    ) -> std::tuple<string, float, float, float, float>
{
  return {t1, ux1, uy1, uz1, w1};
}

//--------------------------------------------------
// Monte Carlo pairing of particles
//
// This monolithic object implements an advanced Monte Carlo pairing of:
//
// 1) Single body interactions; here processes proceed as
//      p1 -> p3 + (p4) 
//    where p4 is optional 
//
// 2) Two-body binary interactions; here processes proceed as
//    p1 + p2 -> p3 + p4
//
// The binary pairing proceeds in steps, as outlined in Stern et al. 1995. Most importantly,
// we rely on:
//
// 1) sampling the processes via a virtual interaction channel (sum of all interactions)
//    and then, if the virtual channel is selected, we sample which real process is performed.
//    This avoids calculating the cross sections for each real physical process and instead requires 
//    only the maximum cross section to be calculated once for each particle. Implementation is in
//    comp_pmax()
//
// 2) Automated re-weighting of particles during the interaction. The weights can be fine-tuned via
//    the functions given in ene_weight_funs() (NOTE: these are inverse of values described in Stern95).
//    
template<size_t D>
class Pairing
{
private:
  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<float> uni_dis;

  using InteractionPtr = std::shared_ptr<qed::Interaction>;

  // using raw pointer instead of smart ptrs; it does not take ownership of the object
  // so the container is not deleted when the temporary storage goes out of scope.
  using ConPtr = pic::ParticleContainer<D>* ;
  //using ConPtr = std::weak_ptr< pic::ParticleContainer<D> >; // this could also work
  //using ConPtr = std::reference_wrapper< pic::ParticleContainer<D> >; // this maybe as well

public:

  Timer timer; // internal timer for profiling

  // constructor with incident/target types
  Pairing() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0),
    timer("qed pairing")
  { 
    update_hist_lims(hist_emin, hist_emax, hist_nbin);

    timer.do_print = true;
    timer.verbose = 0;
  }

  //using Tile_map = std::unordered_map<TileID_t, Tileptr>;
  //Tile_map tiles; /// Map with tile_id & tile data

  // one-body single interactions
  std::vector<InteractionPtr> single_interactions;

  // two-body binary interactions
  std::vector<InteractionPtr> binary_interactions;

  // normalization factor for two-body interaction probabilities
  float prob_norm = 1.0f;

  // normalization factor for single-body interaction probabilities
  float prob_norm_onebody = 1.0f;


  // force e- e+ to have unit weights irrespective of weighting functions
  bool force_ep_uni_w = true; 

  // maximum number of particles allowed in a tile; no addition above this
  // NOTE: this means energy is not conserved if we skip prtcl addition
  int max_tile_prtcl_num = 1000000000; 
  int max_tile_phot_num = 1000000000; 

  //--------------------------------------------------
  // optional virtual field component (esp. for 1D sims to mimic varying backgrounds)

  bool use_vir_curvature = false;
  float vir_pitch_ang = 0.0; // sin\alpha/\gamma of the (virtual) curvature pitch angle
                               
  float r_curv = 1.0f; // curvature radius
  float r_gap = 100.0f; // gap length 

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

  double tau_global = 0.0;
  
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
    for(int i=0; i<hist_nbin; i++) hist[i] = 0.0;

    leaked_ene  = 0.0;
    leaked_wsum = 0.0;
    leaked_pnum = 0;

    inj_ene_ph  = 0.0;
    inj_ene_ep  = 0.0;
  }


  //--------------------------------------------------
  // auxiliary containers for Monte Carlo sampling 
  std::vector<float> 
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
  float rand() { return uni_dis(gen); };
  
  // add interactions to internal memory of the class; 
  // done via pointers to handle pybind interface w/ python
  void add_interaction(InteractionPtr iptr)
  {
    assert(iptr); // check that we are not appending nullptr

    //-------------------------------------------------- 
    if( iptr->interaction_order == 1 ){ // single-body interactions

      auto name = iptr->name;
      auto t1 = iptr->t1;
      auto long_name = name + "_" + t1;

      single_interactions.push_back(iptr);

    //-------------------------------------------------- 
    } else if( iptr->interaction_order == 2 ){ // two-body binary interactions

      auto name = iptr->name;
      auto t1 = iptr->t1;
      auto t2 = iptr->t2;
      auto long_name = name + "_" + t1 + "_" + t2;

      //std::cout << " adding: " << name << " of t1/t2 " << t1 << " " << t2 << std::endl;
      binary_interactions.push_back(iptr);

      // additionall arrays
      info_max_int_cs[name] = 0.0;
      info_int_nums[long_name] = 0.0;
    }

  }

  //--------------------------------------------------
  // check if interaction list is empty for type t1

  bool is_empty_single_int(string& t1) 
  {
    int i=0;
    for(auto iptr : single_interactions){
      if(t1 == iptr->t1) i += 1;
    }
    return (i > 0) ? false : true; 
  }

  bool is_empty_binary_int(string& t1) 
  {
    int i=0;
    for(auto iptr : binary_interactions){
      if(t1 == iptr->t1) i += 1;
    }
    return (i > 0) ? false : true; 
  }


  // compute maximum partial interaction rates for each process 
  // that LP of type t1 and energy of e1 can experience.
  void comp_pmax(string t1, float e1, 
                 std::map<std::string, pic::ParticleContainer<D>*>& cons)
  {

    //size_t n_ints = interactions.size(); // get number of interactions
    probs.clear(); // maximum probabillity
    wsums.clear(); // sum over possible target weights
    cmaxs.clear(); // maximum cross section
    jmins.clear(); // minimum indices of prtcls that can particiapte in initeraction
    jmaxs.clear(); // maximum -||-
    ids.clear(); // internal id of the interaction in the storage
    //faccs.clear(); // w-weigthed average target's accumulation factor // NOTE: not needed; now in wsums
                 
    //std::fill( probs.begin(), probs.end(), 0); // maximum probabillity
    //std::fill( wsums.begin(), wsums.end(), 0); // sum over possible target weights
    //std::fill( cmaxs.begin(), cmaxs.end(), 0); // maximum cross section
    //std::fill( jmins.begin(), jmins.end(), 0); // minimum indices of prtcls that can particiapte in initeraction
    //std::fill( jmaxs.begin(), jmaxs.end(), 0); // maximum -||-
    //std::fill(   ids.begin(),   ids.end(), 0); // internal id of the interaction in the storage


    size_t id = 0;
    for(auto iptr : binary_interactions){

      if(t1 == iptr->t1)
      {
        const auto& t2 = iptr->t2;      // get target
        const auto  con_tar = cons[t2]; // target container
        //const size_t N2 = con_tar->size(); // total number of particles

        // skip containers with zero targets
        if(con_tar->eneArr.size() == 0) continue;

        const float cross_max = iptr->cross_section; // maximum cross section (including x2 for head-on collisions)

        // NOTE: assumes that target distribution remains static for the duration of the time step.
        // In that case, no LP changes energy and the limits finding is ok.
        // This assumption is valid in the limit of small time step dt << mean time between interactions

        //#find the min/max interaction energies and the corresponding array indices jmin/jmax
        auto [emin, emax] = iptr->get_minmax_ene(t1, t2, e1);

        //--------------------------------------------------
        // double counting prevention since we only consider targets with energy more than incident particle

        // TODO
        // require e1 < e2 = e_tar
        // e_tar = [emin, emax]
        // therefore, if e1 > emin we need to set emin' = max(e1, emin)
        // also, if e1 > emax -> interaction is not possible
        emin = std::max(e1, emin); 
          
        // ver2: other way around; require e1 > e2 = etarget
        // e_tar = [emin, emax]
        //
        // therefore emax' = min(e1, emax)
        // and if emin > e1 --> not possible
        //if(emin > e1){
        //  id++; // remembering to increase counter nevertheless
        //  continue; // no possible targets to interact with
        //}
        //emax = std::min(e1, emax);
        //--------------------------------------------------

        //if(emax > con_tar->eneArr[0]) jmin = 0;
        //if(emin > con_tar->eneArr[0]) jmax = 0;

        // NOTE: assuming reverse order here since ene is sorted in decreasing order
        // NOTE: assumes that sort_in_rev_energy is called in the container; this creates eneArr
        //int jmin = toolbox::find_rev_sorted_nearest( con_tar->eneArr, emax );
        //int jmax = toolbox::find_rev_sorted_nearest( con_tar->eneArr, emin );
        //int jme  = toolbox::find_rev_sorted_nearest( con_tar->eneArr, e1 );

        //size_t jmin = toolbox::find_rev_sorted_nearest(    con_tar->eneArr, emax );
        //size_t jmax = toolbox::revfind_rev_sorted_nearest( con_tar->eneArr, emin );

        // profiled to be fastest
        size_t jmin = toolbox::find_rev_sorted_nearest_algo2( con_tar->eneArr, emax );
        size_t jmax = toolbox::find_rev_sorted_nearest_algo2( con_tar->eneArr, emin );

        //std::cout << " efirst last: " << con_tar->eneArr[0] << " " << con_tar->eneArr[N2] << std::endl;
        //std::cout << "jminjmax " << jmin << " " << jmin2 << " _ " << jmax << " " << jmax2 << " s " << con_tar->eneArr.size() << " e " << emax << " " << emin << std::endl;

        //std::cout << "jmin " << jmin << " " << jmin2 << " s " << con_tar->eneArr.size() << " e " << emax << std::endl;
        //std::cout << "jmax " << jmax << " " << jmax2 << " s " << con_tar->eneArr.size() << " e " << emin << std::endl;
        //assert(jmin == static_cast<int>(jmin2));
        //assert(jmax == static_cast<int>(jmax2));


        // TEST remove energy check
        //jmin = 0;
        //jmax = N2;
                
        // ver1
        // total weight of prtcls between jmin-jmax
        // ver1: calculates total weight for every particle; value evolves dynamically
        //float wsum = 0.0f;
        //for(size_t j=jmin; j<jmax; j++) wsum += con_tar->wgt( j ); // sum( w[jmin:jmax] )

        // ver2
        // NOTE switching to this version as it is faster
        // assumes static targets; value calculated in the beginnig of loop
        float wsum2 = 0.0f;
          
        //--------------------------------------------------
        if(! iptr->do_accumulate ){ // normal mode; no accumulation

          if(jmin < jmax) { // in the opposite case arrays dont span a range and so wsum2 = 0
            float wsum_min = jmin == 0 ? 0.0f : con_tar->wgtCumArr[jmin];
            wsum2 = con_tar->wgtCumArr[jmax-1] - wsum_min;
          }

        //--------------------------------------------------
        } else { // accumulate interactions; effectively reduces weight
          float wprev = jmin == 0 ? 0.0 : con_tar->wgtCumArr[jmin-1];
          for(size_t j=jmin; j<jmax; j++) {
            
            // weight between [j-1, j]
            //w = con_tar->wgt(j); // real weight
            float w = con_tar->wgtCumArr[j] - wprev; // calc via cum array 
            wprev = con_tar->wgtCumArr[j];

            float e2 = con_tar->eneArr[j];
            auto [f1,f2] = iptr->accumulate(t1, e1, t2, e2);
            float f = f1*f2; //std::max(f1,f2);  // TODO is max ok here? or product?

            wsum2 += w/f; 
          }
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
        //    // or max // ver2
        //    //facc = std::max(facc, f);                        

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
        //// debug; thest how much wsum and wsum2 deviate
        ////std::cout << "wsum " <<  wsum << " " << wsum2 << " jminmax " << jmin << " " << jmax << std::endl;
        ////assert( std::abs( wsum - wsum2 ) < EPS );
           
        //--------------------------------------------------
        //total sum of target weights; TODO can this be removed?
        //float wtot = 0.0f;
        //for(size_t j=0; j<N2; j++) wtot += con_tar->wgt( j ); // sum( w )
        //if(wtot <= 0.0f) wtot = 1.0f; // guard for NaNs

        //--------------------------------------------------
        // maximum partial interaction rate
        float par_int_rate = 2.0f*cross_max * wsum2; // factor 2 comes from v_rel = 2c

        // update/accept interaction only if it has non-zero prob to occur
        if(par_int_rate > 0.0f) {

          probs.push_back(par_int_rate);
          cmaxs.push_back(cross_max);
          wsums.push_back(wsum2);
          jmins.push_back(jmin);
          jmaxs.push_back(jmax);
          ids  .push_back(id); // update interaction numbering

          // add stuff 
          //probs[id] = par_int_rate;
          //cmaxs[id] = cross_max;
          //wsums[id] = wsum2;
          //jmins[id] = jmin;
          //jmaxs[id] = jmax;
          //ids[id]   = id; // update interaction numbering

          //faccs.push_back(facc);
          //
          // debug prints
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
    const size_t N = probs.size();
    //double ptot = toolbox::sum(probs); // TODO
    float ptot = 0.0;
    for(size_t i=0; i<N; i++) ptot += probs[i];

    // calculate cumulative sum of all interactions
    ManVec<float> cumsum_int_probs;
    cumsum_int_probs.resize(N);

    cumsum_int_probs[0] = probs[0];
    for(size_t i=1; i<N; i++) cumsum_int_probs[i] = cumsum_int_probs[i-1] + probs[i];

    // normalize
    for(size_t i=0; i<N; i++) cumsum_int_probs[i] = cumsum_int_probs[i]/ptot;

    // draw random interaction from cumulative distribution
    int reti = toolbox::sample_prob(cumsum_int_probs, rand() );

    // get corresponding id
    return reti;
  }


  // energy-dependent weight adapation functions
  //
  // value corresponds to number of particles (of type t) produced with given energy x
  // e.g., constant value means that every interactions splits the particle into that many pieces
  float ene_weight_funs(std::string t, float x) 
  {

    // standard unit weights
    //if(       t == "ph") { return 1.0; 
    //} else if(t == "e-") { return 1.0; 
    //} else if(t == "e+") { return 1.0; 
    //}
      
    //// photon emphasis
    //if(       t == "ph") { return std::pow(x/0.01f, 0.1f); 
    //} else if(t == "e-") { return 1.0; 
    //} else if(t == "e+") { return 1.0; 
    //}

    //// photon emphasis
    if(       t == "ph") { return std::pow(x/0.01f, 0.04f); 
    } else if(t == "e-") { return 1.0; 
    } else if(t == "e+") { return 1.0; 
    }

    assert(false);
  }


  //--------------------------------------------------
  void solve_twobody(pic::Tile<D>& tile)
  {

    timer.start(); // start profiling block


    // build pointer map of types to containers; used as a helper to access particle types
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //--------------------------------------------------
    // call pre-iteration functions to update internal arrays 
    //for(auto&& con : tile.containers) {
    //  con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
    //}

    // keep this ordering; initialization of arrays assumes this way of calling the functions
    // NOTE: cannot move this inside the loop because particle removal assumes that indices remain static
    timer.start_comp("sort_ene");
    for(auto&& con : tile.containers) con.sort_in_rev_energy();
    timer.stop_comp("sort_ene");

    timer.start_comp("upd_cum_arr");
    for(auto&& con : tile.containers) con.update_cumulative_arrays();
    timer.stop_comp("upd_cum_arr");

    //--------------------------------------------------
    // collect statistics for bookkeeping
    std::map<std::string, int> info_prtcl_num;
    for(auto&& con : tile.containers) {
      auto t1 = con.type;
      info_prtcl_num[t1] = con.size();
    }

    const auto mins = tile.mins;
    const auto maxs = tile.maxs;

    //--------------------------------------------------

    // initialize temp variable storages
    float wmin; //, wmax, prob;
    float prob_upd3, prob_upd4, prob_kill3, prob_kill4;

    //--------------------------------------------------
    // loop over incident types

    // ver2: random shuffled indices of containers; makes iteration order random in every step
    //std::vector<std::string> t_inis = { "e-", "e+", "ph" };
    //std::shuffle(std::begin(t_inis), std::end(t_inis), gen);
    //for(auto t1 : t_inis)
    //{
    //  if(is_empty(t1)) continue; // no interactions with incident type t1
    //  auto con1 = cons[t1];


    // ver1: ordered iteration over prtcls
    for(auto&& con1 : tile.containers) 
    {
      auto t1 = con1.type;
      if(is_empty_binary_int(t1)) continue; // no interactions with incident type t1

      //size_t Ntot1 = con1.size();
      size_t Ntot1 = info_prtcl_num[t1]; // read particle number from here; 
                                         // it changes adaptively and routines assume non-evolving arrays

      // ver2
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
      //for(int n1=con1.size()-1; n1>=0; n1--) { // reverse iteration

        //unpack incident 
        auto lx1 = con1.loc(0,n1);
        auto ly1 = con1.loc(1,n1);
        auto lz1 = con1.loc(2,n1);

        auto ux1 = con1.vel(0,n1);
        auto uy1 = con1.vel(1,n1);
        auto uz1 = con1.vel(2,n1);
        auto w1  = con1.wgt(n1);

        auto e1 = con1.get_prtcl_ene(n1);

        if(w1 < EPS) continue; // omit zero-w incidents

        //pre-calculate maximum partial interaction rates
        timer.start_comp("comp_pmax");
        comp_pmax(t1, e1, cons); 
        timer.stop_comp("comp_pmax");

        if(ids.size() == 0) continue; // no targets to interact with 

        // total probability
        // NOTE: maximum interaction rate = sigma_max * w2_sum * w1/prob_norm
        float prob_vir_max = 0.0;
        for(size_t i=0; i<ids.size(); i++) prob_vir_max += 2.0f*cmaxs[i]*wsums[i]; 


        //--------------------------------------------------
        // no interactions allowed outside the gap
        // v0; sharp cutoff
        //if((std::abs(lx1/r_gap) > 1.0)) continue;

        // v1: smoothed drop
        // NOTE: currently this suppresses all interactions. 
        //       Should it only suppress 1) two-photon pair creation or 2) Compton?
        //       Or 3) reduce photon target densities for whatever reaction?
        // REASONING: In theory, yes but, in practice, both Compton and pair-creation 
        //       are suppressed so fast beyond x > H_gap that shouldn't matter. 
        // CAVEAT: This would matter, however, for pair annihilation that should always be
        //       possible. It's rate is however tiny (we hope).
        prob_vir_max *= shape(lx1, r_gap, 20.0f); // tanh profile with delta = 5 cells
        //--------------------------------------------------

        // no targets to interact with
        if(prob_vir_max < EPS) continue;

        //if(prob_vir_max>=1.01){ // some tolerance here
        //  std::cout << " prob_vir_max:" << prob_vir_max << std::endl;
        //  std::cout << " norm:" << prob_norm << std::endl;
        //  for(size_t i=0; i<ids.size(); i++) std::cout << "c: " << cmaxs[i] << " w: " << wsums[i] << " f: " << faccs[i] << std::endl;
        //  //assert(false);
        //}

        // exponential waiting time between interactions
        float t_free = -log( rand() )*prob_norm/(prob_vir_max*w1); //NOTE w1 here
                                                                      
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
          timer.start_comp("draw_proc");
          int i = draw_rand_proc(); // i:th interaction in probs array
          timer.stop_comp("draw_proc");

          // NOTE i = ids[i]; we do not factor this out because in theory ids array could change in size and then this is needed.
          
          //--------------------------------------------------
          // unpack interaction
          auto jmin = jmins[i];
          auto jmax = jmaxs[i];
          auto cmax = cmaxs[i];
          //auto wsum = wsums[i];
          //auto facc_max = faccs[i];

          auto int_id = ids[i];                    // id in global array
          auto iptr = binary_interactions[int_id]; // pointer to interaction 
          auto t2   = iptr->t2;                    // target type
          auto con2 = cons[t2];                    // target container

          if(con2->size()==0) continue;

          //--------------------------------------------------
          // get random target with energy between jmin/jmax
          // propability is proportional to weight of LPs

          timer.start_comp("sample_prob");
          //size_t n2 = toolbox::sample_prob_between(     con2->wgtCumArr, rand(), jmin, jmax);
          size_t n2 = toolbox::sample_prob_between_algo(con2->wgtCumArr, rand(), jmin, jmax);
          timer.stop_comp("sample_prob");

          //std::cout << "n2/3" << n2 << " " << n3 << std::endl;
          //assert(n2 == n3);

          if( (t1 == t2) && (n1 == n2) ) continue; // do not interact with self

          // unpack target
          auto lx2 = con2->loc(0,n2);
          auto ly2 = con2->loc(1,n2);
          auto lz2 = con2->loc(2,n2);

          auto ux2 = con2->vel(0,n2);
          auto uy2 = con2->vel(1,n2);
          auto uz2 = con2->vel(2,n2);
          auto w2  = con2->wgt(  n2);
          auto e2  = con2->get_prtcl_ene(n2);

          if(w2 < EPS) continue; // omit zero-w targets
                                   
          //--------------------------------------------------
          // double counting prevention by considering only e1 < e2 cases
          //if(e1 > e2) continue;
                  
          //avoid double counting by considering only e1 > e2 cases
          //if(e1 < e2) continue;
          //could also put double counting check to comp_pmax
          //--------------------------------------------------
            
          auto lx3 = lx1; //(lx1 + lx2)*0.5; 
          auto ly3 = ly1; //(ly1 + ly2)*0.5; 
          auto lz3 = lz1; //(lz1 + lz2)*0.5; 

          auto lx4 = lx2; //(lx1 + lx2)*0.5; 
          auto ly4 = ly2; //(ly1 + ly2)*0.5; 
          auto lz4 = lz2; //(lz1 + lz2)*0.5; 
                            
#ifdef DEBUG
          //--------------------------------------------------
          // confirm that locations are inside tile

          //size_t loc_flag1 = 0;
          //if(D>= 1 && mins[0]-2 <= lx1 && lx1 <= maxs[0]+2) loc_flag1++;
          //if(D>= 2 && mins[1]-2 <= ly1 && ly1 <= maxs[1]+2) loc_flag1++;
          //if(D>= 3 && mins[2]-2 <= lz1 && lz1 <= maxs[2]+2) loc_flag1++;
          //                         
          //size_t loc_flag2 = 0;
          //if(D>= 1 && mins[0]-2 <= lx2 && lx2 <= maxs[0]+2) loc_flag2++;
          //if(D>= 2 && mins[1]-2 <= ly2 && ly2 <= maxs[1]+2) loc_flag2++;
          //if(D>= 3 && mins[2]-2 <= lz2 && lz2 <= maxs[2]+2) loc_flag2++;

          size_t loc_flag3 = 0;
          if(D>= 1 && mins[0]-3 <= lx3 && lx3 <= maxs[0]+2) loc_flag3++;
          if(D>= 2 && mins[1]-3 <= ly3 && ly3 <= maxs[1]+2) loc_flag3++;
          if(D>= 3 && mins[2]-3 <= lz3 && lz3 <= maxs[2]+2) loc_flag3++;


          size_t loc_flag4 = 0;
          if(D>= 1 && mins[0]-3 <= lx4 && lx4 <= maxs[0]+2) loc_flag4++;
          if(D>= 2 && mins[1]-3 <= ly4 && ly4 <= maxs[1]+2) loc_flag4++;
          if(D>= 3 && mins[2]-3 <= lz4 && lz4 <= maxs[2]+2) loc_flag4++;

          //if(loc_flag1 < D ) {
          //    std::cerr << "PAIRING: loc err 1"
          //      << " x=" << lx1
          //      << " y=" << ly1
          //      << " z=" << lz1
          //      << " ux=" << ux1
          //      << " uy=" << uy1
          //      << " uz=" << uz1
          //      << " w="  << w1
          //      << " minsx:" << mins[0] 
          //      << " minsy:" << mins[1] 
          //      << " minsz:" << mins[2] 
          //      << " maxsx:" << maxs[0] 
          //      << " maxsy:" << maxs[1] 
          //      << " maxsz:" << maxs[2] 
          //      << std::endl;
          //    assert(false);
          //}

          //if(loc_flag2 < D ) {
          //    std::cerr << "PAIRING: loc err 2"
          //      << " x=" << lx2
          //      << " y=" << ly2
          //      << " z=" << lz2
          //      << " ux=" << ux2
          //      << " uy=" << uy2
          //      << " uz=" << uz2
          //      << " w="  << w2
          //      << " minsx:" << mins[0] 
          //      << " minsy:" << mins[1] 
          //      << " minsz:" << mins[2] 
          //      << " maxsx:" << maxs[0] 
          //      << " maxsy:" << maxs[1] 
          //      << " maxsz:" << maxs[2] 
          //      << std::endl;
          //    assert(false);
          //}

          if(loc_flag3 < D ) {
              std::cerr << "PAIRING: loc err 3"
                << " t="  << t1
                << " n="  << n1
                << " x=" << lx3
                << " y=" << ly3
                << " z=" << lz3
                << " ux=" << ux1
                << " uy=" << uy1
                << " uz=" << uz1
                << " w="  << w1;

              if(D >= 1) std::cerr << " minsx:" << mins[0];
              if(D >= 2) std::cerr << " minsy:" << mins[1];
              if(D >= 3) std::cerr << " minsz:" << mins[2];
              if(D >= 1) std::cerr << " maxsx:" << maxs[0];
              if(D >= 2) std::cerr << " maxsy:" << maxs[1];
              if(D >= 3) std::cerr << " maxsz:" << maxs[2];
              std::cerr << std::endl;

              assert(false);
          }

          if(loc_flag4 < D ) {
              std::cerr << "PAIRING: loc err 4"
                << " t="  << t2
                << " n="  << n2
                << " x=" << lx4
                << " y=" << ly4
                << " z=" << lz4
                << " ux=" << ux2
                << " uy=" << uy2
                << " uz=" << uz2
                << " w="  << w2;

              if(D >= 1) std::cerr << " minsx:" << mins[0];
              if(D >= 2) std::cerr << " minsy:" << mins[1];
              if(D >= 3) std::cerr << " minsz:" << mins[2];
              if(D >= 1) std::cerr << " maxsx:" << maxs[0];
              if(D >= 2) std::cerr << " maxsy:" << maxs[1];
              if(D >= 3) std::cerr << " maxsz:" << maxs[2];
              std::cerr << std::endl;

              assert(false);
          }
#endif
          
          //--------------------------------------------------
          
          // real probablity of interaction
          timer.start_comp("comp_cs");
          auto [cm, vrel] = iptr->comp_cross_section(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2 );
          timer.stop_comp("comp_cs");

          // collect max cross section
          float cm_cur = info_max_int_cs[iptr->name];
          info_max_int_cs[iptr->name] = std::max( cm_cur, cm*vrel/2.0f);

          // comparison of interaction to max interaction 
          float prob_vir = cm*vrel/(2.0*cmax);

          // correct average accumulation factor with the real value
          timer.start_comp("acc");
          auto [facc3, facc4] = iptr->accumulate(t1, e1, t2, e2);
          facc3 = iptr->do_accumulate ? facc3 : 1.0f;
          facc4 = iptr->do_accumulate ? facc4 : 1.0f;
          timer.stop_comp("acc");


          // FIXME remove check if sure this works
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


          if(rand() < prob_vir)  // check if this interaction is chosen among the sampled ones
          {
            auto long_name = iptr->name + "_" + t1 + "_" + t2;
            info_int_nums[long_name] += 1;

            // particle values after interaction
            timer.start_comp("dupl_prtcl");
            auto [t3, ux3, uy3, uz3, w3] = duplicate_prtcl(t1, ux1, uy1, uz1, w1);
            auto [t4, ux4, uy4, uz4, w4] = duplicate_prtcl(t2, ux2, uy2, uz2, w2);
            timer.stop_comp("dupl_prtcl");

            // interact and udpate variables in-place
            timer.start_comp("interact");
            iptr->interact( t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4 );
            timer.stop_comp("interact");

            // new energies; NOTE: could use container.m to get the mass
            float m3 = (t3 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
            float m4 = (t4 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
            float e3 = std::sqrt( m3*m3 + ux3*ux3 + uy3*uy3 + uz3*uz3 );
            float e4 = std::sqrt( m4*m4 + ux4*ux4 + uy4*uy4 + uz4*uz4 );

            // weight adaptation; original Stern95 version
            //float fw3 = (e3/e1)*ene_weight_funs(t1, e1)/ene_weight_funs(t3, e3); 
            //float fw4 = (e4/e2)*ene_weight_funs(t2, e2)/ene_weight_funs(t4, e4); 

            // more intuitive version (flipped)
            timer.start_comp("weight_funs");
            float fw3 = ene_weight_funs(t3, e3)/ene_weight_funs(t1, e1);
            float fw4 = ene_weight_funs(t4, e4)/ene_weight_funs(t2, e2);
            timer.stop_comp("weight_funs");

            // limit explosive particle creation
            float n3 = std::min( fw3, 32.0f );
            float n4 = std::min( fw4, 32.0f );

            //--------------------------------------------------
            if(t1 == t3 && t2 == t4) { // scattering interactions

              // redistribute weights among the new copies
              w3 = w1/n3; 
              w4 = w2/n4; 
                
              wmin = min(w3, w4); // minimum available weight
              //wmax = max(w3, w4); // maximum available weight

              //# NOTE these two expressions are equal: w_i = w_j/wmax == wmin/w_i
              //# this is where the algorithm differs from original LP MC method by Stern95;
              //# here, instead of killing the LP we do not update its energy.

              // we can think here that only wmin = min(w1, w2) fraction of the LPs weight is 
              // used in the interaction. If everything is used then (i.e., wmin=w1 for LP1 
              // or wmin=w2 for LP2)) then prob_upd = 1

              // NOTE: update probability is later reduced as if particle weights are w -> w*f_acc 
              //       when accumulation is used

              // NOTE every new particle has a change of prob_upd of being updated so this rejection sample 
              // needs to be calculated using the new weights; then multiple evaluations of the rejection sampling
                

              // probability for the particle energy to be updated 
              // NOTE: equals the prob_upd calculated with parent weights. Alternatively, we could check prob_upd
              // once and add all the copies thereafter; this however leads to more MC sampling noise.
              prob_upd3 = wmin/w3; //w2/wmax;
              prob_upd4 = wmin/w4; //w1/wmax;

              // check against wrong usage of these
              prob_kill3 = -1.0f;
              prob_kill4 = -1.0f;
            } else { // annihilation interactions

              // weight of the new type is the minimum of the two incident particles
              // share mutual weight between outgoing prtcls
              wmin = min(w1, w2); // minimum available weight
              w3 = wmin/n3; 
              w4 = wmin/n4; 

              // probability to kill the parent (since not all of it may not be used)
              prob_kill3 = wmin/w1;
              prob_kill4 = wmin/w2;
                
              // check against wrong usage of these
              prob_upd3 = -1.0f;
              prob_upd4 = -1.0f;
            }


            //--------------------------------------------------
            // branch for rescaling electron-positron plasma to have unitary weights

            // t1
            if(t1 == t3 && t2 == t4) { // scattering interactions

              // t3
              if(force_ep_uni_w && (t3 == "e-" || t3 == "e+") ){ 
                w3 = 1.0f;
                n3 = facc3*w1/w3; // remembering to increase prtcl num w/ facc
                facc3 = 1.0f; // restore facc (since it is taken care of by n3)
              }
                
              //4
              if(force_ep_uni_w && (t4 == "e-" || t4 == "e+") ){
                w4 = 1.0f;
                n4 = facc4*w2/w4; // remembering to increase prtcl num w/ facc
                facc4 = 1.0f; // restore facc
              }

              // remember to recalc prob_upd
              wmin = min(w3, w4); // minimum available weight
              prob_upd3 = wmin/w3; 
              prob_upd4 = wmin/w4; 

            } else { // annihilation interactions
                       
              // t3
              if(force_ep_uni_w && (t3 == "e-" || t3 == "e+") ){
                w3 = 1.0;
                n3 = wmin/w3;
              }

              // t4
              if(force_ep_uni_w && ( t4 == "e-" || t4 == "e+") ){
                w4 = 1.0;
                n4 = wmin/w4;
              }
            }


            //--------------------------------------------------
            // split too big LPs
            //if(w3 > 1.0e5f) {
            //  w3 *= 0.5f;
            //  n3 *= 2.0f;
            //}

            //if(w4 > 1.0e5f) {
            //  w4 *= 0.5f;
            //  n4 *= 2.0f;
            //}
            //--------------------------------------------------


            //--------------------------------------------------
            // accumulation means that we treat the particles weight as w -> w*f_cc
            // but consider the interaction less often as prob -> prob/f_acc
            // to compensate, this means that we need to also genereate n_cop -> ncop/f_acc 
            // less of particles.
            //
            // NOTE: this increases the weight of the new LP3 automatically, i.e., we do not need 
            // to multiply w -> w*f_acc anymore.


            //if( w3*facc3 > 1.0e8f || w4*facc4 > 1.0e8f) {
            //  std::cout << " ERR: pairing" << std::endl;
            //  std::cout << " w3: " << w3 << std::endl;
            //  std::cout << " w4: " << w4 << std::endl;
            //  std::cout << " e3: " << e3 << std::endl;
            //  std::cout << " e4: " << e4 << std::endl;
            //  std::cout << " n3: " << n3 << std::endl;
            //  std::cout << " n4: " << n4 << std::endl;
            //  std::cout << " f3: " << facc3 << std::endl;
            //  std::cout << " f4: " << facc4 << std::endl;
            //  assert(false);
            //}

            //std::cout << " ERR: pairing" << std::endl;
            //std::cout << " ts: " << t1 << " " << t2 << " " << t3 << " " << t4 << std::endl;
            //std::cout << " w4: " << w4 << std::endl;
            //std::cout << " e3: " << e3 << std::endl;
            //std::cout << " e4: " << e4 << std::endl;
            //std::cout << " w3: " << w3 << std::endl;
            //std::cout << " w4: " << w4 << std::endl;
            //std::cout << " e3: " << e3 << std::endl;
            //std::cout << " e4: " << e4 << std::endl;
            //std::cout << " n3: " << n3 << std::endl;
            //std::cout << " n4: " << n4 << std::endl;
            //std::cout << " f3: " << facc3 << std::endl;
            //std::cout << " f4: " << facc4 << std::endl;
            //std::cout << " prob3:  " << prob_upd3 << std::endl;
            //std::cout << " prob4:  " << prob_upd4 << std::endl;
            //std::cout << " prob_k3:" << prob_kill3 << std::endl;
            //std::cout << " prob_k4:" << prob_kill4 << std::endl;


            //-------------------------------------------------- 
            double ncop = 0.0;

            //if(t1 == t3) { // same type before/after interactions; update energy with prob_upd
            if(t1 == t3 && t2 == t4) { // scattering interactions
                             
              // -------------scattering interactions go here----------------
                
              // guard against wrong branching; just a double check, can be removed
              assert(prob_upd3 >= 0.0f);
              assert(prob_upd4 >= 0.0f);


              timer.start_comp("add_sc_prtcl1");
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
                    cons[t1]->wgt(   n1) = w3*facc3;
                  } else {
                    cons[t1]->wgt(n1) = w3;
                  }
                } else {
                  if( rand() < prob_upd3 ) {
                    cons[t1]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w3*facc3); // new ene & w
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
              timer.stop_comp("add_sc_prtcl1");

              timer.start_comp("del_parent1");
              // remove parent prtcl if nothing was added
              if( ncop < EPS ) {
                //cons[t1]->to_other_tiles.push_back( {1,1,1,n1} ); // NOTE: CPU version
                cons[t1]->info(n1) = -1; // mark for deletion via outflow routines
                cons[t1]->wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              }
              timer.stop_comp("del_parent1");

            //--------------------------------------------------
            } else { //# different before/after type; kill parent with a prob_upd

              // annihilation interactions go here
                
              // guard against wrong branching; just a double check, can be removed
              assert(prob_kill3 >= 0.0f);
              assert(prob_kill4 >= 0.0f);

              // TODO are these independent or same draw for prob_kill3
              // i.e., kill parent and create copies or let parent live and no copies?

              bool do_addition = true;
              if(  (t3 == "ph")                  && (info_prtcl_num[t3] > max_tile_phot_num ) ) do_addition = false; // switch off for photons
              if( ((t3 == "e-") || (t3 == "e+")) && (info_prtcl_num[t3] > max_tile_prtcl_num) ) do_addition = false; // switch off for pairs

              if( do_addition ) { // add if we are below tile limit
                                                   
                timer.start_comp("add_prtcl1");
                double z1 = rand();
                while( n3 > z1 + ncop ){
                  // TODO NOTE lx3 here not lx1
                  cons[t3]->add_particle( {{lx3, ly3, lz3}}, {{ux3, uy3, uz3}}, w3); // new ene & w
                  ncop += 1.0;
                }
                timer.stop_comp("add_prtcl1");
              }

              timer.start_comp("del_parent1");
              // kill parent
              if( prob_kill3 > rand() ) {
                //cons[t1]->to_other_tiles.push_back( {1,1,1,n1} ); // NOTE: CPU version
                cons[t1]->info(n1) = -1; // mark for deletion via outflow routines
                cons[t1]->wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              }
              timer.stop_comp("del_parent1");

            } // end of prtcl t1/t3 addition


            //-------------------------------------------------- 
            // add prtcl t2/t4
            ncop = 0.0;

            //if(t2 == t4) { // same type before/after interactions; update energy with prob_upd
            if(t1 == t3 && t2 == t4) { // scattering interactions

              // scattering interactions go here

              timer.start_comp("add_sc_prtcl2");
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
                    cons[t2]->wgt(   n2) = w4*facc4;
                  } else {
                    cons[t2]->wgt(n2) = w4;
                  }
                } else {
                  if( rand() < prob_upd4 ) {
                    cons[t2]->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w4*facc4); // new ene & w
                  } else {
                    cons[t2]->add_particle( {{lx2, ly2, lz2}}, {{ux2, uy2, uz2}}, w4); // new w
                  }
                }

                ncop += 1.0;
              } // end of while
              timer.stop_comp("add_sc_prtcl2");

              timer.start_comp("del_parent2");
              // remove parent prtcl if nothing was added
              if( ncop < EPS ) {
                //cons[t2]->to_other_tiles.push_back( {1,1,1,n2} ); // NOTE: CPU version
                cons[t2]->info(n2) = -1; // mark for deletion via outflow routines
                cons[t2]->wgt(n2) = 0.0f; // make zero wgt so its omitted from loop
              }
              timer.stop_comp("del_parent2");

            //--------------------------------------------------
            } else { //# different before/after type; kill parent with a prob_upd
                       
              // annihilation interactions go her
                
              bool do_addition = true;
              if(  (t4 == "ph")                  && (info_prtcl_num[t4] > max_tile_phot_num ) ) do_addition = false; // switch off for photons
              if( ((t4 == "e-") || (t4 == "e+")) && (info_prtcl_num[t4] > max_tile_prtcl_num) ) do_addition = false; // switch off for pairs

              if( do_addition ) { // add if we are below tile limit
                                                 //
                timer.start_comp("add_prtcl2");
                double z1 = rand();
                while( n4 > z1 + ncop ){
                  //cons[t4]->add_particle( {{lx2, ly2, lz2}}, {{ux4, uy4, uz4}}, w4); // new ene & w
                  // TODO note lx4 here
                  cons[t4]->add_particle( {{lx4, ly4, lz4}}, {{ux4, uy4, uz4}}, w4); // new ene & w
                  ncop += 1.0;
                }
                timer.stop_comp("add_prtcl2");
              }

              timer.start_comp("del_parent2");
              // kill parent
              if( prob_kill4 > rand() ) {
                //cons[t2]->to_other_tiles.push_back( {1,1,1,n2} ); // NOTE: CPU version
                cons[t2]->info(n2) = -1; // mark for deletion
                cons[t2]->wgt(n2) = 0.0f; // make zero wgt so its omitted from loop
              }
              timer.stop_comp("del_parent2");

            } // end of prtcl t1/t3 addition

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

      int num_of_dels = 0;
      #pragma omp simd reduction(+:num_of_dels)
      for(int n=0; n<cons[t1]->size(); n++){
        num_of_dels += cons[t1]->info(n) == -1 ? 1 : 0; // if -1 add one
      }
      info_prtcl_kill[t1] = num_of_dels;
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
    timer.start_comp("del");
    for(auto&& con : tile.containers)
    {
      con.delete_transferred_particles(); // remove annihilated prtcls; this transfer storage 
                                          // is used as a tmp container for storing the indices
    }
    timer.stop_comp("del");


    timer.stop();
    return;
  }


  //--------------------------------------------------
  // one-body single particle interactions
  void solve_onebody(pic::Tile<D>& tile)
  {
    timer.start(); // start profiling block

    // build pointer map of types to containers; used as a helper to access particle types
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //--------------------------------------------------
    // call pre-iteration functions to update internal arrays 
    //for(auto&& con : tile.containers) {
    //  con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
    //}


    //--------------------------------------------------
    // collect statistics for bookkeeping
    std::map<std::string, int> info_prtcl_num;
    for(auto&& con : tile.containers) {
      auto t1 = con.type;
      info_prtcl_num[t1] = con.size();
    }

    const auto mins = tile.mins;
    //const auto maxs = tile.maxs;

    const auto& gs = tile.get_grids(); 

    //--------------------------------------------------
    // initialize temp variable storages
    float ux4, uy4, uz4, w4=0.0; // empty particle created in p1 -> p3 + p4 splitting
    std::string t4;  // type variable for secondary prtcl;


    // ver1: ordered iteration over prtcls
    for(auto&& con1 : tile.containers) 
    {
      auto t1 = con1.type;

      if(is_empty_single_int(t1)) continue; // no interactions with incident type t1


      //--------------------------------------------------
      // select interaction
      // TODO assume that each particle type has only one interaction
      //      in a more general case, we could use the same virtual channel as in binary interactions

      size_t id = 0;
      for(auto iptr : single_interactions)
      {
        if(t1 == iptr->t1) break; // assume only one target
        id += 1;
      }
      auto iptr = single_interactions[id]; // interaction


      //--------------------------------------------------
      // loop over 
      size_t Ntot1 = info_prtcl_num[t1]; // read particle number 
      for(size_t n1=0; n1<Ntot1; n1++) {

        //unpack incident 
        const auto lx1 = con1.loc(0,n1);
        const auto ly1 = con1.loc(1,n1);
        const auto lz1 = con1.loc(2,n1);

        const auto xborn = ly1; //Used only in 1D

        const auto ux1 = con1.vel(0,n1);
        const auto uy1 = con1.vel(1,n1);
        const auto uz1 = con1.vel(2,n1);

        const auto w1  = con1.wgt(n1);
        const auto e1  = con1.get_prtcl_ene(n1);

        if(w1 < EPS) continue; // omit zero-w incidents

        const auto [emin, emax] = iptr->get_minmax_ene(t1, "", e1);

        //std::cout << "emin/emax " << e1 << " " << emin << " " << emax << "\n";

        if( e1 < emin ) continue; // low-energy cutoff
        if( e1 > emax ) continue; // high-energy cutoff

        //--------------------------------------------------
        // v1; active interpolation (via nearest neighbor)
          
        // get E and B field values
        int i=0, j=0, k=0;
        if(D >= 1) i = static_cast<int>(floor(lx1) - mins[0]);
        if(D >= 2) j = static_cast<int>(floor(ly1) - mins[1]);
        if(D >= 3) k = static_cast<int>(floor(lz1) - mins[2]);

        //--------------------------------------------------
        const size_t ind = gs.ex.indx(i,j,k);

        //--------------------------------------------------
        // construct (optional) virtual curvature into the EM fields (for 1D cases)
        // We mimic a real B-field line curvature by inserting the equivalent field perpendicular component to B_y via by_vir parameter.

        // Synchrotron: particle experiences an angle sin\theta = \gamma r_g/R_curv where r_g is the gyroradius and R_curv the field line curvature
        // in its rest frame (hence an additional \gamma factor).

        // Multi-photon annihilation: Photons are emitted mostly to the direction of the electron and as they propagate, have an increasing angle against the curving field line.
        // To mimic  this, we calculate the angle via \sin\theta = sin((x-x_born)/R_curv).
        // The height where the photon is born (x_born) is stored in 1D simulations within the y-location of the particle.

        float by_vir = 0.0f;
        if (use_vir_curvature) {
          if(iptr->name == "multi-phot-ann"){
            by_vir = gs.bx(ind)*std::sin( (lx1 - xborn)/r_curv ); // v0 

            // more complicated v1 that should be numerically more stable and has radius dependency
            //by_vir = gs.bx(ind)*std::abs(lx1 - xborn)/r_curv;        // approximate sin\theta \approx \theta
            //by_vir = by_vir*std::pow(1.0f - std::abs(lx1/r_gap), 2);  // decrease field strength linearly with height

          } else if(iptr->name == "synchrotron") {
            float gam = sqrt(1.0 + ux1*ux1 + uy1*uy1 + uz1*uz1 );
            by_vir = gs.bx(ind)*gam*vir_pitch_ang; // \gamma B_x \sin\alpha
          }

          // turn off pair production for h/L > 1
          // v0: sharp drop
          //if(std::abs(lx1/r_gap) > 1.0) by_vir = 0.0f;

          // v1: smoothed drop
          by_vir *= shape(lx1, r_gap, 20.0f); // tanh profile with delta = 5 cells
        }

        const float ex = gs.ex(ind); 
        const float ey = gs.ey(ind); 
        const float ez = gs.ez(ind); 

        const float bx = gs.bx(ind); 
        const float by = gs.by(ind) + by_vir; 
        const float bz = gs.bz(ind); 

        
        //--------------------------------------------------
        // v2; passive fetching; assumes a call has been made to interp before this function
        // NOTE does not work because this solve_onebody method modifies the arrays with add_prtcl1; 
        //      this causes epart and bpart arrays to not be in sync 

        //auto ex2 = con1.ex(n1); 
        //auto ey2 = con1.ey(n1); 
        //auto ez2 = con1.ez(n1); 

        //auto bx2 = con1.bx(n1); 
        //auto by2 = con1.by(n1); 
        //auto bz2 = con1.bz(n1); 

        //std::cout << " nnbor: " << ex << " " << ey << " " << ez << " " << bx << " " << by << " " << bz << "\n"; 
        //std::cout << " inter: " << ex2 << " " << ey2 << " " << ez2 << " " << bx2 << " " << by2 << " " << bz2 << "\n"; 


        // local optical depth; 
        // NOTE: em field is stored during this call and does not need to be called again in interact()
        timer.start_comp("optical_depth");
        const float tau_int = iptr->comp_optical_depth(
                                  t1, 
                                  ux1, uy1, uz1, 
                                  ex, ey, ez, 
                                  bx, by, bz);
        timer.stop_comp("optical_depth");

        // exponential waiting time between interactions
        const float t_free = -log( rand() )*prob_norm_onebody/tau_int; 

        //std::cout << "oneb: t_free:" << t_free << " tau: " << tau_int << " N:" << prob_norm_onebody << "\n";

        if(t_free < 1.0) // interact
        { 

          // particle values after interaction
          auto [t3, ux3, uy3, uz3, w3] = duplicate_prtcl(t1, ux1, uy1, uz1, w1);

          iptr->wtar2wini = prob_norm_onebody; // inject N1Q normalization into the interaction for possible extra calculations

          timer.start_comp("interact");
          iptr->interact( t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4);
          timer.stop_comp("interact");

          // new energies; NOTE: could use container.m to get the mass
          const float m3 = (t3 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
          const float m4 = (t4 == "ph") ? 0.0f : 1.0f; // particle mass; zero if photon
          const float e3 = std::sqrt( m3*m3 + ux3*ux3 + uy3*uy3 + uz3*uz3 );
          const float e4 = std::sqrt( m4*m4 + ux4*ux4 + uy4*uy4 + uz4*uz4 );

          timer.start_comp("weight_funs");
          // NOTE both are compared to the same parent t1 
          float fw3 = ene_weight_funs(t3, e3)/ene_weight_funs(t1, e1); // possible re-weighting of the parent particle 
          float fw4 = ene_weight_funs(t4, e4)/ene_weight_funs(t1, e1); // re-weighting of the secondary particle 
          timer.stop_comp("weight_funs");

          // limit explosive particle creation
          float n3 = std::min( fw3, 32.0f );
          float n4 = std::min( fw4, 32.0f );

          //--------------------------------------------------
          // accumulation
          // NOTE: accumulation works differently in onebody interactions than in binary interactions.
          //       Here, we reduce the number of particles created and increase the weight. In the binary interactions,
          //       the probability of the process is reduced and the process itself never occurs if accumulated.
          //       The way accumulation is done here is better for pruning the low-energy synchrotorn photons.

          timer.start_comp("acc");
          auto [facc3, facc4] = iptr->accumulate(t3, e3, t4, e4); // updated parent and emitted prtcl
          facc3 = iptr->do_accumulate ? facc3 : 1.0f;
          facc4 = iptr->do_accumulate ? facc4 : 1.0f;
          timer.stop_comp("acc");

          //n3 = n4/facc3; // NOTE never modify the parent; could be implemented but then need to change also the 
                           //      parent update below.

          n4 = n4/facc4; 

          // NOTE: weight is automatically modified correctly since it is updated based on number of copies

          //--------------------------------------------------
          if(t1 == t3){ // single-body emission 

            w3 = w1/n3; // new parent particle weight
            w4 = w1/n4; // emitted prtcl inherits weight from parent

            w4 = w4*iptr->wtar2wini; //add a growth factor calculated inside the processes' interact routine
                          
            // NOTE: we keep location the same
            con1.vel(0,n1) = ux3;
            con1.vel(1,n1) = uy3;
            con1.vel(2,n1) = uz3;
            // NOTE assume that weight w3 does not change; therefore, no need to add via MC  routine
            // TODO ignoring this weight change is correct only when onebody interactions include synch and multiphotann
            //      if new processes are added (that can have a parent particle that is not electron/positron) then this 
            //      needs re-updating.


            bool do_addition = true;
            if(  (t4 == "ph")                  && (info_prtcl_num[t4] > max_tile_phot_num)  ) do_addition = false; // switch off for photons
            if( ((t4 == "e-") || (t4 == "e+")) && (info_prtcl_num[t4] > max_tile_prtcl_num) ) do_addition = false; // switch off for pairs

            if( do_addition ) { // add if we are below tile limit
                                                                
              // add prtcl 4
              timer.start_comp("add_ems_prtcls");
              float ncop = 0.0;
              float z1 = rand();

              //saving the x-value of the created photon for 1D calculation if virtual curv used
              float ly1vir = use_vir_curvature ? lx1 : ly1;

              //if(use_vir_curvature && lx1 > r_gap) z1 = n4 + 0.5; // prevent emission beyond gap size

              while(n4 > z1 + ncop) {
                cons[t4]->add_particle( {{lx1, ly1vir, lz1}}, {{ux4, uy4, uz4}}, w4);
                ncop += 1.0;
              }
              timer.stop_comp("add_ems_prtcls");
            }

          //--------------------------------------------------
          } else { // single-body annihilation into t3 and t4 pair

              if(force_ep_uni_w && (t3 == "e-" || t3 == "e+") ){ 
                w3 = 1.0f;
                n3 = w1/w3; // remembering to increase prtcl num w/ facc
              }

              if(force_ep_uni_w && (t4 == "e-" || t4 == "e+") ){
                w4 = 1.0f;
                n4 = w1/w4; // NOTE w1 here since parent is same for both t3 and t4
              }

              bool do_addition = true;
              if( ((t3 == "e-") || (t3 == "e+")) && (info_prtcl_num[t3] > max_tile_prtcl_num) ) do_addition = false; // switch off for pairs
              if( ((t4 == "e-") || (t4 == "e+")) && (info_prtcl_num[t4] > max_tile_prtcl_num) ) do_addition = false; // switch off for pairs
              if(  (t3 == "ph")                  && (info_prtcl_num[t3] > max_tile_phot_num ) ) do_addition = false; // switch off for photons
              if(  (t4 == "ph")                  && (info_prtcl_num[t4] > max_tile_phot_num ) ) do_addition = false; // switch off for photons
                                                                                                                       
                                                                                                                       
              if( do_addition ) { // add if we are below tile limit
                
                // add new particle t3 and t4; particles are assumed to be identical
                timer.start_comp("add_ann_prtcls");
                float ncop = 0.0;
                float z1 = rand();
                while(n4 > z1 + ncop) {
                  cons[t3]->add_particle( {{lx1, ly1, lz1}}, {{ux3, uy3, uz3}}, w3); 
                  cons[t4]->add_particle( {{lx1, ly1, lz1}}, {{ux4, uy4, uz4}}, w4); 
                  ncop += 1.0;
                }
                timer.stop_comp("add_ann_prtcls");
              }

              // remove old parent particle t1
              timer.start_comp("del_parent");
              cons[t1]->info(n1) = -1; //to_other_tiles.push_back( {1,1,1,n1} ); // NOTE: CPU version
              cons[t1]->wgt(n1) = 0.0f; // make zero wgt so its omitted from loop
              timer.stop_comp("del_parent");
          }

        } // if interact
      } // end over n1 prtcl loop
    } // end of con1 loop


    //--------------------------------------------------
      
    timer.start_comp("del");
    for(auto&& con : tile.containers)
    {
      con.delete_transferred_particles(); // remove annihilated prtcls; this transfer storage 
                                          // is used as a tmp container for storing the indices
    }
    timer.stop_comp("del");


    timer.stop(); // start profiling block
  }




  //--------------------------------------------------
  // normalize container of type t1
  void rescale(pic::Tile<D>& tile, string& t1, double f_kill)
  {

    // build pointer map of types to containers; used as a helper to access particle types
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //cons[t1]->to_other_tiles.clear(); // clear book keeping array

    size_t N1 = cons[t1]->size(); // read particle number from here; 

    // total weight = sum(ws)
    float wtot = 0.0;
    for(size_t n1=0; n1<N1; n1++) wtot += cons[t1]->wgt(n1);

    // TODO it is not energy conserving to select particles equally w/o w-weighting

    // loop over particles
    float w1, zeta, prob_kill;
    for(size_t n1=0; n1<N1; n1++) {
      w1  = cons[t1]->wgt(n1);

      prob_kill = 1.0f - 1.0f/f_kill;

      zeta = rand();
      if( zeta < prob_kill) {
        cons[t1]->info(n1) = -1; //to_other_tiles.push_back( {1,1,1,n1} ); // NOTE: CPU deletion version
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
      float temp_inj, 
      float wph_inj,
      float Nph_inj) 
  {
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    auto mins = tile.mins;
    auto maxs = tile.maxs;

    float x_min = D >= 1 ? mins[0] : 0.0f;
    float y_min = D >= 2 ? mins[1] : 0.0f;
    float z_min = D >= 3 ? mins[2] : 0.0f;

    float lenx = D >= 1 ? (maxs[0] - mins[0]) : 0.0f;
    float leny = D >= 2 ? (maxs[1] - mins[1]) : 0.0f;
    float lenz = D >= 3 ? (maxs[2] - mins[2]) : 0.0f;

    // inject Nph_inj photons
    float vx, vy, vz, xi;
    float ux, uy, uz, xinj;
    float xloc, yloc, zloc;

    float ncop = 0.0f;
    float z1 = rand();

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
      float xi1 = rand();
      float xi2 = rand();
      float xi3 = rand();
      float xi4 = rand();
    
      float xi, jj, fsum;
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
      float slope, 
      float pmin,
      float pmax,
      float w_inj,
      float N_inj) 
  {

    //const float pmin = 10.0f;
    //const float pmax = 100.0f;

    assert(pmin > 1.0f); // pmin interpreted as gamma 
    assert(pmax > 1.0f); // pmax interpreted as gamma 

    float pminp = std::pow(pmin, 1.0f+slope);  // pmin^(1-g)
    float pmaxp = std::pow(pmax, 1.0f+slope);  // pmax^(1-g)

    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    auto mins = tile.mins;
    auto maxs = tile.maxs;

    float x_min = D >= 1 ? mins[0] : 0.0f;
    float y_min = D >= 2 ? mins[1] : 0.0f;
    float z_min = D >= 3 ? mins[2] : 0.0f;

    float lenx = D >= 1 ? (maxs[0] - mins[0]) : 0.0f;
    float leny = D >= 2 ? (maxs[1] - mins[1]) : 0.0f;
    float lenz = D >= 3 ? (maxs[2] - mins[2]) : 0.0f;

    // inject N_inj pairs
    float vx, vy, vz, xi;
    float ux, uy, uz, ginj, pinj;
    float xloc, yloc, zloc;

    float ncop = 0.0f;
    float z1 = rand();

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



  void comp_tau(
      pic::Tile<D>& tile, 
      double w2tau_units
      )
  {

    // build pointer map of types to containers; used as a helper to access particle types
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //--------------------------------------------------
    // pair number density
    float wsum_ep = 0.0f;
    float wvrel   = 0.0f;
    float w1, beta, gam;

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
    //float tauT  = wvrel*w2tau_units; // \tau_T = \sigma_T <v_rel> wsum
    float tauT = wsum_ep*w2tau_units; // \tau_T = \sigma_T wsum

    //std::cout << "escape" << std::endl;
    //std::cout << "   tauT: " << tauT << std::endl;
    //std::cout << "   tauT2:" << tauT2 << std::endl;
    //std::cout << "   wsum: " << wsum_ep << std::endl;
    //std::cout << "   N_-:  " << cons["e-"]->size() << std::endl;
    //std::cout << "   N_+:  " << cons["e+"]->size() << std::endl;
    //std::cout << "   N_w:  " << w2tau_units << std::endl;

    // increase global tau measure
    tau_global += tauT; 

    return;
  }

  //--------------------------------------------------
  // photon escape from the box w/ escape probability formalism
  void leak_photons(
      pic::Tile<D>& tile, 
      double tc_per_dt,
      double tau_ext
      )
  {

    // build pointer map of types to containers; used as a helper to access particle types
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //--------------------------------------------------
    // NOTE two variants; TODO which one is more correct for escape prob. formalism?
      
    //float tauT  = wvrel*w2tau_units; // \tau_T = \sigma_T <v_rel> wsum
    float tauT = tau_global; // global optical depth of whole domain

    if(tau_ext > 0.0) tauT = tau_ext; // use external tau if given

    //--------------------------------------------------
    std::string t1 = "ph";
    size_t Nx = cons[t1]->size(); // read particle number from here; 
                                  //
    //cons[t1]->to_other_tiles.clear(); // clear book keeping array

    float w, x, f, sKN; //, P_esc;
    for(size_t n1=0; n1<Nx; n1++) {
      w = cons[t1]->wgt(n1);
      x = cons[t1]->get_prtcl_ene(n1);

      // Klein-Nishina cross section
      if (x < 0.05f) { //Taylor expansion of KN
        sKN = (1.0f-(x*(20.0f+x*(133.0f*x-52.0f)))*0.1f);
      } else { // full cross section
        sKN = 0.75f*( 
          ( (1.0f + x)/(x*x*x))*( 2.0f*x*(1.0f + x)/(1.0f + 2.0f*x) - std::log(1.0f + 2.0f*x)) 
            + std::log(1.0f + 2.0f*x)/(2.0f*x)  - (1.0f + 3.0f*x)/pow(1.0f + 2.0f*x, 2) );
      }

      //sKN = 1.0f; // Thomson cross-section

      // empirical escape probability function to account for forward scattering pile-up
      // see Lightman \& Zdarskiaki 1987
      if(        x <= 0.1) { f = 1.0f;
      } else if (x >  1.0) { f = 0.0f;
      } else               { f = (1.0f - x)/0.9f; }

      //(c/R)*dt = dt/t_c
      //P_esc = dt_per_tc/( 0.75f + 0.188f*tauT*f*sKN ); // sphere
      float t_esc = ( 1.0f + tauT*f*sKN ); // slab
      //float t_esc = ( 1.0f + tauT*f*sKN + (1.0f-f)*2.886f ); // slab w/ asymptic scaling to 5/sqrt(3)
                                                               // this mimics pair-production opacity
                                                               // asymptotic solution to slab geometry
                                                               // with rad. transf. when tau -> infty
      //P_esc = t_esc/dt_per_tc; // P = R/c / dt for tau -> 0 

      // photon has an escape probability rate of P_esc = 1.0/(t_c*t_esc) to escape.
      // Probability is p_esc = P_esc*dt = dt/(t_c*t_esc) 

      //std::cout << "esc:" << 1.0f/tc_per_dt/t_esc << " " << t_esc << " " << tc_per_dt << std::endl;

      if(sKN > 1.0f || sKN < 0.0f){
        std::cout << "ERR: leak ph KN" << std::endl;
        std::cout << "  x:" << x << std::endl;
        std::cout << "SKN:" << sKN << std::endl;
        std::cout<< 1.3333333f*(1.0f-(x*(20.0f+x*(133.0f*x-52.0f)))*0.1f) << std::endl;

        sKN = 0.75f*( 
          ( (1.0f + x)/x*x*x)*( 2.0f*x*(1.0f + x)/(1.0f + 2.0f*x) - std::log(1.0f + 2.0f*x)) 
            + std::log(1.0f + 2.0f*x)/(2.0f*x)  - (1.0f + 3.0f*x)/pow(1.0f + 2.0f*x, 2) );
        std::cout << "KN:" << sKN << std::endl;
      }


      if( 1.0f/tc_per_dt/t_esc > rand() ) {
        leaked_ene  += x*w;
        leaked_wsum += w;
        leaked_pnum += 1;

        // book keeping of escaped flux
        add_to_histogram(x, w);

        cons[t1]->info(n1) = -1; //to_other_tiles.push_back( {1,1,1,n1} ); // NOTE: CPU deletion version
        cons[t1]->wgt(n1) = 0.0f; 
      }

    }

    //// remove annihilated prtcls; this transfer storage 
    cons[t1]->delete_transferred_particles(); 

    return;
  }


};



} // end of namespace qed
