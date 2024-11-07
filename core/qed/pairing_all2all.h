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

#include "core/qed/interactions/interaction.h"


namespace qed {
  using std::string;
  using std::tuple;
  using std::min;
  using std::max;



//--------------------------------------------------
// Standard binary all2all pairing of particles
//
// This is a simplified (and naive) pairing object for two-body binary interactions.
// It performs the pairing by calculating all the possible combinations, and hence does
// all-to-all calculation. 
//
// For real (non-test) applications, see the full Monte Carlo implmementation in pairing.h


template<size_t D>
class PairingAll2All
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

  // constructor with incident/target types
  PairingAll2All() :
    gen(42), // gen(rd() ) 
    uni_dis(0.0, 1.0)
  { 
    update_hist_lims(hist_emin, hist_emax, hist_nbin);
  }

  //using Tile_map = std::unordered_map<TileID_t, Tileptr>;
  //Tile_map tiles; /// Map with tile_id & tile data

  // one-body single interactions
  //std::vector<InteractionPtr> single_interactions;

  // two-body binary interactions
  std::vector<InteractionPtr> binary_interactions;

  // normalization factor for two-body interaction probabilities
  float prob_norm = 1.0f;

  // normalization factor for single-body interaction probabilities
  //float prob_norm_onebody = 1.0f;

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

      assert(false); // TODO not implemented
      //single_interactions.push_back(iptr);

    //-------------------------------------------------- 
    } else if( iptr->interaction_order == 2 ){ // two-body binary interactions

      auto name = iptr->name;
      auto t1 = iptr->t1;
      auto t2 = iptr->t2;
      auto long_name = name + "_" + t1 + "_" + t2;

      //std::cout << " adding: " << name << " of t1/t2 " << t1 << " " << t2 << std::endl;
      binary_interactions.push_back(iptr);
    }

  }

  //--------------------------------------------------
  // check if interaction list is empty for type t1

  //bool is_empty_single_int(string& t1) 
  //{
  //  int i=0;
  //  for(auto iptr : single_interactions){
  //    if(t1 == iptr->t1) i += 1;
  //  }
  //  return (i > 0) ? false : true; 
  //}

  bool is_empty_binary_int(string& t1) 
  {
    int i=0;
    for(auto iptr : binary_interactions){
      if(t1 == iptr->t1) i += 1;
    }
    return (i > 0) ? false : true; 
  }

  // duplicate particle info into fresh variables
  inline auto duplicate_prtcl(
      string t1, float ux1, float uy1, float uz1, float w1
      ) -> std::tuple<string, float, float, float, float>
  {
    return {t1, ux1, uy1, uz1, w1};
  }


  //--------------------------------------------------
  //all-to-all binary comparison of particles and all processes
  // NOTE: this is very expensive...
  void solve_twobody(pic::Tile<D>& tile)
  {
      
    // build pointer map of types to containers; used as a helper to access particle tyeps
    std::map<std::string, ConPtr> cons;
    for(auto&& con : tile.containers) cons.emplace(con.type, &con );

    //--------------------------------------------------
    // call pre-iteration functions to update internal arrays 
    for(auto&& con : tile.containers)
    {
      con.sort_in_rev_energy();
      //con.update_cumulative_arrays();
      //con.to_other_tiles.clear(); // empty tmp container; we store killed particles here
    }

    //--------------------------------------------------
    // loop over interactions
    for(auto iptr : binary_interactions){

      //--------------------------------------------------
      // loop over incident types
      for(auto&& con1 : tile.containers)
      {
        auto t1 = con1.type;
        if(is_empty_binary_int(t1)) continue; // no interactions with incident type t1

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
            for(size_t n1=0; n1<con1.size(); n1++) { 

              // loop over targets
              //for(int n2=con2.size()-1; n2>=0; n2--) {
              for(size_t n2=0; n2<con2.size(); n2++) { 

                // NOTE: incident needs to be unpacked in the innermost loop, since 
                // some interactions modify its value during the iteration
                  
                //unpack incident 
                auto lx1 = con1.loc(0,n1);
                auto ly1 = con1.loc(1,n1);
                auto lz1 = con1.loc(2,n1);

                auto ux1 = con1.vel(0,n1);
                auto uy1 = con1.vel(1,n1);
                auto uz1 = con1.vel(2,n1);
                auto w1  = con1.wgt(n1);

                //e1  = con1.eneArr[n1]; 
                //auto e1  = con1.get_prtcl_ene(n1);

                if(w1 < EPS) continue; // omit zero-w incidents

                // unpack target
                auto lx2 = con2.loc(0,n2);
                auto ly2 = con2.loc(1,n2);
                auto lz2 = con2.loc(2,n2);

                auto ux2 = con2.vel(0,n2);
                auto uy2 = con2.vel(1,n2);
                auto uz2 = con2.vel(2,n2);
                auto w2  = con2.wgt(  n2);

                //auto e2 = con2.get_prtcl_ene(n2);

                if(w2 < EPS) continue; // omit zero-w targets
                                         

                //--------------------------------------------------
                //avoid double counting by considering only e1 < e2 cases
                //if(e1 > e2) continue;

                // interaction cross section
                auto [cm, vrel] = iptr->comp_cross_section(t1, ux1, uy1, uz1,  t2, ux2, uy2, uz2 );
                
                // interaction probability
                //auto wmin = min(w1, w2);
                auto wmax = max(w1, w2);
                auto prob = cm*vrel*w1*w2;
                // NOTE: difference of all2all scheme is here where prob depends on w1*w2

                // exponential waiting time between interactions
                float t_free = -log( rand() )*prob_norm/prob;

                //-------------------------------------------------- 
                if(t_free < 1.0){

                  // particle values after interaction
                  auto [t3, ux3, uy3, uz3, w3] = duplicate_prtcl(t1, ux1, uy1, uz1, w1);
                  auto [t4, ux4, uy4, uz4, w4] = duplicate_prtcl(t2, ux2, uy2, uz2, w2);

                  // interact and udpate variables in-place
                  iptr->interact( t3, ux3, uy3, uz3,  t4, ux4, uy4, uz4);

                  auto p_ini = w2/wmax;
                  auto p_tar = w1/wmax;

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
                      con1.info(n1) = -1; //to_other_tiles.push_back( {1,1,1,n1} ); // NOTE: CPU version
                      con1.wgt(n1) = 0.0f; // make zero wgt so its omitted from loop

                      // add new
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
                      con2.info(n2) = -1; //to_other_tiles.push_back( {1,1,1,n2} ); // NOTE: CPU version
                      con2.wgt(n2) = 0.0f; // make zero wgt so its omitted from loop

                      // add_particle
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
  // normalize container of type t1
  void rescale(pic::Tile<D>& tile, string& t1, double f_kill)
  {

    // build pointer map of types to containers; used as a helper to access particle tyeps
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

    // build pointer map of types to containers; used as a helper to access particle tyeps
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

    // build pointer map of types to containers; used as a helper to access particle tyeps
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
