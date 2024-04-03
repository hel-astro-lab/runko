#pragma once

#include <chrono>
#include <vector>
#include <numeric>
#include <map>
#include <string>

#define FMT_HEADER_ONLY
#include "../fmt/include/fmt/format.h"
#include "../fmt/include/fmt/ranges.h"

// NOTE not defining this will cause the timer to not do anything
// #define USE_INTERNAL_TIMER


// helper script to return if map contains key
template<typename T, typename S>
inline bool contains(std::map<T, S>& map, T val)
{
  auto it = map.find(val);
  return it != map.end();
}

//--------------------------------------------------
namespace math {
    
// sum
template<typename T>
inline T sum(std::vector<T>& v){
  return std::accumulate(v.begin(), v.end(), T(0));
}

// mean
template<typename T>
inline T mean(std::vector<T>& v){
  return sum(v)/v.size();
}

// std
template<typename T>
inline T std(std::vector<T>& v){
  T avg = mean(v);
  T sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), T(0));
  return std::sqrt(sq_sum/v.size() - avg*avg);
}

} // end of ns math
//--------------------------------------------------



class Timer {

  public:

  // printing options
  bool do_print = true;
  int verbose = 0;

  using clock        = std::chrono::high_resolution_clock;
  using time_point_t = std::chrono::time_point<clock>;
  using time_vec_t   = std::vector<time_point_t>;
  using dict_t       = std::map< std::string, time_vec_t>;
  using duration_t   = std::chrono::nanoseconds;

  // dictionary of time components
  dict_t components_beg;
  dict_t components_end;

  time_point_t s1; // start of main block session
  time_point_t s2; // end of main block session

  duration_t duration_on{0}; // duration that the timer has been on
  int num_cycles_measured = 0; // number of times the block has been switched on/off

  std::string timer_name; 

  //--------------------------------------------------
  Timer(std::string timer_name) : timer_name(timer_name)
  { }

  void start()
  {
#if defined(USE_INTERNAL_TIMER)
    s1 = clock::now();
#endif
  }

  void stop()
  {
#if defined(USE_INTERNAL_TIMER)
    const auto s2 = clock::now();

    // add elapsed time to the internal clock
    duration_t dur{s2 - s1};
    duration_on += dur; 

    num_cycles_measured++; // increase clock counter
#endif
  }


  //--------------------------------------------------
  std::string start_comp(std::string name)
  {
#if defined(USE_INTERNAL_TIMER)
    time_point_t t0 = clock::now(); // get time

    // create new comp if not there
    if( !contains(components_beg, name) ) components_beg[name] = time_vec_t();
     
    // add time
    components_beg[name].push_back(t0);
#endif
    return name;
  }

  void stop_comp(std::string name)
  {
#if defined(USE_INTERNAL_TIMER)
    time_point_t t0 = clock::now(); // get time
    
    // create new comp if not there
    if( !contains(components_end, name) ) components_end[name] = time_vec_t();

    components_end[name].push_back(t0);
#endif
  }

  //--------------------------------------------------

  void comp_stats()
  {
#if defined(USE_INTERNAL_TIMER)

    double t0tot = duration_on.count();

    if(do_print) fmt::print("------------------------------------------------------------------------\n");
    if(do_print) fmt::print(" {:<}: {:8.5f} s ({:>.8} ns)      | cycles: {:>5d} \n", 
        timer_name, t0tot/1.0e9, t0tot, num_cycles_measured);


    double totper = 0.0; // check that percentage sums up to 100
                         //
    for(auto const & [name, beg] : components_beg) 
    {
      if(beg.empty()) continue;
      auto& end = components_end[name];
      if(end.empty()) continue;


      // cast to nanoseconds
      const int N = beg.size();
      std::vector<double> ts(N);
      for(int i=0; i<N; i++) ts[i] = std::chrono::nanoseconds( end[i] - beg[i] ).count();
      
      double tavg = math::mean(ts);
      int    cnts = ts.size();
      //double tstd = math::std(ts);

      double t0 = math::sum(ts);
      double relt = 100.0*t0/t0tot;
      totper += relt;

      //fmt::print("--- {:>20} {} \n", name, ts);

      if(        tavg < 1e-5*1e9 ){

        if(do_print) fmt::print("--- {:<20}   {:6.3f}%  |  time: {:8.5f} mus/{:<8d}  ({:>.4g} ns)\n",
                   name, relt, tavg*1.0e6/1.0e9, cnts, tavg);

      } else if( tavg < 1.0e-2*1e9 ){

        if(do_print) fmt::print("--- {:<20}   {:6.3f}%  |  time: {:8.5f} ms /{:<8d}  ({:>.4g} ns)\n",
                   name, relt, tavg*1.0e3/1.0e9, cnts, tavg);

      } else {

        if(do_print) fmt::print("--- {:<20}   {:6.3f}%  |  time: {:8.5f} s  /{:<8d}  ({:>.4g} ns)\n",
                   name, relt, tavg/1.0e9, cnts, tavg);

      }


    } // end of components loop

    if(do_print) fmt::print("                        += {:6.3f}% \n", totper);
    if(do_print) fmt::print("------------------------------------------------------------------------\n");
#endif
  }


  void clear()
  {
    duration_on = duration_t{0};
    num_cycles_measured = 0;

    for(auto const & [key, val] : components_beg) components_beg[key] = time_vec_t();
    for(auto const & [key, val] : components_end) components_end[key] = time_vec_t();
  }


};

