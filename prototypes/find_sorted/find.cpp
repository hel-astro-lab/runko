
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <cassert>

#define FMT_HEADER_ONLY
#include "../fmt/include/fmt/format.h"
#include "../fmt/include/fmt/ranges.h"


#define USE_INTERNAL_TIMER
#include "../timer/timer.h"

//#include <immintrin.h>


// TODO this will work with c++-20 for any container that satisfies common range
//size_t iterate_container(const std::ranges::common_range auto& container) 
//{ 
//  for(const auto& item : container) { /* do something */}:
//}


// forward find of reverse sorted array
// linear search that is O(N)
template <typename T>
inline size_t find_rev_sorted_nearest( const std::vector<T>& arr, const T val) {

  const size_t arr_size = arr.size();
  if(arr_size == 0) return 0;

  // value is between array elements
  size_t i=0u;
  for(; i<arr_size; i++) if( val > arr[i] ) return i;

  // last element
  return arr_size; 
}


// reverse find of reverse sorted array
// linear search that is O(N)
template <typename T>
inline size_t revfind_rev_sorted_nearest_short( const std::vector<T>& arr, const T val) {

  if(arr.size() == 0) return 0;

  // special case where we would have to wrap all around the array
  if(val > arr[0u]) return 0u;

  const size_t N = arr.size();
  if(val <= arr[N-1]) return N;

  for(size_t i=N; i>=0u; i--){
    //fmt::print("    m3: iterating {} comparing {} < {}\n",i, val, arr[i]);
    if( val <= arr[i] ) return i+1;
  }

  return 0u;
}

// reverse find of reverse sorted array
// uses binary search that is O(log(N))
template <typename T>
inline size_t find_rev_sorted_nearest_algo( const std::vector<T>& arr, const T val) {
    auto it = std::lower_bound( arr.cbegin(), arr.cend(), val, std::greater_equal<T>() );
    return it - arr.begin();
}



// Sample between [imin, imax[
template <typename T>
size_t inline sample_prob_between( 
    const std::vector<T>& ws, 
    const T val, 
    const size_t imin, 
    const size_t imax) 
{
    T wmin, wmax;
    if(imin == 0u) {
        wmin = static_cast<T>(0); // # equal to ws[-1]
    } else {
        wmin = ws[imin-1u];
    }
    wmax = ws[imax-1u];
    //const wminmax = wmax - wmin;
    const T vminmax = val*(wmax - wmin) + wmin;

    for(size_t j=imin; j<imax; j++){
        //if( val < (ws[j] - wmin)/(wmax - wmin) ) return j;
        //if( val < (ws[j] - wmin)/wminmax ) return j;
        if( vminmax < ws[j] ) return j;
    }
    return imax;
}

// Sample between [imin, imax[
template <typename T>
size_t inline sample_prob_between_algo( 
    const std::vector<T>& ws, 
    const T val, 
    const size_t imin, 
    const size_t imax) 
{
    T wmin, wmax;
    if(imin == 0u) {
        wmin = static_cast<T>(0); // # equal to ws[-1]
    } else {
        wmin = ws[imin-1u];
    }
    wmax = ws[imax-1u];
    const T vminmax = val*(wmax - wmin) + wmin;

    // search in the range [imin, imax] using binary search
    auto beg = std::next(ws.cbegin(), imin);
    auto end = std::next(ws.cbegin(), imax);
    auto it = std::upper_bound(beg, end, vminmax);

    return it - ws.cbegin();
}


// adapted from https://en.algorithmica.org/hpc/data-structures/binary-search/
//size_t lower_bound(
//    const std::vector<float>& base,
//    float x) 
//{
//
//    int *base = t, len = n;
//    while (len > 1) {
//        int half = len / 2;
//        len -= half;
//        __builtin_prefetch(&base[len / 2 - 1]);
//        __builtin_prefetch(&base[half + len / 2 - 1]);
//        base += (base[half - 1] < x) * half;
//    }
//    return *base;
//}

// based on https://en.algorithmica.org/hpc/data-structures/binary-search/
int lower_bound1(
    std::vector<float>& t,
    const float x) 
{

  const size_t N = t.size();
  if(x < t[N-1]) return N;

  float *base = t.data();
  int len = N;
  while (len > 1) {
        int half = len / 2;

        base += (base[half - 1] >= x) * half; // will be replaced with a "cmov"
        len -= half;
    }
    //return *base; // value
    return base - &t[0]; // index
}


int lower_bound2(
    std::vector<float>& t,
    const float x) 
{

  const size_t N = t.size();
  if(x < t[N-1]) return N;

  float *base = t.data();
  int len = N;
  while (len > 1) {
        int half = len / 2;
        len -= half;

        __builtin_prefetch(&base[len / 2 - 1]);
        __builtin_prefetch(&base[half + len / 2 - 1]);
        base += (base[half - 1] >= x) * half; // will be replaced with a "cmov"
    }
    //return *base; // value
    return base - &t[0]; // index
}


// eytzinger + prefetching
int lower_bound3(
    std::vector<float>& arr,
    const float x) 
{
  const size_t n = arr.size();
  if(x < arr[n-1]) return n;

  float* t = arr.data();
  int k = 1;
  while (k <= n) {
    __builtin_prefetch(t + k*16);
    k = 2 * k + (t[k] < x);
  }
  k >>= __builtin_ffs(~k);
  return k;
}


int main ()
{


  std::vector<float> arrsmall{10,9,8,7,6,5,4,3,2,1,0};

  fmt::print("{}\n", arrsmall);
  std::vector<float> vs{20., 10., 9.0, 9.5, 5.0, 2.0, 1.4, 0.1, 0.0, -0.5};

  int j1,j2,j3,j4;
  for(float v : vs){
    fmt::print(" value {:f} ---------------\n", v);

    j1 = find_rev_sorted_nearest<float>(arrsmall, v);
    j2 = revfind_rev_sorted_nearest_short<float>(arrsmall, v);
    j3 = find_rev_sorted_nearest_algo(arrsmall, v);
    j4 = lower_bound1(arrsmall, v);
    j4 = lower_bound3(arrsmall, v);

    fmt::print("j {} {} {} {}\n",j1,j2,j3,j4);
  }


  //std::vector<float> arr{0,0,1,2,3,4,5,5,5,6,7,10,10,11,12,13,14,19,20,20};
  //fmt::print("{}\n", arr);

  //std::vector<float> vs{0, 0.1, 0.2, 0.5, 0.8, 0.99, 1.0};

  //int i1, i2;
  //for(float v : vs){
  //  fmt::print(" value {:f} ---------------\n", v);

  //  i1 = sample_prob_between(     arr, v, 0, arr.size());
  //  i2 = sample_prob_between_algo(arr, v, 0, arr.size());

  //  fmt::print("v {} | j {} {}\n",v, i1, i2);
  //}



  //--------------------------------------------------
  //profiling

  Timer timer("search");


  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<float> uni_dis(0.0f, 1.0f);

  const size_t N = 2000u;
  std::vector<float> arr(N);
  //std::vector<__m256> arr(N);
  for(size_t i=0; i<N; i++) arr[i] = 1.0e3*uni_dis(gen);


  const size_t N_samples = 1000000u;
  std::vector<float> vals(N_samples);
  for(size_t i=0; i<N_samples; i++) vals[i] = 1.05e3*uni_dis(gen) - 25.0f;


  std::sort(arr.begin(), arr.end(), std::greater<float>() );
  //fmt::print("{}", arr);

  std::vector<size_t> ind1(N_samples); // results
  std::vector<size_t> ind2(N_samples); // results
  std::vector<size_t> ind3(N_samples); // results
  std::vector<size_t> ind4(N_samples); // results
                                 

  timer.start();


  timer.start_comp("linear");
  for(size_t i=0; i<N_samples; i++){
    ind1[i] = find_rev_sorted_nearest(arr, vals[i]);
  }
  timer.stop_comp("linear");


  timer.start_comp("std::algo");
  for(size_t i=0; i<N_samples; i++){
    ind2[i] = find_rev_sorted_nearest_algo(arr, vals[i]);
  }
  timer.stop_comp("std::algo");


  timer.start_comp("algor1");
  for(size_t i=0; i<N_samples; i++){
    ind3[i] = lower_bound1(arr, vals[i]);
  }
  timer.stop_comp("algor1");

  timer.start_comp("algor2");
  for(size_t i=0; i<N_samples; i++){
    ind4[i] = lower_bound2(arr, vals[i]);
  }
  timer.stop_comp("algor2");

  timer.stop();
  timer.comp_stats();
  timer.clear();


  fmt::print("arr stats are minmax {} {}\n", arr[0], arr[arr.size()-1]);

  for(size_t i=0; i<N_samples; i++){
    if(ind1[i] != ind2[i]) fmt::print("error {} vs {} at {} with val {}\n", ind1[i], ind2[i], i, vals[i]);
    if(ind1[i] != ind3[i]) fmt::print("error {} vs {} at {} with val {}\n", ind1[i], ind3[i], i, vals[i]);
    if(ind1[i] != ind4[i]) fmt::print("error {} vs {} at {} with val {}\n", ind1[i], ind4[i], i, vals[i]);

    assert(ind1[i] == ind2[i]);
    assert(ind1[i] == ind3[i]);
    assert(ind1[i] == ind4[i]);
  }




}



