#include <vector>
#include <algorithm>

template<class T>
bool has_elem(std::vector<T>& v, T x){
  if(std::find(v.begin(), v.end(), x) != v.end()) {
      return true;
  } else {
      return false;
  }
}

