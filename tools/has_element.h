// Copyright 2025 - 2026, Joonas Nättilä and the hel-astro-lab contributors
// SPDX-License-Identifier: GPL-3.0-or-later

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

