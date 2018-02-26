#pragma once

namespace toolbox {


/// Signum of value
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

/*
 * Defining these typically leads to unambigious compiler linking errors
template <typename T> int sign(T& val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> int sign(const T& val) {
    return (T(0) < val) - (val < T(0));
}
*/


}
