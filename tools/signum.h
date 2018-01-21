namespace toolbox {


/// Signum of value
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


}
