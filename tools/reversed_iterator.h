#pragma once


// --- Reversed iterable
// see https://stackoverflow.com/questions/8542591/c11-reverse-range-based-for-loop
using std::rbegin;
using std::rend;


template <typename T>
struct reversion_wrapper { T& iterable; };

template <typename T>
auto begin (reversion_wrapper<T> w) { return rbegin(w.iterable); }

template <typename T>
auto end (reversion_wrapper<T> w) { return rend(w.iterable); }

template <typename T>
reversion_wrapper<T> reverse (T&& iterable) { return { iterable }; }


