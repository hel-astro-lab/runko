// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <atomic>
#include <cmath>
#include <type_traits>

#include "tyvi/backend.h"

#ifdef TYVI_BACKEND_HIP
#include "hip/hip_runtime.h"
#include "thrust/memory.h"
#endif

namespace sstd {

// =====================================================================
// raw_ref — identity on CPU, thrust::raw_reference_cast on HIP
// =====================================================================

template<typename T>
constexpr auto&
raw_ref(T& x) {
#ifdef TYVI_BACKEND_CPU
    return x;
#elif defined(TYVI_BACKEND_HIP)
    return thrust::raw_reference_cast(x);
#endif
}

// =====================================================================
// atomic_add — direct += on CPU, unsafeAtomicAdd on HIP
// =====================================================================

template<typename T>
constexpr void
atomic_add(T* addr, T val) {
#ifdef TYVI_BACKEND_CPU
    std::atomic_ref<T>(*addr).fetch_add(val, std::memory_order_relaxed);
#elif defined(TYVI_BACKEND_HIP)
    ::unsafeAtomicAdd(addr, val);
#endif
}

// =====================================================================
// atomic_max / atomic_min — for non-negative floating point values
// =====================================================================
// On HIP the bit pattern of a non-negative IEEE number orders like an
// unsigned integer of the same width, so the integer atomics serve.

template<typename T>
constexpr void
atomic_max(T* addr, T val) {
#ifdef TYVI_BACKEND_CPU
    std::atomic_ref<T> r(*addr);
    auto old = r.load(std::memory_order_relaxed);
    while (old < val and not r.compare_exchange_weak(old, val, std::memory_order_relaxed)) {}
#elif defined(TYVI_BACKEND_HIP)
    if constexpr (std::is_same_v<T, float>) {
        ::atomicMax(reinterpret_cast<unsigned int*>(addr), __float_as_uint(val));
    } else {
        ::atomicMax(reinterpret_cast<unsigned long long*>(addr),
                    static_cast<unsigned long long>(__double_as_longlong(val)));
    }
#endif
}

template<typename T>
constexpr void
atomic_min(T* addr, T val) {
#ifdef TYVI_BACKEND_CPU
    std::atomic_ref<T> r(*addr);
    auto old = r.load(std::memory_order_relaxed);
    while (old > val and not r.compare_exchange_weak(old, val, std::memory_order_relaxed)) {}
#elif defined(TYVI_BACKEND_HIP)
    if constexpr (std::is_same_v<T, float>) {
        ::atomicMin(reinterpret_cast<unsigned int*>(addr), __float_as_uint(val));
    } else {
        ::atomicMin(reinterpret_cast<unsigned long long*>(addr),
                    static_cast<unsigned long long>(__double_as_longlong(val)));
    }
#endif
}

// =====================================================================
// sin / cos — constexpr type-dispatched wrappers
// =====================================================================

template<typename T>
constexpr auto
sin(const T x) {
    if constexpr (std::is_same_v<T, float>) {
        return ::sinf(x);
    } else {
        return ::sin(x);
    }
}

template<typename T>
constexpr auto
cos(const T x) {
    if constexpr (std::is_same_v<T, float>) {
        return ::cosf(x);
    } else {
        return ::cos(x);
    }
}

// =====================================================================
// sqrt — constexpr type-dispatched wrapper
// =====================================================================

template<typename T>
constexpr auto
sqrt(const T x) {
    if constexpr (std::is_same_v<T, float>) {
        return ::sqrtf(x);
    } else {
        return ::sqrt(x);
    }
}

// =====================================================================
// min / max — SIMD-friendly value-based min/max
// =====================================================================
// Using ternary instead of std::min/std::max to avoid reference
// semantics that can inhibit auto-vectorization.

template<typename T>
constexpr auto
min(const T a, const T b) {
    return a < b ? a : b;
}

template<typename T>
constexpr auto
max(const T a, const T b) {
    return a > b ? a : b;
}

// =====================================================================
// clamp — branchless clamp to [lo, hi]
// =====================================================================

template<typename T>
constexpr auto
clamp(const T x, const T lo, const T hi) {
    return min(max(x, lo), hi);
}

// =====================================================================
// abs — SIMD-friendly absolute value
// =====================================================================

template<typename T>
constexpr auto
abs(const T x) {
    return x < T{ 0 } ? -x : x;
}

// =====================================================================
// fmod — constexpr type-dispatched wrapper
// =====================================================================

template<typename T>
constexpr auto
fmod(const T a, const T b) {
    if constexpr (std::is_same_v<T, float>) {
        return ::fmodf(a, b);
    } else {
        return ::fmod(a, b);
    }
}

// =====================================================================
// floor — SIMD-friendly floor (FRINTM on ARM NEON)
// =====================================================================

template<typename T>
constexpr auto
floor(const T x) {
    if constexpr (std::is_same_v<T, float>) {
        return ::floorf(x);
    } else {
        return ::floor(x);
    }
}

// =====================================================================
// ceil — SIMD-friendly ceil (FRINTP on ARM NEON)
// =====================================================================

template<typename T>
constexpr auto
ceil(const T x) {
    if constexpr (std::is_same_v<T, float>) {
        return ::ceilf(x);
    } else {
        return ::ceil(x);
    }
}

// =====================================================================
// log10 — constexpr type-dispatched wrapper
// =====================================================================

template<typename T>
constexpr auto
log10(const T x) {
    if constexpr (std::is_same_v<T, float>) {
        return ::log10f(x);
    } else {
        return ::log10(x);
    }
}


// =====================================================================
// exp / log / log1p / expm1 / tanh / atanh / cbrt / pow — type-dispatched
// =====================================================================

#define RUNKO_SSTD_UNARY(name)                       \
  template<typename T>                               \
  constexpr auto name(const T x) {                   \
    if constexpr (std::is_same_v<T, float>) {        \
      return ::name##f(x);                           \
    } else {                                         \
      return ::name(x);                              \
    }                                                \
  }

RUNKO_SSTD_UNARY(exp)
RUNKO_SSTD_UNARY(log)
RUNKO_SSTD_UNARY(log1p)
RUNKO_SSTD_UNARY(expm1)
RUNKO_SSTD_UNARY(tanh)
RUNKO_SSTD_UNARY(atanh)
RUNKO_SSTD_UNARY(cbrt)

#undef RUNKO_SSTD_UNARY

template<typename T>
constexpr auto
pow(const T x, const T y) {
    if constexpr (std::is_same_v<T, float>) {
        return ::powf(x, y);
    } else {
        return ::pow(x, y);
    }
}


/// Signum of value
template <typename T> int sign(const T& val) {
    return (T(0) < val) - (val < T(0));
}

} // namespace sstd
