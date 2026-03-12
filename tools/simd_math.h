#pragma once

#include <atomic>
#include <cmath>
#include <cstdint>
#include <type_traits>

#include "tyvi/backend.h"

// =====================================================================
// runko::index_t — 32-bit index type to avoid 64-bit ↔ float promotion
// that halves SIMD throughput (size_t 64-bit → int → double × float).
// =====================================================================

namespace runko {
using index_t = uint32_t;
}  // namespace runko

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
inline void
atomic_add(T* addr, T val) {
#ifdef TYVI_BACKEND_CPU
    std::atomic_ref<T>(*addr).fetch_add(val, std::memory_order_relaxed);
#elif defined(TYVI_BACKEND_HIP)
    ::unsafeAtomicAdd(addr, val);
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

} // namespace sstd
