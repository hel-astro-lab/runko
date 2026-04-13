#pragma once

#include <cstdint>
#include "tyvi/backend.h"

#if defined(TYVI_BACKEND_HIP)
#include <hip/hip_runtime.h>
#endif

namespace runko {

/// Return GPU device memory in use (kB), or -1 when no GPU backend.
inline int64_t get_gpu_mem_kB()
{
#if defined(TYVI_BACKEND_HIP)
  size_t free_bytes = 0, total_bytes = 0;
  if (hipMemGetInfo(&free_bytes, &total_bytes) == hipSuccess) {
    return static_cast<int64_t>((total_bytes - free_bytes) / 1024);
  }
  return -1;
#else
  return -1;
#endif
}

} // namespace runko
