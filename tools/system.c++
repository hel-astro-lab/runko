#include "tools/system.h"

#include <cstdlib>
#include <string_view>

[[nodiscard]]
bool
  toolbox::system_supports_gpu_aware_mpi()
{
  static constexpr auto mpich_gpu_support_name   = "MPICH_GPU_SUPPORT_ENABLED";
  static constexpr auto force_gpu_aware_mpi_name = "RUNKO_FORCE_GPU_AWARE_MPI";


  static bool enable_gpu_aware_mpi = [] {
    if(const auto force_gpu_aware_mpi_var = std::getenv(force_gpu_aware_mpi_name)) {
      if(std::string_view { force_gpu_aware_mpi_var } == "1") { return true; }
    }
    if(const auto mpich_gpu_support_var = std::getenv(mpich_gpu_support_name)) {
      if(std::string_view { mpich_gpu_support_var } == "1") { return true; }
    }

    return false;
  }();

  return enable_gpu_aware_mpi;
}
