// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/runtime.h"

#include <iostream>
#include <mpi.h>
#include <pika/init.hpp>
#include <pika/mpi.hpp>
#include <print>
#include <stdexcept>
#include <string_view>

namespace runko {

RuntimeInstance::RuntimeInstance()
{
  int mpi_flag;
  if(MPI_Initialized(&mpi_flag) != MPI_SUCCESS) {
    throw std::runtime_error { "MPI_Initialized failed!" };
  }
  if(not static_cast<bool>(mpi_flag)) {
    int provided {}, preferred = pika::mpi::experimental::get_preferred_thread_mode();

    MPI_Init_thread(nullptr, nullptr, preferred, &provided);
    if(provided != preferred) {
      throw std::runtime_error {
        "Provided level of thread support in MPI is not as requested."
      };
    }
    this->initialized_mpi_ = true;
  }

  if(pika::is_runtime_initialized()) {
    throw std::runtime_error {
      "Trying to initialize runko runtime while pika has already been initialized."
    };
  }
  pika::start(0, nullptr);
  pika::suspend();
}

RuntimeInstance::~RuntimeInstance()
{
  if(pika::is_runtime_initialized()) {
    // For some reason we have to resume before finalizing pika.
    pika::resume();
    pika::finalize();
    pika::stop();
  }

  if(this->initialized_mpi_) {
    if(MPI_SUCCESS != MPI_Finalize()) {
      std::println(std::cerr, "MPI_Finalize failed.");
      std::terminate();
    }
  }
}

RuntimeActivator::RuntimeActivator() { pika::resume(); }

RuntimeActivator::~RuntimeActivator() { pika::suspend(); }
}  // namespace runko
