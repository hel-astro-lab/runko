// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "pika/mpi.hpp"

#include <memory>

namespace runko {

struct unmovable {
  constexpr unmovable()                  = default;
  constexpr ~unmovable()                 = default;
  unmovable(unmovable&&)                 = delete;
  unmovable(const unmovable&)            = delete;
  unmovable& operator=(unmovable&&)      = delete;
  unmovable& operator=(const unmovable&) = delete;
};

/// Raii class for runko runtime.
///
/// Supports multiple concurrent runtime instances.
/// However, only one runtime is initialized
/// at the construction of the first overlapping instance
/// and deinitialized at the destruction of the last overlapping instance.
struct [[nodiscard]] RuntimeInstance : unmovable {
  RuntimeInstance();

private:
  struct impl {
    impl();
    ~impl();

    /// We initialize mpi if it is not already been initialized.
    ///
    /// And if we do so, we finalize it at the end.
    /// In the long term, we would like to use MPI sessions,
    /// so we don't have to worry about who is responsible for MPI.
    /// However, pika does not support them yet.
    bool initialized_mpi_ { false };
  };

  std::shared_ptr<impl> impl_;
};

/// RAII helper for starting and suspending runtime.
struct RuntimeActivator : unmovable {
  pika::mpi::experimental::enable_polling enable_polling {};

  RuntimeActivator();
  ~RuntimeActivator();
};

}  // namespace runko
