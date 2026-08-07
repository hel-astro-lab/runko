// Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "tyvi/actions_ast.h"

namespace runko {

enum class symbol : std::uint8_t { print, println, version, mt_showcase };

[[nodiscard]]
tyvi::actions::sexpr build_stdenv();

}  // namespace runko
