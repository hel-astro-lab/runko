// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/tools/config_parser.h"

#include <optional>
#include <string>

namespace runko {

/// Output directory from argument or from config (io_outdir).
std::string
  resolve_outdir(const toolbox::ConfigParser&, const std::optional<std::string>&);

}  // namespace runko
