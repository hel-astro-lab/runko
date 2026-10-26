// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#include "runko/io/fs_utils.h"

std::string
  runko::resolve_outdir(
    const toolbox::ConfigParser& config,
    const std::optional<std::string>& outdir)
{
  if(outdir) { return outdir.value(); }

  const auto x = config.get<std::string>("io_outdir");
  if(not x or x.value().empty()) { return "runko_output"; }
  if(x.value() == "auto") {
    throw std::logic_error {
      "resolve_outdir: can not resolve `auto` outdir. It should have already been "
      "resolved in python. "
      "Please, submit a bug report about this."
    };
  }
  return x.value();
}
