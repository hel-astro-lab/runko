// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/simulation_context.h"
#include "tyvi/actions_ast.h"

#include <bit>
#include <cstddef>
#include <optional>
#include <string>

namespace runko {

static_assert(
  std::endian::native == std::endian::little,
  "MPI-IO binary format assumes little-endian");

/// Magic number: "RNKO" in ASCII
inline constexpr uint32_t magic = 0x524E4B4F;

/// Binary format version
inline constexpr uint32_t version = 3;

/// Fixed header size in bytes
inline constexpr int32_t header_size = 512;

/// Maximum number of particle species stored
inline constexpr int32_t max_species = 5;

/// Number of electromagnetic field quantities (ex, ey, ez, bx, by, bz, jx, jy, jz)
inline constexpr int32_t num_emf_fields = 9;

/// EMF field names in canonical order
inline constexpr const char* emf_field_names[9] = { "ex", "ey", "ez", "bx", "by",
                                                    "bz", "jx", "jy", "jz" };

/// Number of particle snapshot fields (x,y,z,ux,uy,uz,ex,ey,ez,bx,by,bz)
inline constexpr int32_t num_prtcl_fields = 12;

/// Particle field names in canonical order
inline constexpr const char* prtcl_field_names[12] = { "x",  "y",  "z",  "ux",
                                                       "uy", "uz", "ex", "ey",
                                                       "ez", "bx", "by", "bz" };

/// Maximum number of particle species for spectra
inline constexpr int32_t max_spectra_species = 3;

/// Number of spectral quantities per species (u, beta_x, beta_y, beta_z)
inline constexpr int32_t num_spectra_per_species = 4;

/// Spectral quantity suffixes
inline constexpr const char* spectra_suffixes[4] = { "u", "bx", "by", "bz" };


tyvi::actions::sexpr_sender
  emf_snapshot(runko::simulation_context&, long lap, std::optional<std::string>);

/// Write sampled particles of each species to {outdir}/prtcls_{species}_{lap}.bin.
///
/// Number of sampled particles is read from config (io_n_sampled_prtcls),
/// if it is not given explicitly. Same for outdir (io_outdir).
tyvi::actions::sexpr_sender prtcl_snapshot(
  runko::simulation_context&,
  long lap,
  std::optional<std::string> outdir,
  std::optional<long> n_sampled_prtcls);
tyvi::actions::sexpr_sender spectra_snapshot(runko::simulation_context&, long lap);

}  // namespace runko
