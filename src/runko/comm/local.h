// Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "runko/communication_common.h"
#include "runko/emf/common.h"
#include "tyvi/actions.h"

namespace runko {

tyvi::actions::sexpr_sender comm_local(runko::simulation_context&, runko::comm_mode);

}  // namespace runko
