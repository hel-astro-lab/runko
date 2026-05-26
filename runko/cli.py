# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse
import sys

help_str = """usage: runko <subcommand> [ARGS...]

Invokes <subcommand> with ARGS.

Available subcommands: help
                       weak-scaling"""


def main():

    if len(sys.argv) < 2:
        print(help_str)
        quit(1)

    subcommand = sys.argv[1]

    if subcommand == "help":
        print(help_str)
        quit()
    elif subcommand == "weak-scaling":
        import runko.subcommands.weak_scaling
        runko.subcommands.weak_scaling.main(sys.argv[1:])
    else:
        print(f"error: unregonized subcommand: '{subcommand}'")
        print(help_str)
        quit(1)
