# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse


def print_conf(conf):
    for name in dir(conf):
        if name.startswith("_"):
            continue
        value = getattr(conf, name)
        if type(value) == str:
            print(f"{name} = \"{value}\"")
        else:
            print(f"{name} = {value}")


def main(argv: list[str]):
    parser = argparse.ArgumentParser(prog="runko inspect",
                                     description="Inspect runko output files.")

    parser.add_argument("outdir",
                        type=str,
                        help="Paths to runko output directory.")

    parser.add_argument("--config",
                        action="store_true",
                        help="""
                        Output contents of the configuration file of the given output directory.
                        """)


    args = parser.parse_args(argv[1:])

    outdir = args.outdir

    # Load heavy modules only after argument parsing.
    import runko
    conf = runko.read_config(outdir)
    print_conf(conf)
