import glob
import re
import os
import argparse

from configSetup import Configuration


def parse_input():

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--conf', 
            dest='conf_filename', 
            default=None,
            type=str,
            help='Name of the configuration file (default: None)')

    parser.add_argument('-l', '--lap', 
            dest='lap', 
            default=None,
            type=int,
            help='Specific lap to analyze (default: all that can be found)')

    parser.add_argument('-d', '--dir', 
            dest='fdir', 
            default=None,
            type=str,
            help='Directory to analyze')

    parser.add_argument('-v', '--var', 
            dest='var', 
            default=None,
            type=str,
            help='Variable to analyze')

    args = parser.parse_args()
    #args = parser.parse_known_args()

    if not(args.conf_filename == None):
        conf = Configuration(args.conf_filename)
        fdir = conf.outdir

    if not(args.fdir == None):
        fdir = args.fdir

        #look for conf file inside dir
        if args.conf_filename == None:
            print("looking for Conf file inside {}".format(args.fdir))
            print("TODO: not implemented")


    return conf, fdir, args




