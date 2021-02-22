# -*- coding: utf-8 -*- 

import argparse
import os


# top level input argument parser
def parse_args():

    # create top level parser
    parser = argparse.ArgumentParser(description="Runko plasma simulation framework")

    # configuration file
    parser.add_argument(
        "--conf",
        dest="conf_filename",
        default=None,
        help="Name of the configuration file (default: None)",
    )
    parser.add_argument('-l', '--lap', 
            dest='lap', 
            default=None,
            type=int,
            help='Specific lap to analyze (default: all that can be found)')

    parser.add_argument('-v', '--var', 
            dest='var', 
            default=None,
            type=str,
            help='Variable to analyze')


    # --------------------------------------------------
    # first high level parse
    args = parser.parse_args()

    return args


def create_output_folders(conf):

    if not os.path.exists(conf.outdir):
        os.makedirs(conf.outdir)
    if not os.path.exists(conf.outdir + "/restart"):
        os.makedirs(conf.outdir + "/restart")
    if not os.path.exists(conf.outdir + "/full_output"):
        os.makedirs(conf.outdir + "/full_output")

    # finally, copy current configuration to the folder
    #if True:
    try:
        # only try this because we dont want to require shutils
        import shutil
        shutil.copyfile(
                conf.conf_filename, 
                conf.outdir + "/" + conf.conf_filename)
    except:
        pass


# check for existing restart files; 
# manage io_status object accordingly with current lap, restart directory, etc...
def check_for_restart(conf):

    io_status = {}

    # set flag to conf if we need to initialize or not
    io_status["do_initialization"] = True

    # check if this is the first time and we do not have any restart files
    if not os.path.exists(conf.outdir + "/restart/laps.txt"):
        conf.laprestart = -1  # to avoid next if statement

    # restart from latest file
    io_status["deep_io_switch"] = 0
    io_status["restart_num"] = 0

    if conf.laprestart >= 0:
        io_status["do_initialization"] = False

        # switch between automatic restart and user-defined lap

        # zero signals automatic checking of latest lap
        if conf.laprestart == 0:
            #print("restarting automatically...")

            # get latest restart file from housekeeping file
            with open(conf.outdir + "/restart/laps.txt", "r") as lapfile:
                # lapfile.write("{},{}\n".format(lap, deep_io_switch))
                lines = lapfile.readlines()
                slap, sdeep_io_switch = lines[-1].strip().split(",")

                # keep track of number of restarts by encoding them to deep io switch
                max_ios = 0
                for line in lines:
                    a, b = line.strip().split(",")
                    if int(b) > max_ios:
                        max_ios = int(b)

                io_status["restart_num"] = max_ios + 1

                #TODO think about some way to cycle restart_num to save disk space

                io_status["lap"] = int(slap)
                io_status["deep_io_switch"] = int(sdeep_io_switch)

            io_status["read_lap"] = io_status['deep_io_switch']
            io_status["read_dir"] = conf.outdir + "/restart"

        # >0 lap means that we need to restart from full_output
        elif conf.laprestart > 0:
            io_status["lap"] = conf.laprestart
            io_status["read_lap"] = lap
            io_status["read_dir"] = conf.outdir + "/full_output"


    return io_status
