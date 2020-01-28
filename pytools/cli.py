import argparse

# top level input argument parser
def parse_args():

    #create top level parser
    parser = argparse.ArgumentParser(
            description='Runko plasma simulation framework'
            )

    # configuration file
    parser.add_argument(
            '--conf', 
            dest='conf_filename', 
            default=None,
            help='Name of the configuration file (default: None)')


    #--------------------------------------------------
    # first high level parse
    args = parser.parse_args()







    return args


