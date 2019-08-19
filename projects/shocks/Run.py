import os
import time
from parser import parse_input

def check():
    choice = input("Would you like to run analytical tools for Runko? (Y/n)")
    if (choice == "Y" or choice == "y" or choice == "N" or choice == "n"):
        x = 1
    else:
        print("Invalid input, please enter (Y/n)")
        x = 0
    return choice, x

def checkc():
    cores = input("How many cores do you wish to use?")
    if (cores.isdigit()):
        print("Running with {} cores...".format(cores))
        y = 1
    else:
        print("Invalid input, please enter a number")
        y = 0
    return cores, y

if __name__ == "__main__":
    print("-------------RUNKO---------------")
    conf, fdir, args = parse_input()
    x = 0
    while x == 0:
        choice, x = check()
    y = 0
    while y == 0:
        cores, y = checkc()

    print("------------------------------")
    print("Simulating shock...")
    time.sleep(2)
    os.system("mpirun -n {} python3 pic.py --conf shock.ini".format(cores))
    print("------------------------------")
    print("Simulation complete")
    time.sleep(2)
    if (choice == 'Y' or choice == 'y'):
        print("------------------------------")
        print("Running analytical tools...")
        print("Running Particle spectra...")
        time.sleep(2)
        os.system("python3 prtcl_spec.py --conf shock.ini")
        print("Particle spectra complete...")
        time.sleep(2)
        print("Running Particle Path...")
        time.sleep(2)
        os.system("python3 prtcl_path.py --conf shock.ini")
        print("Particle path complete")
        time.sleep(2)
        print("Running Uniplot...")
        time.sleep(2)
        os.system("python3 UniPlot.py --conf shock.ini")
        print("Uniplot complete...")
        time.sleep(2)
        print("Running Shock velocity...")
        time.sleep(2)
        os.system("python3 RH_shk.py --conf shock.ini")
        print("Shock velocity complete")
        print("Analytic tools complete...")
        print("------------------------------")
        print("Running ffmpeg...")
        time.sleep(2)
        os.system("ffmpeg -r 20 -f image2 -pattern_type glob -i \"{}/shock*?png\" -c:v huffyuv shock.avi".format(conf.outdir))
        print("Video created...")
    else:
        print("Skipping analytical tools...")
    print("Program complete, exiting...")
    time.sleep(2)
