#Import Numpy, Matplotlib and csv reader
import numpy as np
import matplotlib.pyplot as plt
import csv

#Import Savitzky–Golay filter and command line input
from scipy.signal import savgol_filter
from parser import parse_input

#Function to find lines of file - from online
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#Function to find the nearest value in array compared to a given value. Outputs
#element of nearest value - from online
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#Main function
if __name__ == "__main__":
    conf, fdir, args = parse_input()
    frames = int(conf.Nt/conf.interval) #Number of simulation frames
    x_shock = [] #Lists for position of shock and respective lap no.
    t_laps  = []
    
    for i in range(frames+1):
        #Read 1D density values from file
        rhox = []
        rhoy = []
        slap = str(i*conf.interval).rjust(4, '0')
        file = "{}/rho1_{}.csv".format(conf.outdir, slap)

        with open(file.format(slap), "r") as f:
            length = file_len(file)
            reader = csv.reader(f)
            bob = list(reader)
            for l in range(length):
                rhox.append(float(bob[l][0]))
                rhoy.append(float(bob[l][1]))
        f.close()
        
        #Omit First 5 unphysical values
        rhox = rhox[5:]
        rhoy = rhoy[5:]
        
        #Set upstream density to unity
        omp = conf.cfl/conf.c_omp #Plasma frequency
        qe = (omp**2.*conf.gamma)/((conf.ppc*.5)*(1.+abs(conf.me/conf.mi)))
        for conv in range(len(rhoy)):
            rhoy[conv] = rhoy[conv]/(conf.ppc*qe)

        #Apply Savitzky-Golay filter to data
        rhoy_f = savgol_filter(rhoy, 51, 2) #21,2 ... 51,1
        CR = np.amax(rhoy_f) #Compression ratio set as maximum of filtered density

        xs_index = find_nearest(rhoy_f, CR/2) #Midpoint of shock found
        
        #Data for lap appended
        x_shock.append(rhox[xs_index])
        t_laps.append(i*conf.interval)

        print("Frame {} appended".format(i))
    
    #First 15 frames ommitted - Remove unphysical data points
    x_shock = x_shock[15:]
    t_laps = t_laps[15:]
    
    #Lists converted to numpy arrays
    t_laps = np.array(t_laps)
    x_shock = np.array(x_shock)

    #Unit conversion
    t_omp = t_laps*(conf.cfl/conf.c_omp)
    x_shock_sd = x_shock / conf.c_omp

    #Fit straight line to data
    line = np.polyfit(t_omp, x_shock_sd, 1)
    print("The gradient is: {} ".format(line[0]))
    print("The y-intercept is: {}".format(line[1]))
    
    #Plot data and fit
    plt.plot(t_omp, x_shock_sd, '.')
    y = []
    for t in t_omp:
        thing = line[0]*t + line[1]
        y.append(thing)
    
    #Output data and Beta-shock value
    print("beta_shock: {}".format(line[0]*np.sqrt(conf.gamma)))
    plt.plot(t_omp, y, '-')
    plt.show()
