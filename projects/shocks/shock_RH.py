# Import Numpy, Matplotlib and csv reader
import numpy as np
import matplotlib.pyplot as plt
import csv

# Import Savitzky–Golay filter and command line input
from scipy.signal import savgol_filter
from parser import parse_input

# Function to find lines of file - from online
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# Function to find the nearest value in array compared to a given value. Outputs
# element of nearest value - from online
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# Main function
if __name__ == "__main__":
    conf, fdir, args = parse_input()
    omp = conf.cfl / conf.c_omp  # Plasma frequency
    qe = (omp ** 2.0 * conf.gamma) / (conf.ppc)
    frames = int(conf.Nt / conf.interval)  # Number of simulation frames
    x_shock = []  # Lists for position of shock and respective lap no.
    t_laps = []
    CR_upstream_ar = [] #Arrays for finding CR
    CR_downstream_ar = []
    CR_delete_index = []

##############################################################################
    for i in range(15, frames + 1): #First 15 frames unphysical due to shock formation
        # Read 1D density values from file
        rhox = []
        rhoy = []
        slap = str(i * conf.interval).rjust(5, "0")
        file = "{}/rho1_{}.csv".format(conf.outdir, slap)

        with open(file.format(slap), "r") as f:
            length = file_len(file)
            reader = csv.reader(f)
            bob = list(reader)
            for l in range(length):
                rhox.append(float(bob[l][0]))
                rhoy.append(float(bob[l][1]))
        f.close()

        # Set upstream density to unity
        for conv in range(len(rhoy)):
            rhoy[conv] = rhoy[conv] / (conf.ppc * qe)

        # Apply Savitzky-Golay filter to data
        rhoy_f = savgol_filter(rhoy, 51, 3)  # 21,2 ... 51,1
        #print("Data points = {}".format(len(rhoy_f)))
        ###########################

        # This block of code finds the position of the shock, then identifies an upstream position to measure upstream density. Only use if upstream density <> 1
        #ep = 0
        #r = 0
        #while ep != 1:
        #    r += 1
        #    #print("rhox[r] = {} and rhox[0] = {} and rhoy_f[1] = {}".format(rhox[r], rhox[0], rhoy_f[10]))
        #    rec_a = (rhox[r] - rhox[0])*rhoy_f[10]
        #    rhoy_i = rhoy_f[:r]
        #    rhox_i = rhox[:r]
        #    trap_area = np.trapz(rhoy_i, rhox_i)
        #    #print("r is {} and i is {}".format(r,i))
        #    #print("rec_a = {} amd trap_area = {}".format(rec_a, trap_area))
        #    if (rec_a - trap_area) > 2: #2 for Perpendicular
        #        #CR_upstream = rhoy_f[r] #Use where Upstream is not 1
        #        CR_upstream = 1
        #        #print("CR found at: {}".format(r))
        #        ep = 1

        #CR upstream value nominally 1
        CR_upstream = 1
        shock_loc = 2*CR_upstream #True midpoint of shock y-axis

        xs_index = find_nearest(rhoy_f, shock_loc)  # Midpoint and position of shock found

        #This block of code finds the average Downstream density per frame
        CR_downstream = 0.0
        CR_LFS = 20 #20 is good for Perpendicular shocks, 1000 for Parallel
        CR_upper_index = xs_index-CR_LFS #Controls how far left of shock program measures Downstream
        if (CR_upper_index > 10): #Filter out bad frames
            CR_upstream_ar.append(CR_upstream) #This is upstream density, appended only if a downstream density can be found

            x_upper_count = 10 #Starting away from reflecting wall
            while (x_upper_count < CR_upper_index):
                CR_downstream = CR_downstream + rhoy_f[x_upper_count]
                x_upper_count = x_upper_count + 1
                #print("x_upper_count is {}".format(x_upper_count))
                #print("CR_upper_index is {}".format(CR_upper_index))
            #print("upper avg is {}".format(CR_downstream))
            #print("CR_Upper_index is {}".format(CR_upper_index))
            CR_downstream = CR_downstream/(x_upper_count - 10)

            CR_downstream_ar.append(CR_downstream)
        else:
            CR_delete_index.append(i)
            CR_upstream_ar.append(0.0)
            CR_downstream_ar.append(0.0)

        ###########################

        # Data for lap appended
        x_shock.append(rhox[xs_index])
        t_laps.append(i * conf.interval)
        #print(t_laps[i])

        print("Frame {} appended".format(i))
        print("Upstream density is: {}".format(CR_upstream))
        print("x is: {}".format(rhox[xs_index]))

        d = plt.figure()
        ax = d.add_subplot(111)

        ax.set_xlim(0, 800)
        ax.set_ylim(0, 5)

        plt.plot(rhox, rhoy_f)
        plt.axhline(y=shock_loc, color = 'purple') #Horizontal line halfway up shock
        plt.axhline(y=CR_upstream, color = 'yellow') #Measures upstream density

        if (CR_upper_index > 0):
            plt.axhline(y=CR_downstream, color = 'orange') #Average Downstream density

        plt.axvline(x=rhox[xs_index], color = 'green') #Measures x-position of shock
        #plt.axvline(x = rhox[r], color = 'red') #Use for debug if CR upsteam <> 1
        plt.savefig("{}/sav_{}".format(conf.outdir, i))
        plt.close()
###############################################################################

    x_shock = x_shock[:150]
    t_laps = t_laps[:150]
    CR_upstream_ar = CR_upstream_ar[:150]
    CR_downstream_ar = CR_downstream_ar[:150]

    # Lists converted to numpy arrays
    t_laps = np.array(t_laps)
    x_shock = np.array(x_shock)
    CR_upstream_ar = np.array(CR_upstream_ar)
    CR_downstream_ar = np.array(CR_downstream_ar)

    # Remove Zero points
    zero_points = np.where(x_shock < 20.0)
    x_shock = np.delete(x_shock, zero_points[0])
    t_laps = np.delete(t_laps, zero_points[0])
    CR_upstream_ar = np.delete(CR_upstream_ar, zero_points[0])
    CR_downstream_ar = np.delete(CR_downstream_ar, zero_points[0])

    zero_points = np.where(CR_upstream_ar == 0.0)
    x_shock = np.delete(x_shock, zero_points[0])
    t_laps = np.delete(t_laps, zero_points[0])
    CR_upstream_ar = np.delete(CR_upstream_ar, zero_points[0])
    CR_downstream_ar = np.delete(CR_downstream_ar, zero_points[0])

    #Calculate all CR values across simulation
    CR_ar = CR_downstream_ar/CR_upstream_ar

    # Unit conversion
    t_omp = t_laps * conf.cfl/conf.c_omp

    x_shock_sd = x_shock #/ conf.c_omp

    # Fit straight line to data
    line, covmat = np.polyfit(t_omp, x_shock_sd, 1, cov = True)
    grad_uncert = np.sqrt(float(covmat[0,0]))
    inter_uncert = np.sqrt(float(covmat[1,1]))
    print("The gradient is: {} +- {}".format(line[0], grad_uncert))
    print("The y-intercept is: {} +- {}".format(line[1], inter_uncert))

    #Find CR value by averaging later values
    CR_avg = 0.0
    start_CR_count = 1 #Minimum 1
    stop_CR_count = 50 #Varies depending on how long it takes for CR to stabilise, check CR_RH.png to see if this value needs adjusting
    for CR_count in range(start_CR_count, stop_CR_count):
        CR_avg = CR_avg + CR_ar[-CR_count]
        #print(CR_ar[-CR_count])
    CR_avg = CR_avg/(stop_CR_count - start_CR_count)

    # Plot data and fit
    plt.plot(t_omp, x_shock_sd, ".")
    y = []
    for t in t_omp:
        thing = line[0] * t + line[1]
        y.append(thing)

    # Output data, Beta-shock and CR value
    print("beta_shock: {}".format(line[0]))
    plt.plot(t_omp, y, "-")
    plt.xlabel("t_omp")
    plt.ylabel("x_shock")
    plt.savefig("{}/shock_RH.png".format(conf.outdir))

    plt.clf()

    print("CR: {}".format(CR_avg))
    plt.xlabel("t_omp")
    plt.ylabel("Compression Ratio")
    plt.plot(t_omp, CR_ar, 'bo')
    plt.axhline(y = CR_avg, color = 'red')
    plt.savefig("{}/CR_RH.png".format(conf.outdir))
