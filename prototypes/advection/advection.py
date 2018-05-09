import numpy as np
import matplotlib.pyplot as plt
import sys, os


plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)

fig = plt.figure(figsize=(3.54, 5.0)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(3, 1, hspace=0.15)


axs = []
axs.append( plt.subplot(gs[0,0]) )
axs.append( plt.subplot(gs[1,0]) )
axs.append( plt.subplot(gs[2,0]) )

for ax in axs:
    ax.minorticks_on()
    #ax.set_xlim((0.0, 100.0))

axs[0].set_ylabel(r'Y')
axs[1].set_ylabel(r'Z')
axs[2].set_ylabel(r'W')






dx = 0.01
Nx = 400
xx = np.array([i*dx for i in range(400)])
yy1 = np.zeros(Nx)
yy2 = np.zeros(Nx)


def initialize(xx, yy):
    for i,xv in enumerate(xx):
        if(xv >= 1.0) and (xv <= 1.01):
            print("...")
            yy[i] = 1.0

    #yy = np.abs(np.sin(2.0*np.pi*xx/4.0))

    return yy

yy1 = initialize(xx,yy1)
yy2 = initialize(xx,yy2)


axs[0].plot(xx, yy2)
axs[1].plot(xx, yy2)

#yy = np.abs(np.sin(2.0*np.pi*xx/4.0))


vel = 1.0
cfl = 1.0

def upwind1(xx, yy):
    flux = np.zeros(Nx)

    #positive going
    if -vel >= 0.0:
        #print("xxx") #wrong
        for j in range(1, Nx-1):
            flux[j] = cfl*vel*yy[j+1]

    if -vel < 0.0:
        for j in range(1, Nx-1):
            flux[j] = cfl*vel*yy[j]

    return flux


def upwind1F(M0, Mp1):
    if -vel >= 0.0: #  <-----
        return cfl*vel*Mp1
    if -vel < 0.0: #    ---->
        return cfl*vel*M0

    
def central2(xx,yy):
    flux = np.zeros(Nx)

    for j in range(1, Nx-1):
        Mp = yy[j+1] + yy[j]
        Mn = yy[j+1] - yy[j]

        flux[j] = ( vel*cfl*0.5*Mp - 0.5*Mn*(vel*cfl)**2.0 )

    return flux



time = 0.0
for t in range(101):
    ################################################## 
    # scheme1

    flux1 = np.zeros(Nx)
    for j in range(1,Nx-1):
        flux1[j] = upwind1F(yy1[j], yy1[j+1])


    for j in range(1,Nx-1):
        if yy1[j] > 0.01: 
            print(j, flux1[j])

        yy1[j]   -= flux1[j] #U_i+1/2 (outflowing)
        yy1[j+1] += flux1[j] #U_i-1/2 (inflowing)

    #flux1 = upwind1(xx,yy1)
    #for j in range(1,Nx-1):
    #    yy1[j  ] -= flux1[j] #U_i+1/2 (outflowing)
    #    yy1[j+1] += flux1[j] #U_i-1/2 (inflowing)



    ################################################## 
    # scheme2
    flux2 = central2(xx,yy2)
    for j in range(1,Nx-1):
        yy2[j]   -= flux2[j] #U_i+1/2 (outflowing)
        yy2[j+1] += flux2[j] #U_i-1/2 (inflowing)


    ################################################## 
    # plotting
    if (t % 10 == 0.0):

        im1 = np.argmax(yy1)
        im2 = np.argmax(yy2)


        if vel > 0.0:
            realvel = (xx[im1]-1.0)/time
        elif vel < 0.0:
            realvel = (1.0 - xx[im1])/time

        print(time, xx[im1], xx[im2], realvel, realvel/vel)


        axs[0].plot(xx, yy1)
        axs[1].plot(xx, yy2)

        #axs[2].plot(xx, flux1)

        slap = str(t).rjust(5, '0')
        fname = 'adv_{}.png'.format(slap)
        plt.savefig(fname)

        #plt.subplots_adjust(left=0.18, bottom=0.07, right=0.98, top=0.95, wspace=0.0, hspace=0.0)
        #plt.savefig('advection'+str(t)+'.png')
    time += cfl*dx
