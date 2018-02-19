import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os


#--------------------------------------------------
#set up figure
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=7)
plt.rc('ytick', labelsize=7)
plt.rc('axes', labelsize=7)


fig = plt.figure(figsize=(3.54, 3.54)) #single column fig
#fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
gs = plt.GridSpec(1, 1, wspace=0.0)

ax = plt.subplot(gs[0, 0])
ax.minorticks_on()

ax.set_xlabel(r'$k_x$')
ax.set_ylabel(r'$t$')

#--------------------------------------------------
# read simulation data from file
f = h5py.File('out/run.hdf5','r')




#Read field
ex = f['fields/Ex']
#ex = ex[150:, :] #skim off some warm-up phase
#ex = ex[:, 8:340]
print "Ex shape:", np.shape(ex)



#Read simulation values
dt = f['params'].attrs['dt']
dx = f['params'].attrs['dx']
print dt, dx


#--------------------------------------------------
# filter data with hamming windowing
A = np.transpose( ex )
(lines, cols) = np.shape(A)
window = np.hamming(lines).reshape(lines, 1)
A *= window

Ny, Nx = np.shape(A)


#--------------------------------------------------
# Fourier transform in space 
if Nx % 2 == 0:
    nxx = Nx/2 + 1
else:
    nxx = Nx/2

Fourier = np.zeros((nxx, Ny))
for j in range(Ny):
    #print j
    Fourier[:,j] = np.fft.rfft(A[j,:])


(nx, ny) = np.shape(Fourier)
print "shape after transform:", np.shape(Fourier)


# spatial wave vector x component k_x (only half is considered due to Nqyust frequency cut)
#dk = 2.0*np.pi/(nx * dx)
dk = 1.0/(nx*dx)
k = np.arange(nx)*dk
#print "k:"
#print k
k1 = 1
k2 = nx


#temporal guiding vector 
t = np.arange(ny)*dt
#print "t:"
#print t
t1 = 0
t2 = ny


# Change to spectra by considering |F] 
F = np.log10( np.abs(Fourier) ).transpose()
print "min/max:", np.min(F), np.max(F)



K, T = np.meshgrid(k, t)

print "nx=",nx
print "ny=",ny
print np.shape(K)
print np.shape(T)

#slow (but more flexible) pcolormesh that takes guiding grid
#im = ax.pcolormesh(K[w1:w2, k1:k2], W[w1:w2, k1:k2], F[w1:w2, k1:k2],
#            cmap='plasma', 
#            #vmin=ff.min(),
#            #vmax=ff.max(),
#            )


#faster plotting with imshow
im = ax.imshow(F[t1:t2, k1:k2],
        extent=[k[k1], k[k2-1], t[t1], t[t2-1]],
        origin='lower',
        aspect='auto',
        interpolation='nearest',
        cmap='plasma',
        #vmin=np.min(F),
        #vmax=np.max(F)
        vmin= -8.0,
        vmax= -2.0
        )

cax = fig.add_axes([0.12, 0.86, 0.86, 0.03]) 
plt.colorbar(im, cax=cax, orientation='horizontal', label=r'log$_{10} | \mathcal{F}( E_x ) |$')
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')

#ax.set_xlim(0.0, 6.0)
#ax.set_ylim(0.0, 6.0)

##################################################
#analytical relations

wp = 1.0  #plasma frequency 
vth = 1.0 #thermal velocity 


#numerical propagation speed
def w_num(k):
    v= dx/dt
    return k*v
#ax.plot(k, w_num(k), color='black', linestyle='dotted', linewidth=0.8)


#electron plasma frequency
def w_ep(k, wp):
    return np.ones(len(k))*wp

#ax.plot(k, w_ep(k, wp), color='black', linestyle='dotted', linewidth=0.8)

#warm plasma (Langmuir waves)
def w_warmPlasma(k, wp, vth):
    wps = np.ones(len(k))*wp #*2.0
    vths = k*vth #/2.0
    return np.sqrt( wps**2 + 3.0*vths**2 )

#ax.plot(k, w_warmPlasma(k, wp, vth),
#        color='black', 
#        linestyle='dashed', 
#        linewidth=0.8)



plt.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('wt.pdf')
