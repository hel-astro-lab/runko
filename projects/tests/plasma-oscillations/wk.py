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
ax.set_ylabel(r'$\omega$')

#--------------------------------------------------
# read simulation data from file
#f = h5py.File('electrostatic/run_hires.hdf5','r')
#f = h5py.File('electrostatic/run_hires_2nd.hdf5','r')
f = h5py.File('out/run.hdf5','r')


#Read field
ex = f['fields/Ex']
#ex = ex[:, 1500:4000]
print "Ex shape:", np.shape(ex)



#Read simulation values
dt = f['params'].attrs['dt']
dx = f['params'].attrs['dx']


print dt, dx


#-------------------------------------------------- 
# create test data
#A = np.zeros( (100, 1001) )
#print "simulated data shape:", np.shape(A)
#dx = 0.5
#k1 = 5.0
#w1 = 50.0
#t = np.arange(0.0, 1001)*0.05
#x = np.arange(0.0, 100.0)*dx
#
##for i in range(1001):
##    A[:,i] = np.sin(k1*x - w1*t[i]) +  np.random.rand(len(x))*0.001
##    A[:,i] = np.sin(2.0*np.pi*0.15*x + 1.0*t[i]) + np.sin(2.0*np.pi*0.4*x + 3.0*t[i]) + np.random.rand(len(x))*0.001
#
##for j in range(100):
##    A[j,:] += np.sin(- w1*t)
##    A[j,:] += np.sin(k1*x[j] - w1*t)
#
#def wave(x, t):
#    if x == 0.0:
#        return 0.0
#    k1 = 5.0
#
#    #k = 2.0*np.pi/x
#    k = k1/x
#    w = 1.0*k
#    return np.sin( k1*x - w*t )
#
#for i in range(1001):
#    for j in range(100):
#        A[j,i] = wave(x[j], t[i]) + np.random.rand()*0.1
#A = np.transpose(A)

#--------------------------------------------------
# filter data with hamming windowing
A = np.transpose( ex )
(lines, cols) = np.shape(A)
window = np.hamming(lines).reshape(lines, 1)
A *= window


#--------------------------------------------------
# Fourier transform for power spectra

# example 1
#rate = 30.0
#t = np.arange(0, 10, 1/rate)
#x = np.sin(2*np.pi*4*t) + np.sin(2*np.pi*7*t) + np.random.randn(len(t))*0.2
#p = 20*np.log10(np.abs(np.fft.rfft(x)))
#f = np.linspace(0, rate/2, len(p))
#plot(f, p)
#row = A[:,0]
#p = np.log10(np.abs(np.fft.rfft( row )))
#f = np.fft.fftfreq(len(p), d=dx)
#print f, len(f)
#print p, len(p)
#
##f = f[1:len(f)/2]
##p = p[1:len(p)/2]
#ax.plot(f, p)


# example 2
#import scipy.fftpack
#row = A[:,0]
#nx = len(row) # Number of samplepoints
#dx = 1.0 #sample spacing
#x = np.arange(0.0, nx)*dx
##y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
#yf = scipy.fftpack.fft(row)
#xf = np.linspace(0.0, 1.0/(2.0*dx), nx/2)
#ax.plot(xf, 2.0/nx * np.abs(yf[:nx//2]))
#
#plt.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.95, wspace=0.0, hspace=0.0)
#plt.savefig('wk.pdf')
#sys.exit()


#--------------------------------------------------
# Fourier transform in space and time
Fourier = np.fft.rfft2(A)

(ny, nx) = np.shape(Fourier)
print "shape after transform:", np.shape(Fourier)

# spatial wave vector x component k_x (only half is considered due to Nqyust frequency cut)
dk = 2.0*np.pi/(nx * dx)
k = np.arange(nx)*dk
print "k:"
print k
k1 = 1
k2 = nx

# temporal angular frequency \omega
dw = 2.0*np.pi/(ny * dt)
w = np.arange(ny)*dw
print "w:"
print w
w1 = 1
w2 = ny/2


#compute minimum and maximum values
Lx = dx*nx
T  = ny*dt

kmin = 2.0*np.pi/Lx
kmax = 2.0*np.pi/dx

wmin = 2.0*np.pi/T
wmax = 2.0*np.pi/dt

print "x min/max:", dx, Lx
print "t min/max:", dt, T
print "k min/max:", kmin, kmax
print "w min/max:", wmin, wmax

ax.set_xlim(0.0, 40.0)
ax.set_ylim(0.0, 40.0)


# Change to spectra by considering |F] 
F = np.log10( np.abs(Fourier) )
#print "min/max:", np.min(Fourier), np.max(Fourier)
print "min/max:", np.min(F), np.max(F)


K, W = np.meshgrid(k, w)

print "nx=",nx
print "ny=",ny
print np.shape(K)
print np.shape(W)

#slow (but more flexible) pcolormesh that takes guiding grid
#im = ax.pcolormesh(K[w1:w2, k1:k2], W[w1:w2, k1:k2], F[w1:w2, k1:k2],
#            cmap='plasma', 
#            #vmin=ff.min(),
#            #vmax=ff.max(),
#            )

#faster plotting with imshow
im = ax.imshow(F[w1:w2, k1:k2],
        extent=[k[k1], k[k2-1], w[w1], w[w2-1]],
        origin='lower',
        aspect='auto',
        interpolation='nearest',
        cmap='plasma',
        vmin=np.min(F),
        vmax=np.max(F)
        #vmin=-0,
        #vmax=2.0
        )

cax = fig.add_axes([0.12, 0.86, 0.86, 0.03]) 
plt.colorbar(im, cax=cax, orientation='horizontal', label=r'log$_{10} | \mathcal{F}( E_x ) |$')
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')


##################################################
#analytical relations

wp = 1.0 #plasma frequency XXX why x2? 
vth = 0.25 #thermal velocity XXX why /2?

#numerical propagation speed
def w_num(k):
    v= dx/dt
    return k*v
ax.plot(k, w_num(k), color='black', linestyle='dotted', linewidth=0.8)


#electron plasma frequency
def w_ep(k, wp):
    return np.ones(len(k))*wp

ax.plot(k, w_ep(k, wp), color='black', linestyle='dotted', linewidth=0.8)

#warm plasma (Langmuir waves)
def w_warmPlasma(k, wp, vth):
    wps = np.ones(len(k))*wp #*2.0
    vths = k*vth #/2.0
    return np.sqrt( wps**2 + 3.0*vths**2 )

ax.plot(k, w_warmPlasma(k, wp, vth),
        color='black', 
        linestyle='dashed', 
        linewidth=0.8)


plt.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('out/wk.pdf')
