import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os
import matplotlib.ticker as ticker


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

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$t$')

#--------------------------------------------------
# read simulation data from file
f = h5py.File('out/run.hdf5','r')



#Read field
ex = f['fields/Ex']
#ex = ex[:, :340]
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
#A *= window

ny, nx = np.shape(A)
print "shape:", np.shape(A)


# configuration space parameter
x = np.arange(nx)*dx
print "x:"
print x
x1 = 1
x2 = nx

#temporal guiding vector 
t = np.arange(ny)*dt
print "t:"
print t
t1 = 0
t2 = ny


# Change to spectra by considering |F] 
#F = A*A #energy
F = A
print "min/max:", np.min(F), np.max(F)


X, T = np.meshgrid(x, t)

print "nx=",nx
print "ny=",ny
print np.shape(X)
print np.shape(T)

#slow (but more flexible) pcolormesh that takes guiding grid
#im = ax.pcolormesh(K[w1:w2, k1:k2], W[w1:w2, k1:k2], F[w1:w2, k1:k2],
#            cmap='plasma', 
#            #vmin=ff.min(),
#            #vmax=ff.max(),
#            )


vminmax = np.maximum( np.abs(np.min(F)), np.abs(np.max(F)) )

#faster plotting with imshow
im = ax.imshow(F[t1:t2, x1:x2],
        extent=[x[x1], x[x2-1], t[t1], t[t2-1]],
        origin='lower',
        aspect='auto',
        interpolation='nearest',
        cmap='RdYlGn',
        vmin=-vminmax,
        vmax= vminmax
        )


def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)



cax = fig.add_axes([0.12, 0.86, 0.86, 0.03]) 
cbar = plt.colorbar(im, cax=cax, orientation='horizontal', label=r'$E_x$', format=ticker.FuncFormatter(fmt))
#plt.colorbar(im, cax=cax, orientation='horizontal', label=r'log$_{10} | \mathcal{F}( E_x ) |$')
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')

cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=0)
cbar.ax.tick_params(labelsize=5)

#ax.set_xlim(0.0, 6.0)
#ax.set_ylim(0.0, 6.0)


plt.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.85, wspace=0.0, hspace=0.0)
plt.savefig('et.pdf')
