import numpy as np
import matplotlib.pyplot as plt

import h5py as h5
import sys, os


def lap2time(lap, conf):
    L = conf.Nx*conf.NxMesh
    l0 = L/conf.max_mode
    t0 = l0/conf.cfl

    #print("lap: {} t: {} t0: {} L: {} l_0: {}".format(lap,lap/t0, t0, L, l0))
    return lap/t0


def find_min(xx, th):
    i1 = 0
    for i,x in enumerate(xx):
        if x >= th:
            i1 = i
            break
    return i1



def create_pmom_spectra(fname, args):

    bins = np.logspace(
            np.log10(args['xmin']),
            np.log10(args['xmax']),
            args['nbins'])

    #hist  = np.zeros(args['nbins']-1)
    #histx = np.zeros(args['nbins']-1)
    #histy = np.zeros(args['nbins']-1)
    #histz = np.zeros(args['nbins']-1)

    f5F = h5.File(fname,'r')
    nx = f5F['Nx'][()]
    ny = f5F['Ny'][()]
    nz = f5F['Nz'][()]

    #print("reshaping 1D array into multiD with ({} {} {}) {}".format(nx,ny,nz, fname))

    ux   = np.reshape( f5F['vx'][:]   ,(nx,-1)) 
    uy   = np.reshape( f5F['vy'][:]   ,(nx,-1)) 
    uz   = np.reshape( f5F['vz'][:]   ,(nx,-1)) 
    wgt  = np.reshape( f5F['wgt'][:]  ,(nx,-1)) 

    u    = np.sqrt(ux**2 + uy**2 + uz**2)

    hist , edges = np.histogram(u,  bins=bins, weights=wgt)
    histx, edges = np.histogram(ux, bins=bins, weights=wgt)
    histy, edges = np.histogram(uy, bins=bins, weights=wgt)
    histz, edges = np.histogram(uz, bins=bins, weights=wgt)

    return bins[1:], hist, histx, histy, histz




