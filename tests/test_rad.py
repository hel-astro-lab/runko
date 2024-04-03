from mpi4py import MPI
import unittest

import sys
import numpy as np

import pycorgi
import pyrunko.qed as pyqed
import pytools

from numpy import sqrt
from math import pi

#import initialize_pic as init_pic
#import initialize_rad as init_rad

#from initialize_rad import bbodySample
#from initialize_rad import rand3Dloc
#from initialize_rad import rand3Dvel


do_plots = True

try:
    import matplotlib.pyplot as plt
except:
    pass


#make tests deterministic by fixing the RNG seed
np.random.seed(0)


# L2 norm of vector
def norm(x):
    return np.sqrt(np.dot(x,x))
    #return np.sqrt( vec[0]**2 + vec[1]**2 + vec[2]**2 )

#unit vector into direction a x b
def uCross(vecA, vecB):
    vecC = np.cross(vecA, vecB)

    nC = norm(vecC) #lenght of new vec
    if nC == 0.0: return 0.0
    return vecC / nC 

# perform pair annihilation "scattering" event  
def _interact_pairann(zmvec, zpvec):
    zm = norm(zmvec) # electron momenta z_-
    zp = norm(zpvec) # positron momenta z_+

    gamm = sqrt(zm**2 + 1) # electron gamma_-
    gamp = sqrt(zp**2 + 1) # positron gamma_+

    omm = zmvec/zm # direction of electron momenta Omega_-
    omp = zpvec/zp # direction of positron momenta Omega_+
    
    zeta = np.dot(omm, omp) #angle between electron and positron momenta
    s0 = gamm + gamp                        #s0 = x + x1 #sqrt(2q) in CoM
    s  = sqrt(zm**2 + zp**2 + 2*zm*zp*zeta) #s  = sqrt(x**2  + x1**2 + 2*x*x1*mu)
    q  = gamp*gamm - zp*zm*zeta + 1         #q  = x*x1*(1-mu)
    #q = s0**2 - s**2

    svec = zp*omp + zm*omm                  #svec = x*om + x1*om1


    # CoM frame variables; x_c
    bc = -svec/s0 #lorentz transform along CoM velocity vector
    gc = s0/np.sqrt(2*q) # lorentz factor of the boost; CoM vel
    uv = bc*gc # com four vel
    v2 = np.dot(bc,bc) # v_c^2
    vc = np.sqrt(v2)

    # boosted variables in CoM frame; x_cm 
    gcm = np.sqrt(q/2) # prtcl energies are equal in this frame; photons also have this energy
    vcm = np.sqrt(1-1/gcm**2) # v_cm

    # angle between b_cm and b_c; electron and CoM 
    y = (1/vc/vcm)*((gamp - gamm)/(gamp + gamm)) 

    # build new coordinate system along svec
    kvec = bc/norm(bc) # z axis along b_c vector
    jvec = uCross(kvec, omm) # y axis to electron direction  # np.cross(omp, omm)/sqrt(1 - zeta**2)
    ivec = uCross(jvec, kvec) # x axis just orthogonal to others # ( (zp + zm*zeta)*omm - (zm + zp*zeta)*omp )/(s*sqrt(1 - zeta**2))
    M = np.array([ ivec, jvec, kvec ])  


    #--------------------------------------------------
    # draw angles
    niter = 0
    ran1, ran2, ran3 = 0,0,0
    while True:

        ran1 = np.random.rand()
        ran2 = np.random.rand()
        ran3 = np.random.rand()

        # angle between k_cm and b_c; photon and CoM
        z = -1 + 2*ran1
        phi   = 2*pi*ran2

        # angle between k_cm and b_cm; photon and electron
        x = y*z + np.sqrt(1-y**2)*np.sqrt(1-z**2)*np.cos(phi)  

        # four product scalar between electron/positron and primary/secondary photon
        z1 = (gcm**2)*(1 - vcm*x) 
        z2 = (gcm**2)*(1 + vcm*x) 

        # differential cross section angle part; F function 
        F = 0.5*( (z1/z2) + (z2/z1) + 2*( (1/z1) + (1/z2) ) - ( (1/z1) + (1/z2) )**2  )
        F *= 1/((1+vcm)*gcm**2) # normalize to [0,1]

        if F > ran3:
            break # accept angles
        if niter > 1e3:
            break
        niter += 1


    # new photon vectors in CoM frame
    sinz = np.sqrt(1-z**2) # sin z
    omrR = np.array([ sinz*np.cos(phi), sinz*np.sin(phi), z ])
    #omr1R = -1*omrR # other photon has same but opposite direction

    # rotate back to lab frame angles
    omr  = np.dot( np.linalg.inv(M), omrR  ) 
    omr1 = -1*omr # other photon has same but opposite direction # np.dot( np.linalg.inv(M), omr1R ) 


    # boost matrix back to lab frame; constructed from b_c vector
    B = np.array([
                [gc,     -uv[0],                -uv[1],                -uv[2],                ],
                [-uv[0], 1+(gc-1)*bc[0]**2/v2,  (gc-1)*bc[1]*bc[0]/v2, (gc-1)*bc[2]*bc[0]/v2, ],
                [-uv[1], (gc-1)*bc[0]*bc[1]/v2, 1+(gc-1)*bc[1]**2/v2,  (gc-1)*bc[2]*bc[1]/v2, ],
                [-uv[2], (gc-1)*bc[0]*bc[2]/v2, (gc-1)*bc[1]*bc[2]/v2, 1+(gc-1)*bc[2]**2/v2,  ] 
                ])

    # four momenta of photons
    xp  = gcm*np.array([1, omr[0],  omr[1],  omr[2]])
    xp1 = gcm*np.array([1, omr1[0], omr1[1], omr1[2]])
         
    # boost 
    xpp  = np.dot(B, xp)
    xpp1 = np.dot(B, xp1)

    x  = 1*xpp[0]  # energy of primary photon
    x1 = 1*xpp1[0] # energy of secondary photon

    om  = np.array([ xpp[1],  xpp[2],  xpp[3], ])/x
    om1 = np.array([ xpp1[1], xpp1[2], xpp1[3],])/x1


    # test energy conservation
    enec = gamm + gamp - (x + x1)
    momc = zmvec + zpvec - (x*om + x1*om1) 

    nom1 = norm(om)
    nom2 = norm(om1)
    nom1i = norm(omm)
    nom2i = norm(omp)

    # FIXME: remove these debug tests
    t1,t2,t3,t4,t5,t6,t7,t8 = False,False,False,False,False,False,False,False
    if np.abs(enec) > 1e-12:    t1 = True
    if x  < 0.0:                t2 = True
    if x1 < 0.0:                t3 = True
    if np.abs(nom1i)-1 > 1e-12: t4 = True
    if np.abs(nom2i)-1 > 1e-12: t5 = True
    if np.abs(nom1)-1 > 1e-12:  t6 = True
    if np.abs(nom2)-1 > 1e-12:  t7 = True
    if np.abs(np.sum(momc)) > 1e-12: t8 = True

    print('random numbers', ran1, ran2, ran3)

    #if t1 or t2 or t3 or t4 or t5 or t6 or t7 or t8:
    if True:
        print('ERROR PAIR-ANN:')
        print(t1,t2,t3,t4,t5,t6,t7,t8)
        print('zm, zp   ', zmvec, zpvec)
        print('gm, gp   ', gamm, gamp)
        print('z,s0,s,q ', zeta, s0, s, q)
        print('x,x1,enec', x, x1, enec)
        print('n(om/om1)', nom1, nom2)
        print('n(omi)   ', nom1i, nom2i)
        print('momc     ', np.sum(momc), momc, momc/np.sum(momc))
        print('th,phi,om', z, phi, om, om1)

    return x,om, x1,om1



# basic Conf file/class for PiC simulation testing
class Conf:

    Nx = 1
    Ny = 1
    Nz = 1

    oneD = False
    twoD = False
    threeD = False

    NxMesh = 10
    NyMesh = 10
    NzMesh = 1

    xmin = 0.0
    xmax = 10.0

    ymin = 0.0
    ymax = 10.0

    zmin = 0.0
    zmax = 10.0

    cfl = 0.45
    c_omp = 10.0
    ppc = 1 #particles per cell
    ppt = 1 #photons per tile

    dx = 1.0
    dy = 1.0
    dz = 1.0

    gamma_e = 0.0
    gamma_i = 0.0

    me = 1
    mi = 1

    Nspecies = 1

    outdir = "out"

    qe = 1.0
    qi = 1.0

    #def __init__(self):
    #    print("initialized...")

    #update bounding box sizes
    #
    # NOTE: NxMesh = 5 grid looks like this:
    #
    # xmin      xmax
    #  |         |
    #  v         v
    #  |_|_|_|_|_
    #  0 1 2 3 4 5
    #
    def update_bbox(self):
        self.xmin = 0.0
        self.xmax = self.Nx*self.NxMesh

        self.ymin = 0.0
        self.ymax = self.Ny*self.NyMesh

        self.zmin = 0.0
        self.zmax = self.Nz*self.NzMesh


class radiation(unittest.TestCase):

    def skip_test_initialization(self): # NOTE: no photon container defined anymore; old code

        #plt.fig = plt.figure(1, figsize=(3,3))
        #plt.rc('font', family='serif', size=12)
        #plt.rc('xtick')
        #plt.rc('ytick')
        #
        #gs = plt.GridSpec(1, 1)
        #
        #axs = []
        #for ai in range(1):
        #    axs.append( plt.subplot(gs[ai]) )

        conf = Conf()
        conf.twoD = True

        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.Nx = 3
        conf.Ny = 3
        conf.Ny = 1
        conf.ppc = 1
        conf.update_bbox()

        kT = 1.0 #black body photon field temperature

        Nprtcls = conf.NxMesh * conf.NyMesh * conf.NzMesh * conf.ppc

        container = pyqed.PhotonContainer()
        container.reserve(Nprtcls)

        weight = 1.0

        ene_ref = np.zeros(Nprtcls)
        wgt_ref = np.zeros(Nprtcls)
        x0_ref  = np.zeros((Nprtcls,3))
        u0_ref  = np.zeros((Nprtcls,3))

        for ip in range(Nprtcls):
            ene = pytools.rad.sample_blackbody(kT)

            x0 = pytools.rad.rand_3D_loc(conf)
            u0 = pytools.rad.rand_3D_vel(1.0)
        
            container.add_particle(x0, u0, weight, ene)
        
            # add also to reference array
            ene_ref[ip]  = ene
            wgt_ref[ip]  = weight
            x0_ref[ip,:] = x0
            u0_ref[ip,:] = u0


        ene  = container.ene()
        wgt  = container.wgt()

        loc0 = container.loc(0)
        loc1 = container.loc(1)
        loc2 = container.loc(2)

        vel0 = container.vel(0)
        vel1 = container.vel(1)
        vel2 = container.vel(2)

        for ip in range(Nprtcls):
            self.assertAlmostEqual( container.ene()[ip],  ene_ref[ip],  places=5)
            self.assertAlmostEqual( container.wgt()[ip],  wgt_ref[ip],  places=5)

            self.assertAlmostEqual( container.loc(0)[ip], x0_ref[ip,0], places=5)
            self.assertAlmostEqual( container.loc(1)[ip], x0_ref[ip,1], places=5)
            self.assertAlmostEqual( container.loc(2)[ip], x0_ref[ip,2], places=5)

            self.assertAlmostEqual( container.vel(0)[ip], u0_ref[ip,0], places=5)
            self.assertAlmostEqual( container.vel(1)[ip], u0_ref[ip,1], places=5)
            self.assertAlmostEqual( container.vel(2)[ip], u0_ref[ip,2], places=5)

    
    # test that tile supports both pic and rad initialization
    def skip_test_initialization(self): # NOTE: no radiation tile defined anymore; old code

        conf = Conf()
        conf.twoD = True

        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.NzMesh = 1
        conf.Nx = 1
        conf.Ny = 1
        conf.Nz = 1
        conf.ppc = 1
        conf.update_bbox()

        grid = pycorgi.twoD.Grid(conf.Nx, conf.Ny, conf.Nz)
        grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

        c = pyqed.twoD.Tile(conf.NxMesh, conf.NyMesh, conf.NzMesh)
        pytools.pic.initialize_tile(c, (0, 0, 0), grid, conf)
        pytools.rad.initialize_tile(c, (0,0,0), grid, conf)

    def skip_test_blackbody(self): # NOTE: old code; no Photon Container anymore

        if do_plots:
            try:
                plt.fig = plt.figure(1, figsize=(3,3))
                plt.rc('font', family='serif', size=12)
                plt.rc('xtick')
                plt.rc('ytick')
                
                gs = plt.GridSpec(1, 1)
                 
                axs = []
                for ai in range(1):
                    axs.append( plt.subplot(gs[ai]) )
            except:
                pass

        conf = Conf()
        conf.twoD = True

        conf.NxMesh = 3
        conf.NyMesh = 3
        conf.ppc = 1
        conf.update_bbox()

        kTbb = 1.5e-2/511.0 # Photon blackbody temperature in units of mc^2; 2e-4 = 0.1 keV

        #Nprtcls = conf.NxMesh * conf.NyMesh * conf.NzMesh * conf.ppc
        Nprtcls = int(1e4)

        container = pyqed.PhotonContainer()
        container.reserve(Nprtcls)

        weight = 1.0
        for ip in range(Nprtcls):
            ene = pytools.rad.sample_blackbody(kTbb)

            x0 = pytools.rad.rand_3D_loc(conf)
            u0 = pytools.rad.rand_3D_vel(1.0)
        
            container.add_particle(x0, u0, weight, ene)
        

        nbins = 100
        emin = 1.0e-3
        emax = 2.0e+3
        #Nph = container.size()
        Nph = Nprtcls
        self.assertEqual(Nprtcls, Nph)

        if do_plots:
            try:
                axs[0].set_xlim(emin, emax)
                #axs[0].set_ylim(1e-4, 1.2e-0)
                axs[0].minorticks_on()
                axs[0].set_xscale('log')
                axs[0].set_yscale('log')
            except:
                pass

        ene = np.array(container.ene())
        phhist, edges = np.histogram(ene*511., np.logspace(np.log10(emin), np.log10(emax), nbins))

        if do_plots:
            try:
                axs[0].bar(edges[:-1], phhist.astype(np.float32)/Nph, width=np.diff(edges), log=True)# number of photons per log energy
            except:
                pass

        prtcl_sum=np.sum(phhist*edges[:-1])/Nph
        #print("Energy of photons: {}".format(prtcl_sum))

        bbrad_sum = np.sum(3.0e12*2.*edges**3/(np.exp(edges/kTbb)-1.0))
        #print("Blackbody energy: {}".format(bbrad_sum))


        try:
            plt.savefig("blackbody.pdf")
        except:
            pass



    def test_pairann(self):

        #conf = Conf()
        #conf.threeD = True
        #conf.NxMesh = 3
        #conf.NyMesh = 3
        #conf.NzMesh = 3
        #conf.Nx = 1
        #conf.Ny = 1
        #conf.Nz = 1
        #conf.ppc = 1
        #conf.update_bbox()

        #grid = pycorgi.threeD.Grid(conf.Nx, conf.Ny, conf.Nz)
        #grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax, conf.zmin, conf.zmax)

        t1 = 'e-'
        t2 = 'e+'
        ux1, uy1, uz1 = 1.0, 1.1, 1.2
        ux2, uy2, uz2 = 1.1, 2.2, 1.3

        intr = pyqed.PairAnn('e-', 'e+')
        intr.get_minmax_ene('e-', 'e+', 1.0)

        #print('before')
        #print(t1, ux1, uy1, uz1)
        #print(t2, ux2, uy2, uz2)

        cs = intr.comp_cross_section( t1, ux1, uy1, uz1, t2, ux2, uy2, uz2 )
        #intr.interact(           t1, ux1, uy1, uz1, t2, ux2, uy2, uz2 )
        t3, ux3, uy3, uz3, t4, ux4, uy4, uz4 = intr.interact(t1, ux1, uy1, uz1, t2, ux2, uy2, uz2 )

        #print('after')
        #print('cs', cs)
        #print(t3, ux3, uy3, uz3)
        #print(t4, ux4, uy4, uz4)

        # python version
        #a0, avec, b0, bvec = _interact_pairann( np.array([ux1, uy1, uz1]), np.array([ux2, uy2, uz2]) )
        #print('py ver')
        #print('a', a0, avec)
        #print('b', b0, bvec)


    def test_comp(self):

        t1 = 'e-'
        t2 = 'ph'
        ux1, uy1, uz1 = 1.0, 1.1, 1.2
        ux2, uy2, uz2 = 0.01, 0.26, 0.17

        intr = pyqed.Compton('e-', 'ph')
        intr.get_minmax_ene('e-', 'ph', 1.0)

        #print('before')
        #print(t1, ux1, uy1, uz1)
        #print(t2, ux2, uy2, uz2)

        cs = intr.comp_cross_section( t1, ux1, uy1, uz1, t2, ux2, uy2, uz2 )
        #intr.interact(           t1, ux1, uy1, uz1, t2, ux2, uy2, uz2 )
        t3, ux3, uy3, uz3, t4, ux4, uy4, uz4 = intr.interact(t1, ux1, uy1, uz1, t2, ux2, uy2, uz2 )

        #print('after')
        #print('cs', cs)
        #print(t3, ux3, uy3, uz3)
        #print(t4, ux4, uy4, uz4)

        # python version
        #a0, avec, b0, bvec = _interact_photann( np.array([ux1, uy1, uz1]), np.array([ux2, uy2, uz2]) )
        #print('py ver')
        #print('a', a0, avec)
        #print('b', b0, bvec)




