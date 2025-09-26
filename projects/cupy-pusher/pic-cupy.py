"""
PIC simulation to test speed of a cupy pusher.
"""

import matplotlib # This prevents a segfault when importing runko for some reason
import runko
import numpy as np
import sys
import json
import time

rng = np.random.default_rng(seed=42)

logger = runko.runko_logger()
logger.info(f"Logger started.")

config = runko.Configuration(None)

light_crossing_times = 0.01

config.cfl = 0.45
config.Nx = 1024
config.Ny = 64
config.Nz = 64
config.Nt = int(light_crossing_times * config.Nx / config.cfl)
config.Nparticles = 10**5
config.m0 = 1
config.m1 = config.m0
config.outdir = sys.argv[1]

# Plasma motion in lab frame
Gamma = 10
# Plasma magnetisation in co-moving fluid rest frame
sigma = 0
# Cells-per-skindepth perpendicular to bulk motion direction in both rest and lab frames (c_omp = C / OMega Plasma)
c_omp = 128
# Plasma density per species
ppc = 32
# Plasma temperature in co-moving fluid rest frame in units of electron rest mass
temperature = 0.1

# Dependent variables:
# Plasma motion (beta gamma factor)
v_x = -np.sqrt(1-1/Gamma**2)
u_x = Gamma*v_x
# Plasma density (overall, all species per cell)
oppc = 2*ppc
# Particle charge
# My formula:
config.q0 = -config.cfl / c_omp * np.sqrt(config.m0 / oppc) * Gamma**(3/2)
config.q1 = -config.q0
# Magnetic field strength in the lab frame
B_z = config.cfl * np.sqrt(config.m0 * Gamma * sigma * oppc)
# Electric field strength in the lab frame
E_y = v_x*B_z

# End of problem specific configuration
logger.info(f"Starting simulation with {ppc*config.Nx*config.Ny*config.Nz} particles for {ppc*config.Nt} laps..")
logger.info(f"Boxsize is {ppc*config.Nx}x{config.Ny}x{config.Nz}.")

# Configure fields
def get_fields(t, pos):

    E = np.array([[0,E_y,0]])
    B = np.array([[0,0,B_z]])
    
    E = np.repeat(E, len(pos), axis=0)
    B = np.repeat(B, len(pos), axis=0)

    return E, B


logger.info(f"Generating partices")

# Configure particles
def pgen(N):
    # Position
    pos = np.array([rng.random(N)*config.Nx,
                    rng.random(N)*config.Ny,
                    rng.random(N)*config.Nz
    ]).T
    pos = np.concatenate([pos, pos])
    # Velocity
    vel = runko.sample_boosted_juttner_synge(2*N, temperature, gen=rng, Gamma=Gamma, direction='-x')
    vel = np.array(vel).T
    # Mass
    mas = np.concatenate([np.ones(N)*config.m0, np.ones(N)*config.m1])
    # Charge
    q = np.concatenate([np.ones(N)*config.q0, np.ones(N)*config.q1])
    return pos, vel, mas, q

pos, vel, mas, q = pgen(config.Nparticles)

def boris_pusher(pos, vel, m, q, E, B):
    # Charge to mass ratio and c
    c = config.cfl
    qm = q/m
    # Velocity
    vel = vel * c;
    # read particle-specific fields
    e0 = 0.5*qm[:,np.newaxis]*E
    b0 = 0.5*qm[:,np.newaxis]*B
    # first half electric acceleration
    u0 = vel + e0
    # first half magnetic rotation
    ginv = c/np.sqrt(c**2 + np.sum(vel**2,axis=1))
    b0 *= ginv[:,np.newaxis]
    f = 2.0/(1.0 + np.sum(b0,axis=1));
    u1 = (u0 + np.cross(u0, b0))*f[:,np.newaxis]
    # second half of magnetic rotation & electric acceleration
    # Why is f missing here, what is f
    u0 = u0 + np.cross(u1, b0) + e0
    # --------------------------------------------------
    # normalized 4-velocity advance
    vel = u0/c

    # position advance; 
    pos += vel*ginv[:,np.newaxis]*c

    # done
    return pos, vel

def pic_simulation_step(lap):
    t = time.time()
    global pos, vel
    E, B = get_fields(lap*config.cfl, pos)
    pos, vel = boris_pusher(pos, vel, mas, q, E, B)
    t = time.time() - t
    logger.info(f"Lap {lap}/{config.Nt} done in {t:.2f}s, pushed {len(pos)/t:.3e} particles per second.")

t0 = time.time()
for lap in range(config.Nt):
    pic_simulation_step(lap)



# On completion of simulation, dump some of the relevant simulation parameters to a json file
param_dict = config.__dict__
del param_dict['_config_path'] # Unnecessary parameter
# Add physical parameters
param_dict.update({
    'Gamma': Gamma,
    'sigma': sigma,
    'c_omp': c_omp,
    'ppc':   ppc,
    'temperature': temperature,
    'v_x':   v_x,
    'u_x':   u_x,
    'oppc':  oppc,
    'B_z':   B_z,
    'E_y':   E_y,
    'runtime': time.time()-t0,
    'particles-per-sec': len(pos)/(time.time()-t0)
})
with open(config.outdir+'/params.json', 'w') as file:
    file.write(json.dumps(param_dict))