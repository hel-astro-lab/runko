#!/bin/bash

# this is an automation script to compile plasma-toolbox 
# and copy the binaries to the parallel file system platform
# accessibe by all nodes.



# git repo
PLASMA_DIR="/home/n/natj/plasma-toolbox"

# pfs 
PFS_PLASMA_DIR="/pfs/nobackup/home/n/natj/plasma-toolbox"


# setup environment
module load foss
module load HDF5
module load Python/2.7.14


# make
#make -j4


# copy binaries to parallel file system
cp ${PLASMA_DIR}/python/pyplasma.so     ${PFS_PLASMA_DIR}/
cp ${PLASMA_DIR}/python/pypic.so        ${PFS_PLASMA_DIR}/
cp ${PLASMA_DIR}/corgi/pycorgi/corgi.so ${PFS_PLASMA_DIR}/

# copy python scripts 
mkdir ${PFS_PLASMA_DIR}/python
mkdir ${PFS_PLASMA_DIR}/tests
mkdir ${PFS_PLASMA_DIR}/projects
mkdir ${PFS_PLASMA_DIR}/projects/tests
mkdir ${PFS_PLASMA_DIR}/pic
mkdir ${PFS_PLASMA_DIR}/jobs

#python modules
#cp ${PLASMA_DIR}/python/*.py ${PFS_PLASMA_DIR}/python/

#unit tests
#cp ${PLASMA_DIR}/tests/*.py ${PFS_PLASMA_DIR}/tests/

#projects/tests
#cp ${PLASMA_DIR}/projects/tests/*.py ${PFS_PLASMA_DIR}/projects/tests/
#cp ${PLASMA_DIR}/projects/tests/*.ini ${PFS_PLASMA_DIR}/projects/tests/

#copy slurm jobs
#cp ${PLASMA_DIR}/jobs/*.sh ${PFS_PLASMA_DIR}/jobs/

# particle-in-cell
#cp ${PLASMA_DIR}/pic/*.py ${PFS_PLASMA_DIR}/pic/
#cp ${PLASMA_DIR}/pic/*.ini ${PFS_PLASMA_DIR}/pic/
