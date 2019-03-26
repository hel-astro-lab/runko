#!/bin/sh

# script to install mpi for travis
# source: https://raw.githubusercontent.com/ParBLiSS/kmerind/master/conf/travis-install-mpi.sh

set -e

MPI_IMPL="$1"
#os=`uname`
os=$TRAVIS_OS_NAME

case "$os" in
    linux)
        #sudo apt-get update -q
        case "$MPI_IMPL" in
            mpich3)
                #sudo apt-get install -q gfortran libcr0 default-jdk
                #wget -q http://www.cebacad.net/files/mpich/ubuntu/mpich-3.1/mpich_3.1-1ubuntu_amd64.deb
                #sudo dpkg -i ./mpich_3.1-1ubuntu_amd64.deb
                #rm -f ./mpich_3.1-1ubuntu_amd64.deb
                if [ ! -d "$HOME/local/$MPI_IMPL/bin" ]; then
                    wget http://www.mpich.org/static/downloads/3.1.4/mpich-3.1.4.tar.gz
                    tar -xzf mpich-3.1.4.tar.gz
                    cd mpich-3.1.4
                    ./configure --prefix=$HOME/local/$MPI_IMPL --disable-fortran && make -j 4 && make install && make clean
                    cd ../../
                else
                    echo 'Using cached MPICH 3.1.4 directory';
                fi
                ;;
            openmpi40)
                if [ ! -d "$HOME/local/$MPI_IMPL/bin" ]; then
                    mkdir -p openmpi && cd openmpi
                    wget --no-check-certificate http://www.open-mpi.org/software/ompi/v4.0/downloads/openmpi-4.0.0.tar.bz2
                    tar -xjf openmpi-4.0.0.tar.bz2
                    cd openmpi-4.0.0
                    export OMPI_CC=$COMP_CC
                    export OMPI_CXX=$COMP_CXX
                    ./configure CC=$COMP_CC CXX=$COMP_CXX --prefix=$HOME/local/$MPI_IMPL > /dev/null
                    make -j 4 > /dev/null
                    make install > /dev/null
                    make clean > /dev/null
                    cd ../../

                    export PATH=$PATH:$HOME/local/$MPI_IMPL/bin
                    export PATH=$PATH:$HOME/local/$MPI_IMPL/include
                    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/$MPI_IMPL/lib

                    echo 'Using cached OpenMPI 4.0.0 directory';
                fi
                ;;
            *)
                echo "Unknown MPI implementation: $MPI_IMPL"
                exit 1
                ;;
        esac
        ;;
    osx)
        echo "Installation using homebrew..."
        brew install open-mpi
        ;;
    *)
        echo "Unknown operating system: $os"
        exit 1
        ;;
esac
