Clusters
########

Runko has been successfully compiled on many different computing clusters. Module loading and environment variable initialization scripts are provided in `archs/` folder.

As a rule of thumb:

- Most important thing to look at cluster performance-wise is to check the memory per core. 
  - PIC simulations are most often memory limited so this number is the actual limiting factor. 
- Second important thing is the inter connector because MPI communications are costly. 
- Thirdly, modern Intel CPUs vectorize mesh calculations efficiently and so all the filtering operations are more efficient on them.


Rusty
=====

Flatiron Institute's jack-of-all-trades at New York. Heterogeneous cluster:

* 360 nodes = 11500 cores
* 240 nodes of 2x14 Broadwell cores = 6700 cores 

  - 512GB/node (18GB/core)
* 120 nodes of 2x20 Skylake cores = 4800 cores 

  - 768GB/node (19.2GB/core)
* 100GB/s Omnipath 


Beskow
======

https://www.pdc.kth.se/hpc-services/computing-systems/beskow

Swedish PDC Cray XC40 machine:

* 11 cabinets = 515 blades = 2,060 compute nodes 
* 67,456 cores in total, 2 x Intel CPUs per node
* 9 cabinets of Xeon E5-2698v3 Haswell 2.3 GHz CPUs (2x16 cores per node), 

  - 64GB/node (2GB/core)

* 2 cabinets of Xeon E5-2695v4 Broadwell 2.1 GHz CPUs (2x18 cores per node)

  - 128/node (3.6GB/core)

* High speed network Cray Aries (Dragonfly topology)


Kebnekaise
==========

https://www.hpc2n.umu.se/resources/hardware/kebnekaise

Swedish HPC2N Umeå heterogeneous Lenovo machine:

* 15 racks = 602 nodes
* 19,288 cores (of which 2,448 cores are KNL-cores)
* 432 nodes of 2x14 Intel Xeon

  - 128GB/node (4.6GB/core)

* 52 nodes of 2x14 Intel Skylake, 

  - 192GB/node (6.7GB/core)

* Infiniband FDR/EDR


Tegner
======

https://www.pdc.kth.se/hpc-services/computing-systems/tegner-1.737437

Heterogeneous post-processing cluster for Beskow. 

.. note::

    Actual simulation part is not working, only the analysis scripts are functional.


