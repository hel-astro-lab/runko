# Physical Test Suite
--------------------------------------------------

Simple physical plasma tests are presented here. 

## Weibel instability

Counter-flowing relativistic pair plasma is prone to a current filamentation instability. 
In the case that both beams are symmetric this simplifies to a completely transverse growing modes with simple analytical growth rates.

In order perform a PIC simulation with Weibel setup, execute the code with 4 processors (or less) with

```
mpirun -n 4 python3 pic.py --conf weibel/gam10.ini
```

You can then visualize the resulting simulation files with

```
python3 plot2d_var.py --conf weibel/gam10.ini
```

Growth rate can also be compared to analytical estimate. 
First run `gam3.ini`, `gam10.ini`, and `gam100.ini` setups. 
Then, inside `weibel/` directory, execute the analysis script `python3 plot_all_weibel.py`.



















