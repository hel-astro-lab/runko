#!/usr/bin/env bash
 
##!/bin/zsh
 

declare -a conf_arr=(
#"shock_3d.ini"
#"che_3d.ini"
"3d_65536x128x128_ppc2_c25_s3_np8_cfl05_t1e-04_gam10_bx0by0bz1_v9_faraday/3d_shock.ini"
"3d_65536x128x128_ppc2_c25_s6_np8_cfl05_t1e-04_gam10_bx0by0bz1_v9_faraday/3d_shock.ini"
)


# one-time scripts that read all the time slices in one go
declare -a single_scripts_arr=(
"plot_energy_timeline.py"
"plot_pspec_time.py"
)


# lap-specific scripts that process only one time slice at a time
declare -a lap_scripts_arr=(
"plot_shock_2d.py --dark" 
#"plot3d_pyvista_shock.py --var dens" 
#"plot3d_pyvista_shock.py --var cur" 
)



#-------------------------------------------------- 
# Single lap scripts
for c in "${conf_arr[@]}"
do
    echo "conf filename is: $c"
    for s in "${single_scripts_arr[@]}"
    do
        echo "script is: $s"
        python3 $s --conf $c
    done
done

#-------------------------------------------------- 
# multilap scripts

for lap in $(seq 0 1000 20000)
do
    echo "lap is $lap"
    for c in "${conf_arr[@]}"
    do
        echo "conf filename is: $c"
        for s in "${lap_scripts_arr[@]}"
        do
            echo "script is: $s"
            python3 $s --conf $c --lap $lap
        done
    done
done

