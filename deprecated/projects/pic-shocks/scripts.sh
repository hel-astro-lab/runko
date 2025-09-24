#!/usr/bin/env bash
 
##!/bin/zsh
 

declare -a conf_arr=(
#"1dsig3.ini"
#"shock_x256m128_p32np4c10_s3d1e-4_v5/1dsig3.ini"
#"shock_x512m1024_p32np4c25_s3d1e-4_v6/1dsig3.ini"
#"shock_x512m1024_p32np4c25_s3d1e-4_v7/1dsig3.ini"
#"shock_x1024m60_p2np8c25_s3d1e-4_v6/2dsig3.ini"
#"shock_x1024m60_p2np8c25_s3d1e-4_v8/2dsig3.ini"
#"3d-shock_x2048m60_p3np16c25_s3d1e-4_v14/3dsig3.ini"
#"3d-shock_x2048m60_p3np16c25_s3d1e-4g5_v14/3dsig3.ini"
"3dsig1.ini"
)


# one-time scripts that read all the time slices in one go
declare -a single_scripts_arr=(
"plot_dens.py"
"plot_jump_conds.py"    
"plot_upstream_ene.py"    
#"plot_pspec.py"
##"plot_ene.py"    
)


# lap-specific scripts that process only one time slice at a time
declare -a lap_scripts_arr=(
"plot_win_2d_shock.py" 
#"plot3d_pyvista_shock.py --var dens" 
#"plot3d_pyvista_shock.py --var cur" 
#"plot2d.py" 
#
#"plot_fld_2d_panel.py --var rho" 
#"plot_fld_2d_panel.py --var jx" 
#"plot_fld_2d_panel.py --var jy" 
#"plot_fld_2d_panel.py --var ey" 
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

for lap in {0..300000..1000}
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

