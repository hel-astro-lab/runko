#!/bin/zsh
 
# ! /usr/bin/env bash

declare -a conf_arr=(
"2dsig3.ini"
)


# one-time scripts that read all the time slices in one go
declare -a single_scripts_arr=(
#"plot_dens.py"
#"plot_pspec.py"
##"plot_ene.py"    
##"plot_upstream_ene.py"    
)


# lap-specific scripts that process only one time slice at a time
declare -a lap_scripts_arr=(
"plot_win_2d_shock.py" 
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

for lap in {0..10000..100}
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

