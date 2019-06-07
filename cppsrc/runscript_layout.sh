#!/bin/bash
# Run a bunch of jobs in serial order
# This is for the aggregation/orientation model
# ABSOULUTE PATH TO EXECUTABLE
executable=/u/83/tuikkal1/unix/csm/pipeaggrot/cppsrc/run_control_vol_acorrected.exe
# INPUT: executable phi_grid theta_grid diffusion_coeff shear_rate  aspect_ratio dt max_time s1 s2
# classes what_todo visco_0 Np_max
# what_todo::
#-1 = selfcheck
# 0 = new run
# 1 = continue from last state

phi_grid = 20
theta_grid = 20
aspect_ratio = 20
dt = 1e-5
max_time = 1
# direction of simple shear
s1 = 0
s2 = 1
# do not change these
classes = 1
what_todo = 0

# loop1
for diffusion_coeff in  10.0; do
mkdir  name_me_${s1}${s2}_diff_${diffusion_coeff}
cd name_me_${s1}${s2}_diff_${diffusion_coeff}
# loop2
for sr in 0.1 1 10 20 50 100 
do
mkdir sr_${sr}
cd sr_${sr}
$executable $phi_grid $theta_grid $diffusion_coeff $sr $aspect_ratio $dt $max_time
$s1 $s2 $classes $what_todo   > output.dat
grep VISCORAW output.dat | sed s/#VISCORAW// > viscoraw.dat
cd ..
done # loop4
cd ..
done # loop3


# To get the last viscoraw for each configuration e.g. the following oneliner should work
# for folder in sr_*; do echo $folder | sed s/sr_// | tr "\n" " "; tail -1 ${folder}/viscoraw.dat;
# done > last_viscoraws.dat
# Loop over sr folders, print the numerical part of the folder name, print the last line of the folders viscoraw file
# write all output to last_viscoraws.dat file
