#!/bin/bash


for i in `seq 1 9`;  # loop over directories
do
  mkdir -vp ./m$i  # make new directory 
  cp -vr ./input_vals1.jl ./m$i  # copy prototype input file into directory

  # use sed to change the mesh file name to point to the ith mesh
  sed -i "s/s_1/s_$i/g" ./m$i/input_vals*.jl  

done
