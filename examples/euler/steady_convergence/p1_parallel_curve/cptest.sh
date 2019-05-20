#!/bin/bash


for i in `seq 1 3`;  # loop over directories
do
  mkdir -vp ./m$i  # make new directory 
  cp -vr ./input_vals1.jl ./m$i  # copy prototype input file into directory

  # use sed to change the mesh file name to point to the ith mesh
  sed -i "s/p4_1/p4_$i/g" ./m$i/input_vals*.jl  

done
