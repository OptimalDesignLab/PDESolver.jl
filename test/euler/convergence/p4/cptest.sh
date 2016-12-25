#!/bin/bash

#mkdir -v ./m1
#cp -v ./input_vals_vortex3.jl ./m1

for i in `seq 4 5`;
do
  mkdir -vp ./m$i
  cp -vr ./input_vals_vortex3.jl ./m$i
  sed -i "s/vortex3_1/vortex3_$i/g" ./m$i/input_vals_vortex3.jl
#  sed -n "s/vortex2_1/vortex2_$i/gpw ./m$1/input_vals_vortex3.jl" ./m$i/input_vals_vortex3.jl
#  vi ./m$i/input_vals_vortex3.jl
done
