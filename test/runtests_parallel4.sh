#!/bin/bash

jflags=$@
err=0
start_dir=`pwd`

cd ./advection
  mpirun -np 4 julia $jflags ./runtests_parallel4.jl
  tmp=$?
  err=$((err + tmp))

cd $start_dir

cd ./euler
  mpirun -np 4 julia $jflags ./runtests_parallel4.jl
  tmp=$?
  err=$((err + tmp))
cd $start_dir

echo $err

exit $err
