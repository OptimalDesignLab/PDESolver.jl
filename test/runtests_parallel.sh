#!/bin/bash

jflags="--inline=no --code-coverage=user"
#:jflags=""
err=0
start_dir=`pwd`
cd ./Advection
  mpirun -np 2 julia $jflags ./runtests_parallel.jl
  tmp=$?
  sum=$((err + tmp))
cd $start_dir

echo $err

exit $err
