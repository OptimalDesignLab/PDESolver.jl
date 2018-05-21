#!/bin/bash

# run all parallel tests in separate sessions

#jflags="--code-coverage=user"
jflags=$@

err=0
start_dir=`pwd`

cd ./advection
  mpirun -np 2 julia $jflags ./runtests_parallel.jl
  tmp=$?
  err=$((err + tmp))

  if [[ $err -gt 0 ]];
  then
    exit $err
  fi

  mpirun -np 4 julia $jflags ./runtests_parallel4.jl
  tmp=$?
  err=$((err + tmp))

  if [[ $err -gt 0 ]];
  then
    exit $err
  fi

cd $start_dir

cd ./euler
  mpirun -np 2 julia $jflags ./runtests_parallel.jl
  tmp=$?
  err=$((err + tmp))

  if [[ $err -gt 0 ]];
  then
    exit $err
  fi

  mpirun -np 4 julia $jflags ./runtests_parallel4.jl
  tmp=$?
  err=$((err + tmp))

  if [[ $err -gt 0 ]];
  then
    exit $err
  fi

cd $start_dir

echo $err

exit $err
