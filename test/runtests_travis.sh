#!/bin/bash

# this script is used to run the tests on Travis CI.
# In particular, it is used to run different tests as part of the Travis
# build matrix

err=0
# NOTE: these if statements must exactly match the `env` section in .travis.yml
if [[ $TEST_ADVECTION == "1" ]];
then
  cd ./advection
  julia ./runtests.jl TAG_SHORTTEST
  err=$(( err + $?))

  mpirun -np 2 julia ./runtests_parallel2.jl TAG_SHORTTEST
  err=$(( err + $?))

  mpirun -np 4 julia ./runtests_parallel4.jl TAG_SHORTTEST
  err=$(( err + $?))

  cd ..
fi

if [[ $TEST_EULER == "1" ]];
then
  cd ./euler
  julia ./runtests.jl TAG_SHORTTEST
  err=$(( err + $?))

  mpirun -np 2 julia ./runtests_parallel2.jl TAG_SHORTTEST
  err=$(( err + $?))

  mpirun -np 4 julia ./runtests_parallel4.jl TAG_SHORTTEST
  err=$(( err + $?))

  # do the viscous tests here because they are short and rely on Euler
  cd ../navier_stokes
  julia ./runtests.jl TAG_SHORTTEST
  err=$(( err + $?))

  cd ..
fi

if [[ $TEST_ELLIPTIC == "1" ]];
then
  cd ./elliptic
  julia ./runtests.jl TAG_SHORTTEST
  err=$(( err + $?))

  #mpirun -np 2 julia ./runtests_parallel2.jl TAG_SHORTTEST
  #err=$(( err + $?))

  #mpirun -np 4 julia ./runtests_parallel4.jl TAG_SHORTTEST
  #err=$(( err + $?))

  cd ..
fi

if [[ $TEST_SIMPLEODE == "1" ]];
then
  cd ./simpleODE
  julia ./runtests.jl TAG_SHORTTEST
  err=$(( err + $?))

  cd ..
fi



exit $err
