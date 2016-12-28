#!/bin/bash

# run all serial and parallel tests in single session

jj=julia

jflags=$@  # take flags as command line arguments
start_dir=`pwd`

err=0

cd ./euler
  $jj $jflags ./runtests.jl
  tmp=$?
  err=$((err + tmp))
cd $start_dir

cd ./advection
  $jj $jflags ./runtests.jl
  tmp=$?
  err=$((err + tmp))
cd $start_dir

cd ./simpleODE/
  $jj $jflags ./runtests.jl
  tmp=$?
  err=$((err + tmp))
cd $start_dir

./runtests_parallel.sh $jflags
tmp=$?
err=$((err + tmp))

echo " "
echo "runtests.sh retval: $err"

exit $err





