#!/bin/bash

# run all serial tests in separate sessions

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

echo " "
echo "runtests_serial.sh retval: $err"

exit $err





