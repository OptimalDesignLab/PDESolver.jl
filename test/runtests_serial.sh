#!/bin/bash

# run all serial tests in separate sessions

jj=julia

tags=$@

#jflags=$@  # take flags as command line arguments
start_dir=`pwd`

err=0

cd ./euler
  $jj $jflags ./runtests.jl $tags
  tmp=$?
  err=$((err + tmp))
cd $start_dir

cd ./advection
  $jj $jflags ./runtests.jl $tags
  tmp=$?
  err=$((err + tmp))
cd $start_dir

cd ./simpleODE/
  $jj $jflags ./runtests.jl $tags
  tmp=$?
  err=$((err + tmp))
cd $start_dir

cd ./elliptic
  $jj $jflags ./runtests.jl $tags
  tmp=$?
  err=$((err + tmp))
cd $start_dir

echo " "
echo "runtests_serial.sh retval: $err"

exit $err





