#!/bin/bash

# run all tests

jj=julia
jflags=$@  # take flags as command line arguments
start_dir=`pwd`

err=0
$jj $jflags ./runtests.jl
tmp=$?
err=$((err + tmp))

cd ./Advection
$jj $jflags ./runtests_advection.jl
tmp=$?
err=$((err + tmp))

cd $start_dir
./runtests_parallel.sh $jflags
tmp=$?
err=$((err + tmp))

exit $err




