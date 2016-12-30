#!/bin/bash

# run all parallel tests

jj=julia

jflags=$@  # take flags as command line arguments
start_dir=`pwd`

err=0

# serial tests, all in a single session
$jj $jflags ./runtests.jl
tmp=$?
err=$((err + tmp))

# 2 processor tests, all in a single session
mpirun -np 2 $jj $jflags ./runtests_parallel2.jl
tmp=$?
err=$((err + tmp))

# 4 processor tests, all in a single session
mpirun -np 4 $jj $jflags ./runtests_parallel4.jl
tmp=$?
err=$((err + tmp))

echo " "
echo "runtests_fast.sh retval: $err"

exit $err





