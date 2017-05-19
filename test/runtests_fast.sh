#!/bin/bash

# run all parallel tests

jj=julia

jflags=$@  # take flags as command line arguments
if [ $# == 0 ]; then
  tags="tag_shorttest"
else
  tags=$@
fi

start_dir=`pwd`

err=0

# serial tests, all in a single session
$jj ./runtests.jl $tags
tmp=$?
err=$((err + tmp))

# 2 processor tests, all in a single session
mpirun -np 2 $jj ./runtests_parallel2.jl $tags
tmp=$?
err=$((err + tmp))

# 4 processor tests, all in a single session
mpirun -np 4 $jj ./runtests_parallel4.jl $tags
tmp=$?
err=$((err + tmp))

echo " "
echo "runtests_fast.sh retval: $err"

exit $err





