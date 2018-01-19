#!/bin/bash

# run all parallel tests
# $1 describes the stopping behavior.  If not supplied, the script runs to completion.
# if $1 == 1, the script stopes at the first error
jj=julia
# jj=/usr/local/julia/0.4.7/bin/julia

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
if [[ $err -gt 0 ]];
then
  exit $err
fi

echo "after serial tests, err = $err"

# 2 processor tests, all in a single session
mpirun -np 2 $jj ./runtests_parallel2.jl $tags
tmp=$?
err=$((err + tmp))

if [[ $err -gt 0 ]];
then
  exit $err
fi

echo "after parallel2 tests, err = $err"

# 4 processor tests, all in a single session
mpirun -np 4 $jj ./runtests_parallel4.jl $tags
tmp=$?
err=$((err + tmp))

if [[ $err -gt 0 ]];
then
  exit $err
fi

echo "after parallel4 tests err = $err"


echo " "
echo "runtests_fast.sh retval: $err"

exit $err





