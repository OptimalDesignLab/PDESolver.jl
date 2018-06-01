#!/bin/bash

# run all serial and parallel tests in separate sessions

jj=julia

tags=$@

start_dir=`pwd`

err=0

cd ./euler
$jj $jflags ./runtests.jl
tmp=$?
err=$((err + tmp))

if [[ $err -gt 0 ]];
then
    exit $err
fi

cd $start_dir

cd ./advection
$jj $jflags ./runtests.jl
tmp=$?
err=$((err + tmp))
if [[ $err -gt 0 ]];
then
    exit $err
fi

cd $start_dir

cd ./elliptic
$jj $jflags ./runtests.jl
tmp=$?
err=$((err + tmp))
if [[ $err -gt 0 ]];
then
    exit $err
fi

cd $start_dir

cd ./simpleODE/
$jj $jflags ./runtests.jl
tmp=$?
err=$((err + tmp))
if [[ $err -gt 0 ]];
then
    exit $err
fi
cd $start_dir

./runtests_parallel.sh $jflags
tmp=$?
err=$((err + tmp))

echo " "
echo "runtests.sh retval: $err"

exit $err





