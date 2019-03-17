#!/bin/bash

# run as many tests in parallel as possible
# for the serial tests, a file named fout is created in each directory
# containing the stderr and stdout of the tests in that directory
# For parallel tests, the file is named fout2 or fout4, for the 2 process
# and 4 process tests, respectively
# If tests failed, a message is printed and the script exits with the return
# status of the process that failed

jj=julia

tags=$@

start_dir=`pwd`
declare -a pids
names=("navier_stokes" "advection" "elliptic" "simpleODE" "euler")
err=0

cd $start_dir
cd ./navier_stokes
$jj $jflags ./runtests.jl &> fout &
#echo "hello world NS" 2>&1 > fout &
pids[0]=$!


cd $start_dir
cd ./advection
$jj $jflags ./runtests.jl &> fout &
#echo "hello world NS" 2>&1 > fout &
pids[1]=$!


cd $start_dir
cd ./elliptic
$jj $jflags ./runtests.jl &> fout &
#echo "hello world elliptic" 2>&1 > fout &
pids[2]=$!

cd $start_dir
cd ./simpleODE/
$jj $jflags ./runtests.jl &> fout &
#echo "hello world simpleODE" 2>&1 > fout &
pids[3]=$!

cd $start_dir
cd ./euler
$jj $jflags ./runtests.jl &> fout &
#echo "hello world Euler" 2>&1 > fout &
pids[4]=$!



for i in `seq 0 4`;
do
  wait ${pids[$i]}
  err=$?
  if [[ $err -gt 0 ]];
  then
      echo "serial ${names[$i]} tests returned exit status: $err"
      exit $err  # does this kill all remaining jobs?
  else
      echo "serial ${names[$i]} tests pass"
  fi
done


# 2 process tests
names2=("advection" "euler")

cd $start_dir
cd ./advection
mpirun -np 2 julia $jflags ./runtests_parallel2.jl &> fout2 &
pids[0]=$!

cd $start_dir
cd ./euler
mpirun -np 2 julia $jflags ./runtests_parallel2.jl &> fout2 &
pids[1]=$1

for i in `seq 0 1`;
do
  wait ${pids[$i]}
  err=$?
  if [[ $err -gt 0 ]];
  then
      echo "2 process ${names2[$i]} tests returned exit status: $err"
      exit $err  # does this kill all remaining jobs?
  else
      echo "2 process ${names2[$i]} tests pass"
  fi
done


# 4 process tests
cd $start_dir
cd ./advection
mpirun -np 4 julia $jflags ./runtests_parallel4.jl &> fout4
err=$?

if [[ $err -gt 0 ]];
then
  echo "4 process advection tests returned exit status: $err"
  exit $err
else
  echo "4 process advection tests pass"
fi


cd $start_dir
cd ./euler
mpirun -np 4 julia $jflags ./runtests_parallel4.jl &> fout4
err=$?

if [[ $err -gt 0 ]];
then
  echo "4 process euler tests returned exit status: $err"
  exit $err
else
  echo "4 process euler tests pass"
fi







