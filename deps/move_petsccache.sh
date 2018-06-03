#!/bin/bash

# used by .travis.yml to move cached files out of the way during the build
# and put them back when it is donea

# $1 = "remove": move files out of the way
# $1 = "restore": move files back to original location

PETSC_DIR1=$HOME/.julia/v0.6/PETSc2/deps/petsc-3.7.6
PETSC_ARCH1=arch-linux2-c-debug
pth=$PETSC_DIR1/$PETSC_ARCH1

# move files away
if [[ $1 == "remove" ]];
then
  
  if file $pth/lib/libpetsc.so ;  # if cache loaded successfully
  then
    export PETSC_DIR=$PETSC_DIR1
    export PETSC_ARCH=$PETSC_ARCH1
    ls $pth/lib
    mv -v $pth/lib $HOME/petsc_lib
  fi

  # do this in either case because the directory gets created even if
  # the cache does not exist
  rm -r $HOME/.julia/v0.6/PETSc2
fi


# put them back
if [[ $1 == "restore" ]];
then

  if file $HOME/petsc_lib/libpetsc.so ;  # if cache loaded successfully
  then
    mkdir -vp $pth
    mv -v $HOME/petsc_lib/ $pth/lib
  fi
fi

