#!/bin/bash

# This script must be run from the /docs directory
# It will build the documentation using Documenter and then deploy it to
# github 

start_dir=`pwd`
julia ./make.jl

if [ $? != 0 ]; then
  echo "Error: Building docs failed, aborting ..."
  exit 1
fi

commit_msg=`git log -n 1 --abbrev-commit --oneline`

# this incantation should work on OS X and Linux
tdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tdir'`
echo "using temporary directory = $tdir"

# download another copy of the repo to a temporary directory
cd $tdir

if [ $? != 0 ]; then
  echo "Error: could not cd into temporary directory"
  exit 1
fi

git clone https://github.com/OptimalDesignLab/PDESolver.jl.git
cd ./PDESolver.jl
git checkout gh-pages
rm -r ./*  # remove everything

# copy docs to temporary repo
pres_dir=`pwd`
echo "copying files from $start_dir/build to $pres_dir"
cp -r $start_dir/build/* $pres_dir
git add ./*
git commit -am "building docs from: $commit_msg"

# deploy
echo "deploying to github"
git push origin gh-pages --force

cd $start_dir

rm -rf $tdir  # clean up temporary directory
