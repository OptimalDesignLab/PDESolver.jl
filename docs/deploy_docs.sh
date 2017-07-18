#!/bin/bash

julia ./make.jl
#mkdocs build
mkdocs gh-deploy
