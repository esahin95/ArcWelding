#!/bin/bash

# clean case directory
foamCleanCase
rm -r 0 >/dev/null 2>&1
cp -r 0_orig 0

# create new mesh
blockMesh > log.blockMesh
checkMesh > log.checkMesh

# run solver
miniParticleFoam | tee log.solver
