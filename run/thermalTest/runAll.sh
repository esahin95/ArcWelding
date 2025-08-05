#!/bin/bash

# clean case directory
foamCleanCase
rm -r 0 >/dev/null 2>&1
cp -r 0_orig 0

# create new mesh
blockMesh > log.blockMesh
createPatch -overwrite > log.createPatch
#topoSet > log.topoSet
checkMesh > log.checkMesh

# initialize fields
setFields > log.setFields

# run solver
interThermoFoam | tee log.solver
