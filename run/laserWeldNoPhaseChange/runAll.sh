#!/bin/bash

# clean case directory
foamCleanCase
rm -r 0 >/dev/null 2>&1
cp -r 0_orig 0

# create new mesh
blockMesh > log.blockMesh
checkMesh > log.checkMesh

# initialize fields
setFields > log.setFields

# run solver
decomposePar -cellDist > log.decompose
mpirun -np 4 interThermoFoam -parallel > log.solver
reconstructPar > log.reconstruct

# postprocess
#python3 postProcess.py 

#interThermoFoam > log.solver