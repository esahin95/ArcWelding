#!/bin/bash

# clean case directory
foamCleanCase
rm -r 0 >/dev/null 2>&1
cp -r 0_orig 0

# create new mesh
blockMesh | tee log.blockMesh
createPatch -overwrite 
topoSet | tee log.topoSet
#refineMesh -overwrite | tee log.refineMesh
checkMesh | tee log.checkMesh

# initialize fields
setFields

# run solver
decomposePar -cellDist | tee log.decompose
mpirun -np 4 interFoam -parallel > log.solver
#for i in $(seq 0 3); do touch processor$i/fallingDrop.OpenFOAM; done

reconstructPar | tee log.reconstruct

echo "################################"
tail log.solver
