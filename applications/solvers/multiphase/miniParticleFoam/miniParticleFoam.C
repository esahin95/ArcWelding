/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    miniParticleFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "simpleMatrix.H"
#include "solidParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Info << Foam::Time::controlDictName << endl;
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;

    // use unused variables
    Info<< adjustTimeStep << " , " << maxCo << " , " << maxDeltaT << " , " << checkMeshCourantNo << " , " << moveMeshOuterCorrectors << endl;

    // add cloud 
    solidParticleCloud particles(mesh);

    particle newParticle(mesh, vector(0,0.5,0), label(-1));
    solidParticle newSolidParticle
    (
        mesh, 
        newParticle.coordinates(), 
        newParticle.cell(),
        newParticle.tetFace(),
        newParticle.tetPt(),
        scalar(0.001),
        vector(1.7, 0, 0)
    );
    solidParticle *particlePtr = &newSolidParticle;
    particles.addParticle(particlePtr);
    Info<< "number of particles is " << particles.size() << endl;

    // print matrix coefficients
    while (pimple.run(runTime))
    {
    	runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
	        Info<< "In outer corrector loop" << nl << endl;
            while (pimple.correct())
            {
                Info<< "    In corrector loop" << nl << endl;
                while (pimple.correctNonOrthogonal())
                {
                    Info<< "        In nonorthogonal corrector loop" << nl << endl;
                    Info<< "        Is correctPhi relevant? " << correctPhi << nl << endl;
                }
            }
        }

        // move particles
        particles.move(g);
        Info<< "Particle position is " << newSolidParticle.position() << endl;

        runTime.write();
    }

    #include "writeMatCoeffs.H"

    Info<< "\nSolution: " << TEqn.solve() << endl;
    forAll(T, cellI)
    {
        Info<< T.internalField()[cellI] << " ";
    }
    Info<< endl;
    forAll(T.boundaryField(), patchI)
    {
        Info<< T.boundaryField()[patchI].patch().name() << " ";
        forAll(T.boundaryField()[patchI], faceI)
        {
            Info<< T.boundaryField()[patchI][faceI] << " ";
        }
        Info<< endl;
    }
    
    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
