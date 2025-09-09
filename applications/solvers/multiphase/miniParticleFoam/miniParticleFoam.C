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
#include "laserParticleCloud.H"

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

    // unused variables
    Info<< "UNUSED" << adjustTimeStep << maxCo << maxDeltaT << correctPhi << checkMeshCourantNo << moveMeshOuterCorrectors << "END" << endl;

    // move particles
    while (pimple.run(runTime))
    {
    	runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // move particles
        cloud.move();
        forAllConstIter(Cloud<laserParticle>, cloud, iter)
        {
            Info<< "Particle position is " << iter().position() << endl;
        }

        runTime.write();
    }
    
    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
