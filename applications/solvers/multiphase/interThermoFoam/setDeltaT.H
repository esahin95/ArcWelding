/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    // thermal diffusion number
    const volScalarField::Internal DiNumfv
    (
        fvc::surfaceSum
        (
            mesh.magSf()
          * fvc::interpolate(kappa)
          * mesh.surfaceInterpolation::deltaCoeffs()
        )()() / (mesh.V() * rho * Cp) * runTime.deltaT()
    );
    DiNum = gMax(DiNumfv);
    Info<< "Diffusion number maximum: " << DiNum << endl;
    
    scalar deltaT = min
    (
        min(maxCo/(CoNum + small), maxAlphaCo/(alphaCoNum + small)),
        maxDi/(DiNum + small)
    ) * runTime.deltaTValue();
    deltaT = min(deltaT, fvModels.maxDeltaT());
    deltaT = min
    (
        min(deltaT, runTime.deltaTValue() + 0.1*deltaT),
        1.2*runTime.deltaTValue()
    );
    runTime.setDeltaT(min(deltaT, maxDeltaT));

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
