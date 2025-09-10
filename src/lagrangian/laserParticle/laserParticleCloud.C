/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "laserParticleCloud.H"
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laserParticleCloud, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laserParticleCloud::laserParticleCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    bool readFields
)
:
    Cloud<laserParticle>(mesh, cloudName, false),
    mesh_(mesh),
    laserProperties_
    (
        IOobject
        (
            "laserProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    centre_(laserProperties_.lookupOrDefault("centre", average(mesh_.C().primitiveField()))),
    radius_(laserProperties_.lookupOrDefault("radius", 1.0)),
    direction_(laserProperties_.lookupOrDefault("direction", vector(0.0, -1.0, 0.0))),
    nRays_(laserProperties_.lookupOrDefault("nRays", 0)),
    power_(laserProperties_.lookupOrDefault("power", 1.0))
{
    if (readFields)
    {
        laserParticle::readFields(*this);
    }

    // normalize direction
    direction_ = normalised(direction_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::laserParticleCloud::move()
{
    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");

    interpolationCellPoint<scalar> rhoInterp(rho);

    laserParticle::trackingData td(*this, rhoInterp);

    Cloud<laserParticle>::move(*this, td, mesh_.time().deltaTValue());
}

Foam::scalar Foam::laserParticleCloud::regenerate()
{
    Random rng(12345);

    const meshSearch searchEngine(mesh_);

    // construct basis
    const vector t1 = normalised(perpendicular(direction_));
    const vector t2 = normalised(t1 ^ direction_);
    
    for (int i=0; i< nRays_; i++)
    {
        // sample position
        const scalar r(radius_ * rng.scalar01());
        const scalar w(2*Foam::constant::mathematical::pi * rng.scalar01());
        const vector d(r*Foam::sin(w)*t1 + r*Foam::cos(w)*t2 + centre_);
        
        // distributed power 
        const scalar p(power_ / nRays_);

        // add particle
        if (searchEngine.findCell(d) != -1)
        {
            addParticle(new laserParticle(mesh_, d, p, direction_));
        }
    }
    Info<< "number of initial particles is " << size() << endl;
    return size();
}


// ************************************************************************* //
