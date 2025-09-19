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

\*---------------------------------------------------------------------------*/

#include "laserParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laserParticle, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::laserParticle::move
(
    Cloud<laserParticle>& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    while (td.keepParticle && td.sendToProc == -1 && stepFraction() < 1)
    {
        if (debug)
        {
            Info<< "Time = " << mesh().time().timeName()
                << " trackTime = " << trackTime
                << " stepFraction() = " << stepFraction() << endl;
        }

        // save current data
        td.append(position(), index(), d());

        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*trackTime*U_, f, cloud, td);

        // crossed interface
        const tetIndices tetIs = this->currentTetIndices();
        scalar alpha1c = td.alpha1Interp().interpolate(this->coordinates(), tetIs);
        if (alpha1c > 0.5)
        {
            td.keepParticle = false;
            td.lPower(cell()) += d();
            td.append(position(), index(), d());
        }


    }

    return td.keepParticle;
}


void Foam::laserParticle::hitWallPatch(Cloud<laserParticle>& cloud, trackingData&)
{
    const vector nw = normal();

    const scalar Un = U_ & nw;

    if (Un > 0)
    {
        U_ -= 2.0*Un*nw;
    }
}

void Foam::laserParticle::transformProperties(const transformer& transform)
{
    particle::transformProperties(transform);
    U_ = transform.transform(U_);
}


// ************************************************************************* //
