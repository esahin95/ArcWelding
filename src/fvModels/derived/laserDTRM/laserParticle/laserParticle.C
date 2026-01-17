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

    while (td.keepParticle && td.sendToProc == -1)
    {
        // Current position
        vector p0(position());
        scalar a0 = td.alpha1Interp().interpolate(coordinates(), currentTetIndices());
        if (mesh().time().writeTime())
        {
            td.append(p0, trackIndex_, d_);
        }

        // Track to next face
        trackToAndHitFace(U_, 0.0, cloud, td);

        // New position
        vector p1(position());
        scalar a1 = td.alpha1Interp().interpolate(coordinates(), currentTetIndices());
        if (mesh().time().writeTime())
        {
            td.append(p1, trackIndex_, d_);
        }

        // Check if interface was crossed
        if (a1 < 0.5 && a0 > 0.5)
        {                                    
            // Check Reflection
            const vector nHatc = normalised(td.nHatInterp().interpolate(coordinates(), currentTetIndices()));
            scalar un(U_ & nHatc);
            if (un < 0.0)
            {
                // Reflection of ray direction
                U_ -= 2.0 * un * nHatc;
                reflections_++;
            
                // Absorption of ray power
                td.addToPower(a_ * d_, cell());
                d_ -= a_ * d_;
            }

            // Delete if out of power or maximum reflections reached
            if (d_ / d0_ < 0.05 || reflections_ >= maxReflections_)
            {
                td.addToPower(d_, cell());
                td.keepParticle = false;                    
            }
        }
    }

    return td.keepParticle;
}

void Foam::laserParticle::hitWallPatch(Cloud<laserParticle>& cloud, trackingData& td)
{
    td.keepParticle = false;
}

// ************************************************************************* //