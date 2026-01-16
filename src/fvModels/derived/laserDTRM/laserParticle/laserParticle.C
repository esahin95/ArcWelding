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

        // Current position needed for backtracking
        vector p0(position());
        scalar a0 = td.alpha1Interp().interpolate(coordinates(), currentTetIndices());

        // Track to face
        if (mesh().time().writeTime())
        {
            td.append(p0, trackIndex_, d_);
        }
        trackToFace(U_, 0.0);

        // Check if interface was crossed
        const scalar alimit(0.5);
        scalar a1 = td.alpha1Interp().interpolate(coordinates(), currentTetIndices());
        if (a1 < alimit && a0 > alimit)
        {                        
            // Bisection search to backtrack interface position
            vector p1(position());
            scalar am(a1);
            vector di;
            label it;
            for (it=0; it<10; it++)
            {
                if (mag(am - alimit) < 1e-2) 
                {
                    break;
                }

                // Displacement to midpoint
                di = 0.5 * (p0 + p1) - position();

                // Track to midpoint
                trackToFace(di, 0.0);
                am = td.alpha1Interp().interpolate(coordinates(), currentTetIndices());

                // Update brackets 
                if (am > alimit)
                {
                    p0 = position();
                }
                else 
                {
                    p1 = position();
                }
            }

            // Reflection of ray direction
            const vector nHatc = normalised(td.nHatInterp().interpolate(coordinates(), currentTetIndices()));
            td.addToPower(a_*d_, cell());
            U_ -= 2.0 * (U_ & nHatc) * nHatc;
            reflections_++;
            
            // Absorption of ray power
            d_ -= a_ * d_;

            // Delete if out of power or maximum reflections reached
            if (d_ / d0_ < 0.05 || reflections_ >= maxReflections_)
            {
                td.addToPower(d_, cell());
                td.keepParticle = false;                    
            }

            // Track back to face
            if (it > 0)
            {
                if (mesh().time().writeTime())
                {
                    td.append(position(), trackIndex_, d_);
                }
                trackToFace(U_, 0.0);
            }
        }

        // Change owner or hit patch
        trackToAndHitFace(U_, 0.0, cloud, td);
        if (mesh().time().writeTime())
        {
            td.append(position(), trackIndex_, d_);
        }
    }

    return td.keepParticle;
}


void Foam::laserParticle::hitWallPatch(Cloud<laserParticle>& cloud, trackingData& td)
{
    td.keepParticle = false;
}

void Foam::laserParticle::transformProperties(const transformer& transform)
{
    particle::transformProperties(transform);
    U_ = transform.transform(U_);
}


// ************************************************************************* //
