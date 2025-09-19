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

Class
    Foam::laserParticleCloud

Description
    A Cloud of solid particles

SourceFiles
    laserParticleCloudI.H
    laserParticleCloud.C
    laserParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef laserParticleCloud_H
#define laserParticleCloud_H

#include "Cloud.H"
#include "laserParticle.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     Class laserParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class laserParticleCloud
:
    public Cloud<laserParticle>
{
    // Private Data

        const fvMesh& mesh_;

        IOdictionary laserProperties_;
        vector centre_;
        scalar radius_;
        vector direction_;
        scalar nRays_;
        scalar power_;

public:

    //- Runtime type information
    TypeName("laserParticleCloud");


    // Constructors

        //- Construct given mesh
        laserParticleCloud
        (
            const fvMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );

        //- Disallow default bitwise copy construction
        laserParticleCloud(const laserParticleCloud&) = delete;


    // Member Functions

        // Access
            inline const fvMesh& mesh() const;

            inline const vector& direction() const;
            inline const vector& centre() const;
            inline scalar radius() const;
            inline scalar power() const;
            inline scalar nRays() const;
            //inline const volScalarField& volLaserPower() const;

        // Edit

            //- Move the particles
            void move();

            //- Generate laser particles from properties
            scalar regenerate();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laserParticleCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "laserParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
