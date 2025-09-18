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

Class
    Foam::laserParticle

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    laserParticleI.H
    laserParticle.C
    laserParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef laserParticle_H
#define laserParticle_H

#include "particle.H"
#include "Cloud.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class laserParticleCloud;


// Forward declaration of friend functions and operators

class laserParticle;

Ostream& operator<<(Ostream&, const laserParticle&);


/*---------------------------------------------------------------------------*\
                        Class laserParticle Declaration
\*---------------------------------------------------------------------------*/

class laserParticle
:
    public particle
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- power
        scalar d_;

        //- direction
        vector U_;


public:

    friend class Cloud<laserParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& alpha1Interp_;
    public:

        // Constructors

            inline trackingData
            (
                const laserParticleCloud& spc,
                const interpolationCellPoint<scalar>& alpha1Interp
            );


        // Member Functions

            inline const interpolationCellPoint<scalar>& alpha1Interp() const;
    };


    //- Runtime type information
    TypeName("laserParticle");


    // Constructors

        //- Construct from components
        inline laserParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const scalar d,
            const vector& U
        );

        //- Construct from position
        inline laserParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const scalar d,
            const vector& U
        );

        //- Construct from Istream
        laserParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new laserParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<laserParticle> operator()(Istream& is) const
            {
                return autoPtr<laserParticle>
                (
                    new laserParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return power
            inline scalar d() const;

            //- Return direction
            inline const vector& U() const;


        // Tracking

            //- Move
            bool move(laserParticleCloud&, trackingData&, const scalar);


        // Patch interactions

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(laserParticleCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation
            virtual void transformProperties(const transformer&);


    // I-O

        static void readFields(Cloud<laserParticle>& c);

        static void writeFields(const Cloud<laserParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const laserParticle&);
};


template<>
inline bool contiguous<laserParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "laserParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
