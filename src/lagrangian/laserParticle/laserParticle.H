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
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

        //- Power
        scalar d_;

        //- Initial power
        scalar d0_;

        //- Direction
        vector U_;

        //- Track index
        label trackIndex_;

public:

    friend class Cloud<laserParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& alpha1Interp_;
            const interpolationCellPoint<vector>& nHatInterp_;

        // reference to heat source field

            volScalarField& lPower_;

        // streamline data

            DynamicField<point>& allPositions_;
            DynamicField<label>& allTracks_;
            DynamicField<scalar>& allPowers_;

    public:

        // Constructors

            inline trackingData
            (
                const Cloud<laserParticle>& spc,
                const interpolationCellPoint<scalar>& alpha1Interp,
                const interpolationCellPoint<vector>& nHatInterp,
                volScalarField& lPower,
                DynamicField<point>& allPositions,
                DynamicField<label>& allTracks,
                DynamicField<scalar>& allPowers
            );


        // Member Functions

            inline const interpolationCellPoint<scalar>& alpha1Interp() const;
            inline const interpolationCellPoint<vector>& nHatInterp() const;

            inline scalar& lPower(label cellI);

            inline void append(vector position, label trackIndex, scalar d);
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
            const vector& U,
            const label trackIndex
        );

        //- Construct from position
        inline laserParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const scalar d,
            const vector& U,
            const label cellI,
            const label trackIndex
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

            //- Return track index
            inline label index() const;


        // Tracking

            //- Move
            bool move(Cloud<laserParticle>&, trackingData&, const scalar);


        // Patch interactions

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(Cloud<laserParticle>& cloud, trackingData& td);

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
