/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::fv::laserDTRM

Description
    Model for applying a heat source. Requires the power Q to be specified

Usage
    Example usage:
    \verbatim
    laserDTRM
    {
        type            laserDTRM;

        Q               1e6;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef laserDTRM_H
#define laserDTRM_H

#include "fvModel.H"
#include "Function1.H"
#include "Cloud.H"
#include "setWriter.H"
#include "laserParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
//class meshSearch;
class Random;

namespace fv
{

/*---------------------------------------------------------------------------*\
                Class laserDTRM Declaration
\*---------------------------------------------------------------------------*/

class laserDTRM
:
    public fvModel
{
    // Private data

        //- The heat source
        scalar Q_;
        scalar nRays_;
        scalar radius_;
        scalar sigma_;
        vector centre_;
        vector direction_;
        vector t1_;
        vector t2_;


        //- The lagrangian cloud
        mutable Cloud<laserParticle> cloud_;

        // volume heat source
        mutable volScalarField lPower_;

        //- Current time index
        mutable label curTimeIndex_;

        //- File writer
        autoPtr<setWriter> formatterPtr_;

        //- Random number generator
        mutable Random rng_;

    // Private member functions

        //- Non-virtual read
        void readCoeffs();

        //- Update power deposition
        void update() const;


public:

    //- Runtime type information
    TypeName("laserDTRM");


    // Constructors

        //- Construct from dictionary
        laserDTRM
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~laserDTRM();


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Source term to energy equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Source term to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
