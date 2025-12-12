/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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
    Foam::immiscibleIncompressibleTwoPhaseThermoMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    immiscibleIncompressibleTwoPhaseThermoMixture.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleIncompressibleTwoPhaseThermoMixture_H
#define immiscibleIncompressibleTwoPhaseThermoMixture_H

#include "incompressibleTwoPhaseMixture.H"
#include "interfaceThermoProperties.H"
#include "simpleThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class immiscibleIncompressibleTwoPhaseThermoMixture Declaration
\*---------------------------------------------------------------------------*/

class immiscibleIncompressibleTwoPhaseThermoMixture
:
    public incompressibleTwoPhaseMixture,
    public interfaceThermoProperties
{

    autoPtr<simpleThermoModel> thermo1_;
    autoPtr<simpleThermoModel> thermo2_;

public:

    // Constructors

        //- Construct from components
        immiscibleIncompressibleTwoPhaseThermoMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~immiscibleIncompressibleTwoPhaseThermoMixture()
    {}


    // Member Functions

        //- Correct the phase and interface properties
        virtual void correct()
        {
            incompressibleTwoPhaseMixture::correct();
            interfaceThermoProperties::correct();
            thermo1_->correct();
            thermo2_->correct();
        }

        //- Read base phaseProperties dictionary
        virtual bool read();

        tmp<volScalarField> Cp() const;

        tmp<volScalarField> kappa() const;

        tmp<volScalarField> beta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
