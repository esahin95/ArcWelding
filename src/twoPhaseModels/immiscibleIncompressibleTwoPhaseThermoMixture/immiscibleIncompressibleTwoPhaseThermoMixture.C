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

\*---------------------------------------------------------------------------*/

#include "immiscibleIncompressibleTwoPhaseThermoMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immiscibleIncompressibleTwoPhaseThermoMixture::
immiscibleIncompressibleTwoPhaseThermoMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleTwoPhaseMixture(U, phi),
    interfaceThermoProperties(alpha1(), alpha2(), U, *this),
    thermo1_(new simpleThermoModel(U.mesh(), this->phase1Name())),
    thermo2_(new simpleThermoModel(U.mesh(), this->phase2Name()))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::immiscibleIncompressibleTwoPhaseThermoMixture::Cp() const 
{
    return alpha1() * thermo1_->Cp() + (1 - alpha1()) * thermo2_->Cp();
}

Foam::tmp<Foam::volScalarField>
Foam::immiscibleIncompressibleTwoPhaseThermoMixture::kappa() const 
{
    return alpha1() * thermo1_->kappa() + (1 - alpha1()) * thermo2_->kappa();
}

Foam::tmp<Foam::volScalarField>
Foam::immiscibleIncompressibleTwoPhaseThermoMixture::beta() const 
{
    return alpha1() * thermo1_->beta() + (1 - alpha1()) * thermo2_->beta();
}


bool Foam::immiscibleIncompressibleTwoPhaseThermoMixture::read()
{
    return
        incompressibleTwoPhaseMixture::read()
     && interfaceThermoProperties::read();
}


// ************************************************************************* //
