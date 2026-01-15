/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "incompressibleThreePhaseThermoMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleThreePhaseThermoMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(alpha1_*rho1_ + alpha2_*rho2_ + alpha3_*rho3_);
}

void Foam::incompressibleThreePhaseThermoMixture::calcThermo()
{
    thermo1_->correct();
    thermo2_->correct();
    thermo3_->correct();

    // Update private field data
    Cp_ == 
    (
          alpha1_ * thermo1_->Cp() 
        + alpha2_ * thermo2_->Cp() 
        + alpha3_ * thermo3_->Cp()
    );

    kappa_ == 
    (
          alpha1_ * thermo1_->kappa() 
        + alpha2_ * thermo2_->kappa() 
        + alpha3_ * thermo3_->kappa()
    );

    beta_ == 
    (
          alpha1_ * thermo1_->beta() 
        + alpha2_ * thermo2_->beta() 
        + alpha3_ * thermo3_->beta()
    );

    alphaSolid_ == 
    (
          alpha1_ * thermo1_->alphaSolid() 
        + alpha2_ * thermo2_->alphaSolid() 
        + alpha3_ * thermo3_->alphaSolid()
    );
    alphaSolid_ == min(scalar(1.0), max(scalar(0.0), alphaSolid_));

    
    L_ == 
    (
          alpha1_ * thermo1_->L() * thermo1_->alphaSolid() 
        + alpha2_ * thermo2_->L() * thermo2_->alphaSolid() 
        + alpha3_ * thermo3_->L() * thermo3_->alphaSolid()
    );
    L_ == max(dimensionedScalar(L_.dimensions(), 0.0), L_);
    
    /*
    rhoL_ ==
    (
          alpha1_ * rho1_ * thermo1_->L() * thermo1_->alphaSolid() 
        + alpha2_ * rho2_ * thermo2_->L() * thermo2_->alphaSolid() 
        + alpha3_ * rho3_ * thermo3_->L() * thermo3_->alphaSolid()
    );
    */
    
    qL_ ==
    (
          alpha1_ * thermo1_->L() * rho1_ * fvc::ddt(thermo1_->alphaSolid()())
        + alpha2_ * thermo2_->L() * rho2_ * fvc::ddt(thermo2_->alphaSolid()())
        + alpha3_ * thermo3_->L() * rho3_ * fvc::ddt(thermo3_->alphaSolid()())
    );
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleThreePhaseThermoMixture::incompressibleThreePhaseThermoMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phase1Name_(wordList(lookup("phases"))[0]),
    phase2Name_(wordList(lookup("phases"))[1]),
    phase3Name_(wordList(lookup("phases"))[2]),

    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha3_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U.time().timeName(),
            U.db()
        ),
        U.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),

    nuModel1_(viscosityModel::New(U.mesh(), phase1Name_)),
    nuModel2_(viscosityModel::New(U.mesh(), phase2Name_)),
    nuModel3_(viscosityModel::New(U.mesh(), phase3Name_)),

    rho1_("rho", dimDensity, nuModel1_()),
    rho2_("rho", dimDensity, nuModel2_()),
    rho3_("rho", dimDensity, nuModel3_()),

    thermo1_(new simpleThermoModel(U.mesh(), phase1Name_)),
    thermo2_(new simpleThermoModel(U.mesh(), phase2Name_)),
    thermo3_(new simpleThermoModel(U.mesh(), phase3Name_)),

    Cp_
    (
        IOobject
        (
            "Cp",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(thermo1_->Cp()().dimensions(), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    kappa_
    (
        IOobject
        (
            "kappa",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(thermo1_->kappa()().dimensions(), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    beta_
    (
        IOobject
        (
            "beta",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(thermo1_->beta()().dimensions(), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    alphaSolid_
    (
        IOobject
        (
            "alphaSolid",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(thermo1_->alphaSolid()().dimensions(), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    L_
    (
        IOobject
        (
            "L",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(thermo1_->L().dimensions(), 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    qL_
    (
        IOobject
        (
            "qL",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(thermo1_->L().dimensions() * rho1_.dimensions() / dimTime, 0.0),
        calculatedFvPatchScalarField::typeName
    )

{
    alpha3_ == 1.0 - alpha1_ - alpha2_;
    calcNu();
    calcThermo();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseThermoMixture::mu() const
{
    return volScalarField::New
    (
        "mu",
        alpha1_*rho1_*nuModel1_->nu()
      + alpha2_*rho2_*nuModel2_->nu()
      + alpha3_*rho3_*nuModel3_->nu()
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseThermoMixture::muf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return surfaceScalarField::New
    (
        "mu",
        alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
      + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
      + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseThermoMixture::nuf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return surfaceScalarField::New
    (
        "nu",
        (
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
          + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
        )/(alpha1f*rho1_ + alpha2f*rho2_ + alpha3f*rho3_)
    );
}

bool Foam::incompressibleThreePhaseThermoMixture::read()
{
    if (regIOobject::read())
    {
        nuModel1_->lookup("rho") >> rho1_;
        nuModel2_->lookup("rho") >> rho2_;
        nuModel3_->lookup("rho") >> rho3_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
