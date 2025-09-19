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

\*---------------------------------------------------------------------------*/

#include "laserDTRM.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "Scale.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "meshSearch.H"
#include "fvcGrad.H"
#include "laserParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(laserDTRM, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        laserDTRM,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::laserDTRM::readCoeffs()
{}

void Foam::fv::laserDTRM::update() const
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }
    lPower_.oldTime();
    lPower_ == dimensionedScalar(lPower_.dimensions(), Zero);
    
    Random rng(12345);

    const meshSearch searchEngine(mesh());

    // construct basis
    const vector t1 = normalised(perpendicular(direction_));
    const vector t2 = normalised(t1 ^ direction_);
    
    for (int i=0; i< nRays_; i++)
    {
        // sample position
        const scalar r(radius_ * rng.scalar01());
        const scalar w(2*Foam::constant::mathematical::pi * rng.scalar01());
        const vector d(r*Foam::sin(w)*t1 + r*Foam::cos(w)*t2 + centre_);
        
        // distributed power 
        const scalar p(Q_ / nRays_);

        // add particle
        label cellI = searchEngine.findCell(d);
        if (cellI != -1)
        {
            cloud_.addParticle(new laserParticle(mesh(), d, p, direction_, cellI));
        }
    }
    Info<< "number of initial particles is " << cloud_.size() << endl;

    const volScalarField& alpha1 = mesh().lookupObject<const volScalarField>("alpha1");

    interpolationCellPoint<scalar> alpha1Interp(alpha1);
    interpolationCellPoint<vector> nHatInterp(fvc::grad(alpha1));

    laserParticle::trackingData td
    (
        cloud_, 
        alpha1Interp, 
        nHatInterp,
        lPower_
    );

    /*
    while (cloud_.size() > 0)
    {
        cloud_.move(cloud_, td, mesh().time().deltaTValue());
    }
    */
    cloud_.move(cloud_, td, 1.0);
    Info<<"Total Laser Power Deposited in Field "<< Foam::sum(lPower_)<<endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::laserDTRM::laserDTRM
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    Q_(dict.lookupOrDefault("Q", 0.0)),
    nRays_(dict.lookupOrDefault("nRays", 0)),
    radius_(dict.lookupOrDefault("radius", 1.0)),
    centre_(dict.lookupOrDefault("centre", mesh().bounds().midpoint())),
    direction_(dict.lookupOrDefault("direction", vector(0.0, -1.0, 0.0))),
    cloud_(mesh, "cloudDTRM", IDLList<laserParticle>()),
    lPower_
    (
        IOobject
        (
            IOobject::groupName(name, "power"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1,-1,-3,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    curTimeIndex_(-1)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::laserDTRM::~laserDTRM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::laserDTRM::addSupFields() const
{
    return wordList({"T"});
}


void Foam::fv::laserDTRM::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    update();

    eqn.source() -= lPower_.primitiveField();
}


void Foam::fv::laserDTRM::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    addSup(eqn, fieldName);
}


bool Foam::fv::laserDTRM::movePoints()
{
    return true;
}


void Foam::fv::laserDTRM::topoChange(const polyTopoChangeMap& map)
{}


void Foam::fv::laserDTRM::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::laserDTRM::distribute(const polyDistributionMap& map)
{}


bool Foam::fv::laserDTRM::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
