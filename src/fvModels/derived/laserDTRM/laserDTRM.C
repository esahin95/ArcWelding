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

#include "Pstream.H"
#include "laserDTRM.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcGrad.H"
#include "laserParticle.H"

namespace Foam
{

template<class Type>
void gatherAndFlatten(DynamicField<Type>& field)
{
    List<List<Type>> gatheredField(Pstream::nProcs());
    gatheredField[Pstream::myProcNo()] = field;
    Pstream::gatherList(gatheredField);

    field = ListListOps::combine<List<Type>>
    (
        gatheredField, 
        accessOp<List<Type>>()
    );
}

}

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
    lPower_.storePrevIter();
    lPower_ == dimensionedScalar(lPower_.dimensions(), Zero);
    
    // populate cloud
    vector d;
    scalar r;
    for (int i=0; i < nRays_; i++)
    {        
        // sample position
        do 
        {
            d = t1_ * rng_.scalarAB(-radius_, radius_) + t2_ * rng_.scalarAB(-radius_, radius_);
            r = mag(d);
        } 
        while ( r>radius_ || rng_.scalar01() > Foam::exp(-sqr(r / sigma_) / 2.0) );

        // add particle
        d += centre_;
        label cellI = mesh().findCell(d);
        if (cellI != -1)
        {
            cloud_.addParticle(new laserParticle(mesh(), d, q_, direction_, cellI, i, a_));
        }
        
        // warning if particles are outside mesh
        if (returnReduce(cellI, maxOp<label>()) == -1)
        {
            Info << "Cannot find owner cell for position " << d << endl;
        }
    }
    const volScalarField& alpha1 = mesh().lookupObject<const volScalarField>("alpha.steel");

    interpolationCellPoint<scalar> alpha1Interp(alpha1);
    interpolationCellPoint<vector> nHatInterp(fvc::grad(alpha1));
    DynamicField<point> allPositions;
    DynamicField<label> allTracks;
    DynamicField<scalar> allPowers;

    laserParticle::trackingData td
    (
        cloud_, 
        alpha1Interp, 
        nHatInterp,
        lPower_,
        allPositions,
        allTracks,
        allPowers
    );

    // ray tracing
    cloud_.move(cloud_, td, 10.0);
    Info<< "remaining particles: " << cloud_.size() << endl;
    
    // relax power deposition
    lPower_.relax(relax_);
    Info<<"Total Laser Power Deposited in Field: "<< gSum(lPower_) << endl;

    // write rays to vtk
    if (mesh().time().writeTime())
    {
        if (Pstream::parRun())
        {
            gatherAndFlatten(allPositions);
            gatherAndFlatten(allTracks);
            gatherAndFlatten(allPowers);
        }
        
        if (Pstream::master())
        {
            formatterPtr_->write
            (
                outputPath_,
                mesh().time().timeName(),
                coordSet(allTracks, word::null, allPositions),
                "Power",
                allPowers
            );
        }
    }

    // update time level
    curTimeIndex_ = mesh().time().timeIndex();
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
    sigma_(dict.lookupOrDefault("sigma", radius_ * 10.0)),
    centre_(dict.lookupOrDefault("centre", mesh().bounds().midpoint())),
    direction_(dict.lookupOrDefault("direction", vector(0.0, -1.0, 0.0))),
    t1_(normalised(perpendicular(direction_))),
    t2_(normalised(t1_ ^ direction_)),
    q_(Q_ / max(nRays_, 1)),
    relax_(dict.lookupOrDefault("relax", 1.0)),
    a_(dict.lookupOrDefault("absorption", 1.0)),
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
    curTimeIndex_(-1),
    formatterPtr_(setWriter::New("vtk", dict)),
    outputPath_(mesh().time().globalPath()/"postProcessing"/name),
    rng_(Random(Foam::clock::getTime()))
{
    readCoeffs();

    // create directory
    mkDir(outputPath_);
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
    // update heat deposition with ray tracing
    update();

    // add power source
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
