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

#include "spheresToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(spheresToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, spheresToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::spheresToCell::combine(topoSet& set, const bool add) const
{
    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, cellI)
    {
        bool centreInSphere = false;

        forAll(centres_, sphereI)
        {
            scalar offset = magSqr(centres_[sphereI] - ctrs[cellI]);
            if (offset <= radii2_[sphereI])
            {
                centreInSphere = true;
                break;
            }
        }

        addOrDelete(set, cellI, centreInSphere);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spheresToCell::spheresToCell
(
    const polyMesh& mesh,
    const DynamicField<point>& centres,
    const DynamicField<scalar>& radii2,
    const scalar offset
)
:
    topoSetSource(mesh),
    centres_(centres),
    radii2_(radii2),
    offset_(offset)
{}


Foam::spheresToCell::spheresToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    centres_(),
    radii2_(),
    offset_(dict.lookupOrDefault("offset", 0.0))
{
    // Raw list of data
    List<scalar> data(dict.lookup("data"));
    
    // Read dynamic fields
    for(int i = 0; i < data.size(); i += 5)
    {
        centres_.append(point(data[i], data[i+1], data[i+2]));
        radii2_.append(sqr(scalar(data[i+3]) + offset_));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spheresToCell::~spheresToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::spheresToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with centre within any of " << centres_.size() << " spheres" << endl;

        combine(set, true);
    }
    else
    {
        Info<< "    Given action = " << action << " not supported" << endl;
    }
}


// ************************************************************************* //
