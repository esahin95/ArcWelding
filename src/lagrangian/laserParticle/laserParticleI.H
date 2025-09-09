/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::laserParticle::trackingData::trackingData
(
    const laserParticleCloud& spc,
    const interpolationCellPoint<scalar>& rhoInterp
)
:
    particle::trackingData(spc),
    rhoInterp_(rhoInterp)
{}


inline Foam::laserParticle::laserParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const scalar d,
    const vector& U
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    d_(d),
    U_(U)
{}

inline Foam::laserParticle::laserParticle
(
    const polyMesh& mesh,
    const vector& position,
    const scalar d,
    const vector& U
)
:
    particle(mesh, position, label(-1)),
    d_(d),
    U_(U)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::laserParticle::trackingData::rhoInterp() const
{
    return rhoInterp_;
}

inline Foam::scalar Foam::laserParticle::d() const
{
    return d_;
}

inline const Foam::vector& Foam::laserParticle::U() const
{
    return U_;
}


// ************************************************************************* //
