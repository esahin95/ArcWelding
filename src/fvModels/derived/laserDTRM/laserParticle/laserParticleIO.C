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

#include "laserParticle.H"
#include "IOstreams.H"

const std::size_t Foam::laserParticle::sizeofFields_
(
    sizeof(laserParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laserParticle::laserParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            d_ = readScalar(is);
            is >> d0_ 
               >> U_ 
               >> trackIndex_ 
               >> a_ 
               >> maxReflections_ 
               >> reflections_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&d_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check("laserParticle::laserParticle(Istream&)");
}

Foam::Ostream& Foam::operator<<(Ostream& os, const laserParticle& p)
{
    
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.d0_
            << token::SPACE << p.U_
            << token::SPACE << p.trackIndex_
            << token::SPACE << p.a_
            << token::SPACE << p.maxReflections_
            << token::SPACE << p.reflections_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
            laserParticle::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const laserParticle&)");

    return os;
}

// ************************************************************************* //