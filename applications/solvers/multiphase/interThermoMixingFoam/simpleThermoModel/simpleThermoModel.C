#include "simpleThermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleThermoModel::simpleThermoModel(const fvMesh& mesh, const word& group)
:
    physicalProperties(mesh, group),
    mesh_(mesh),
    Cp0_("Cp", dimensionSet(0, 2, -2, -1, 0), *this),
    kappa0_("kappa", dimensionSet(1, 1, -3, -1, 0), *this),
    beta0_("beta", dimensionSet(0, 0, 0, -1, 0), *this),
    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", group),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        Cp0_
    ),
    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", group),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        kappa0_
    ),
    beta_
    (
        IOobject
        (
            IOobject::groupName("beta", group),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        beta0_
    ),
    L_("L", dimensionSet(0, 2, -2, 0, 0), this->lookupOrDefault("L", 0.0)),
    alphaSolidPtr_(Function1<scalar>::New("alphaSolid", *this)),
    alphaSolid_ 
    (
        IOobject 
        (
            IOobject::groupName("alphaSolid", group),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    )
{
    correct();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::simpleThermoModel::Cp() const
{
    return Cp_;
}

Foam::tmp<Foam::volScalarField>
Foam::simpleThermoModel::kappa() const
{
    return kappa_;
}

Foam::tmp<Foam::volScalarField>
Foam::simpleThermoModel::beta() const
{
    return beta_;
}

Foam::tmp<Foam::volScalarField>
Foam::simpleThermoModel::alphaSolid() const
{
    return alphaSolid_;
}

bool Foam::simpleThermoModel::read()
{
    return physicalProperties::read();
}

void Foam::simpleThermoModel::correct()
{
    // temperature reference
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const volScalarField::Boundary& TBf = T.boundaryField();

    // correct internal fields
    alphaSolid_.field() = alphaSolidPtr_->value(T.field());

    // correct boundary fields
    volScalarField::Boundary& alphaSolidBf = alphaSolid_.boundaryFieldRef();
    forAll(alphaSolidBf, patchI)
    {
        alphaSolidBf[patchI] = alphaSolidPtr_->value(TBf[patchI]);
    }
}