#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading phaseProperties\n" << endl;
immiscibleIncompressibleTwoPhaseThermoMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);
rho.oldTime();

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

volScalarField& alpha1(mixture.alpha1());
mesh.schemes().setFluxRequired(alpha1.name());


#include "createAlphaFluxes.H"

// Construct incompressible turbulence model
incompressibleInterPhaseTransportModel turbulence
(
    U,
    phi,
    alphaPhi1,
    mixture
);


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// --- Additions

Info<< "Reading thermophysical properties\n" << endl;

IOdictionary phaseProperties
(
    IOobject
    (
        "phaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary phase1Dict 
(
    IOobject
    (
        IOobject::groupName
        (
            physicalProperties::typeName,
            mixture.phase1Name()
        ),
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary phase2Dict 
(
    IOobject
    (
        IOobject::groupName
        (
            physicalProperties::typeName,
            mixture.phase2Name()
        ),
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// temperature reference
const dimensionedScalar T0("T0", dimTemperature, pimple.dict().lookupOrDefault<scalar>("T0", 300.0));
Info<< "Set temperature reference to " << T0.value() << endl;

const dimensionedScalar Cp1
(
    "Cp1",
    dimensionSet(0, 2, -2, -1, 0),
    phase1Dict.lookup("Cp")
);
const dimensionedScalar Cp2
(
    "Cp2",
    dimensionSet(0, 2, -2, -1, 0),
    phase2Dict.lookup("Cp")
);
const dimensionedScalar kappa1
(
    "kappa1",
    dimensionSet(1, 1, -3, -1, 0),
    phase1Dict.lookup("kappa")
);
const dimensionedScalar kappa2
(
    "kappa2",
    dimensionSet(1, 1, -3, -1, 0),
    phase2Dict.lookup("kappa")
);
const dimensionedScalar beta1
(
    "beta1",
    dimensionSet(0, 0, 0, -1, 0),
    phase1Dict.lookupOrDefault("beta", 0.0)
);
const dimensionedScalar beta2
(
    "beta2",
    dimensionSet(0, 0, 0, -1, 0),
    phase2Dict.lookupOrDefault("beta", 0.0)
);

const dimensionedScalar Tliq
(
    "Tliq",
    dimensionSet(0, 0, 0, 1, 0),
    phase1Dict.lookupOrDefault("Tliq", scalar(1000.0))
);
const dimensionedScalar Tsol
(
    "Tsol",
    dimensionSet(0, 0, 0, 1, 0),
    phase1Dict.lookupOrDefault("Tsol", Tliq.value() - scalar(50.0))
);
const dimensionedScalar L
(
    "L",
    dimensionSet(0, 2, -2, 0, 0),
    phase1Dict.lookupOrDefault("L", 0.0)
);
const dimensionedScalar DC 
(
    "DC",
    dimensionSet(1,-3,-1,0,0),
    phase1Dict.lookupOrDefault("DC", 1.0e+10)
);
Info<< "Using solid liquid phase change with:\n" 
    << "Tsol: " << Tsol << endl
    << "Tliq: " << Tliq << endl
    << "L: " << L << endl
    << "DC: " << DC << endl;

const dimensionedScalar Tmid
(
    "Tmid",
    scalar(0.5) * (Tliq + Tsol)
);

const dimensionedScalar fa
(
    "fa",
    scalar(4.0) / (Tliq - Tsol)
);

const dimensionedScalar fb
(
    "fb",
    fa / Foam::sqrt(constant::mathematical::pi)
);

volScalarField fliq
(
    "fliq",
    scalar(0.5) * (Foam::erf(fa*(T - Tmid)) + scalar(1.0))
);

volScalarField rhoCp
(
    "rhoCp", 
    mixture.rho() * (mixture.alpha1() * Cp1 + mixture.alpha2() * Cp2)
);
rhoCp.oldTime();

surfaceScalarField rhoPhiCp
(
    "rhoPhiCp",
    fvc::interpolate(rhoCp) * phi
);

volScalarField Cp 
(
    "Cp",
    mixture.alpha1() * Cp1 + mixture.alpha2() * Cp2
);
volScalarField CpApp
(
    "Cp",
    scalar(0.0) * Cp
);

volScalarField kappa 
(
    "kappa",
    mixture.alpha1() * kappa1 + mixture.alpha2() * kappa2
);

volScalarField beta 
(
    "beta",
    mixture.alpha1() * beta1 + mixture.alpha2() * beta2
);

surfaceScalarField corrf
(
    "corrf",
    scalar(2.0) * fvc::interpolate(mixture.rho()) / (mixture.rho1() + mixture.rho2())
);

surfaceScalarField snGradGh
(
    "snGradGh",
    g & mesh.Sf() / mesh.magSf()
);

/*
volVectorField tmpVectorField
(
    IOobject
    (
        "tmpVectorField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(mixture.surfaceTensionForce() * corrf * mesh.magSf() * 0.0)
);

volScalarField tmpScalarField
(
    IOobject
    (
        "tmpScalarField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fliq
);
*/

// --- End