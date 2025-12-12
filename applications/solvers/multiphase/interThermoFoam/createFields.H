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

// temperature reference
const dimensionedScalar T0("T0", dimTemperature, pimple.dict().lookupOrDefault<scalar>("T0", 300.0));
Info<< "Set temperature reference to " << T0.value() << endl;

Info<< "Reading additional phase properties\n" << endl;
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

// solidification melting
const autoPtr<Function1<scalar>> alphaSolidT(Function1<scalar>::New("alphaSolid", phaseProperties));

const dimensionedScalar L
(
    "L",
    dimensionSet(0, 2, -2, 0, 0),
    phaseProperties.lookupOrDefault("L", 0.0)
);

const dimensionedScalar DC 
(
    "DC",
    dimensionSet(1,-3,-1,0,0),
    phaseProperties.lookupOrDefault("DC", 0.0)
);

volScalarField alphaSolid
(
    IOobject
    (
        "alphaSolid",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);
alphaSolid.field() = alphaSolidT->value(T.field());
alphaSolid.correctBoundaryConditions();

volScalarField rhoCp
(
    "rhoCp", 
    mixture.rho() * mixture.Cp()
);
rhoCp.oldTime();

surfaceScalarField rhoPhiCp
(
    "rhoPhiCp",
    fvc::interpolate(rhoCp) * phi
);

surfaceScalarField corrf
(
    "corrf",
    scalar(2.0) * fvc::interpolate(mixture.rho()) / (mixture.rho1() + mixture.rho2())
);

// --- End