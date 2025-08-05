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

#include "createPhi.H"


Info<< "Reading phaseProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

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

// Testing ------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

volScalarField tmpScl
(
    IOobject
    (
        "tmpScl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("tmpScl", dimless, 0.0)
);

volVectorField tmpVec
(
    IOobject
    (
        "tmpVec",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("tmpVec", dimless, vector::zero),
    extrapolatedCalculatedFvPatchVectorField::typeName
);

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

// Additions ----------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

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

const dimensionedScalar dsigmadT
(
    "dsigmadT",
    dimensionSet(1, 0, -2, -1, 0),
    phaseProperties.lookupOrDefault("dsigmadT", 0.0)
);

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

Info<< "Creating thermophysical fields\n" << endl;
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

// temperature reference
const dimensionedScalar T0("T0", dimTemperature, pimple.dict().lookupOrDefault<scalar>("T0", 300.0));
Info<< "Set temperature reference to " << T0.value() << endl;

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

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// surface normal gradient of g & h for boussinesq approximation
const surfaceScalarField snGradGh
(
    "snGradGh",
    g & mesh.Sf() / mesh.magSf()
);


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
