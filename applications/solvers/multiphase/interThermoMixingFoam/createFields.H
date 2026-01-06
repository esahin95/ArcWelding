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

immiscibleIncompressibleThreePhaseThermoMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& rho3 = mixture.rho3();

dimensionedScalar D23("D23", dimViscosity, mixture);

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
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();

// Need to store rhoCp for ddt(rhoCp, T)
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    rho * mixture.Cp() 
);
rhoCp.oldTime();


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);


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
mesh.schemes().setFluxRequired(alpha2.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// temperature reference
const dimensionedScalar T0("T0", dimTemperature, pimple.dict().lookupOrDefault<scalar>("T0", 300.0));
Info<< "Set temperature reference to " << T0.value() << endl;

// numerical parameter for solidification
const dimensionedScalar DC("DC", dimless, pimple.dict().lookupOrDefault<scalar>("DC", 0.0));
Info<< "Set permeability to " << DC.value() << endl;

Info<<"Create solid fraction field"<<endl;
volScalarField alphaSolid
(
    IOobject 
    (
        "alphaSolid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.alphaSolid()
);

surfaceScalarField rhoPhiCp
(
    "rhoPhiCp",
    fvc::interpolate(rhoCp) * phi
);