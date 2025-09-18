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

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0,1,-1,0,0), Vector<double>(1.0, 0.0, 0.0)),
    "fixedValue"
); 

// density for laser particles
volScalarField alpha1
(
    IOobject
    (
        "alpha1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha1", dimensionSet(1,-3,0,0,0), scalar(0.0)),
    "zeroGradient"
); 
forAll(alpha1, cellI)
{
    const scalar x = mesh.C()[cellI].x();
    const scalar y = mesh.C()[cellI].y();
    alpha1[cellI] = y < x ? 1.0 : 0.0;
}
alpha1.correctBoundaryConditions();

// define laser particles
//Info<< "Adding cloud of laser particles" << endl;
//laserParticleCloud cloud(mesh);
//cloud.regenerate();

surfaceScalarField phi("phi", fvc::interpolate(U) & mesh.Sf());

Info<< "Reading diffusivity k\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
dimensionedScalar k
{
    "k",
    dimensionSet(0,2,-1,0,0),
    transportProperties.lookupOrDefault("k", 1.0)
};
Info<< "k = " << k << endl;

#include "createFvModels.H"