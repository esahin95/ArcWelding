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

forAll(U.boundaryField(), patchI)
{
    forAll(U.boundaryField()[patchI], faceI)
    {
        Info<< patchI << " : " << U.boundaryField()[patchI][faceI] << endl;
    }
}

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

fvScalarMatrix TEqn
(
    -fvm::laplacian(k, T)
);