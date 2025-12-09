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

autoPtr<Function1<scalar>> fliqPtr(Function1<scalar>::New("fliq", transportProperties));
Info<<fliqPtr()<<endl;

/*
forAll(U, cellI)
{
    Info<< cellI << " : " << U[cellI] << endl;
}

forAll(U.boundaryField(), patchI)
{
    forAll(U.boundaryField()[patchI], faceI)
    {
        Info<< patchI << " : " << U.boundaryField()[patchI][faceI] << endl;
    }
}
*/

/*
forAll(T, cellI)
{
    scalar x = mesh.C()[cellI][0];
    T[cellI] = fliqPtr->value(x);
    Info<< "( " << x << " , " << T[cellI] << " , " << fliqPtr->value(x)  << " )" << " ";
}
Info<<endl;
*/

//Info << fliqPtr->value(mesh.C() & vector(1.0,0.0,0.0)) << endl;
T.primitiveFieldRef() = fliqPtr->value(mesh.C() & vector(1.0,0.0,0.0));
T.correctBoundaryConditions();
Info<<T<<endl;

/*
forAll(T.boundaryField(), patchI)
{
    forAll(T.boundaryField()[patchI], faceI)
    {
        Info<< T.boundaryField()[patchI].patch().name() << " : " << T.boundaryField()[patchI][faceI] << endl;
    }
}
*/

surfaceScalarField phi("phi", fvc::interpolate(U) & mesh.Sf());

fvScalarMatrix TEqn
(
    -fvm::laplacian(k, T)
);