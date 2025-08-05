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

volScalarField Ttmp
(
    IOobject
    (
        "Ttmp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ttmp", dimTemperature, 50.0),
    "fixedValue"
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

forAll(T.internalField(), celli)
{
    T[celli] = 1.0;
}
T.correctBoundaryConditions();

Ttmp == T;
forAll(Ttmp, celli)
{
    Info<< celli << " : " << U[celli] << endl;
}
Ttmp.correctBoundaryConditions();

/*
const volVectorField::Boundary& TBoundaryList = U.boundaryField();
forAll(TBoundaryList, patchi)
{
    const fvPatchVectorField& fieldBoundary = TBoundaryList[patchi];
    Info<< fieldBoundary.patch().name() << endl;
    forAll(fieldBoundary, facei)
    {
        vector faceT = fieldBoundary[facei];
        Info<< patchi << " : " << faceT << endl;
    }
}
*/
forAll(U.boundaryField(), patchI)
{
    forAll(U.boundaryField()[patchI], faceI)
    {
        Info<< patchI << " : " << U.boundaryField()[patchI][faceI] << endl;
    }
}

surfaceScalarField phi("phi", fvc::interpolate(U) & mesh.Sf());

fvScalarMatrix TEqn
(
    //fvm::laplacian(dimensionedScalar("D", dimensionSet(0,0,0,-1,0,0,0), 1.0), T)
    fvm::div(phi, T)
);

Info<< "diagonal original" << endl;
Info<< TEqn.diag() << endl;
Info<< "lower" << endl;
Info<< TEqn.lower() << endl;
Info<< "upper" << endl;
Info<< TEqn.upper() << endl;
Info<< "lowerAddr" << endl;
Info<< TEqn.lduAddr().lowerAddr() << endl;
Info<< "upperAddr" << endl;
Info<< TEqn.lduAddr().upperAddr() << endl;
Info<< "diagonal negSumDiag" << endl;
TEqn.negSumDiag();
Info<< TEqn.diag() << endl;

Info<< "\nLaplacian\n" << endl;

fvScalarMatrix LEqn
(
    fvm::laplacian(dimensionedScalar("D", dimensionSet(0,0,0,-1,0,0,0), 1.0), T)
);

Info<< "diagonal original" << endl;
Info<< LEqn.diag() << endl;
Info<< "lower" << endl;
Info<< LEqn.lower() << endl;
Info<< "upper" << endl;
Info<< LEqn.upper() << endl;
Info<< "lowerAddr" << endl;
Info<< LEqn.lduAddr().lowerAddr() << endl;
Info<< "upperAddr" << endl;
Info<< LEqn.lduAddr().upperAddr() << endl;
Info<< "diagonal negSumDiag" << endl;
TEqn.negSumDiag();
Info<< LEqn.diag() << endl;

Info<< LEqn.internalCoeffs() << endl;
Info<< LEqn.boundaryCoeffs() << endl;