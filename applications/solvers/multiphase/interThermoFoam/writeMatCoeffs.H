forAll(T, cellI)
{
    T[cellI] = scalar(cellI);
}


{
    label NC = mesh.nCells();
    Info<< "number of Cells " << NC << endl;
    simpleMatrix<scalar> A(NC);

    fvScalarMatrix TEqn
    (
        - fvm::laplacian(kappa, T)
    );

    for (label i=0; i<NC; i++)
    {
        A.source()[i] = 0.0;
        for (label j=0; j<NC; j++)
        {
            A[i][j] = 0.0;
        }
    }

    for (label i=0; i<NC; i++)
    {
        A[i][i] = TEqn.diag()[i];
    }

    for (label faceI=0; faceI<TEqn.lduAddr().lowerAddr().size(); faceI++)
    {
        label u = TEqn.lduAddr().upperAddr()[faceI];
        label l = TEqn.lduAddr().lowerAddr()[faceI];
        A[l][u] = TEqn.upper()[faceI];
        A[u][l] = TEqn.lower()[faceI];
    }

    forAll(T.boundaryField(), patchI)
    {
        const fvPatch& patch = T.boundaryField()[patchI].patch();
        forAll(patch, faceI)
        {
            label cellI = patch.faceCells()[faceI];
            A[cellI][cellI] += TEqn.internalCoeffs()[patchI][faceI];
            A.source()[cellI] += TEqn.boundaryCoeffs()[patchI][faceI];
        }
    }

    Info<<"\nCoefficients for Matrix TEqn\n"<< endl;
    Info << "[ ";
    for (label i=0; i<NC; i++)
    {
        for (label j=0; j<NC; j++)
        {
            Info<< A[i][j] << " ";
        }
        Info<< A.source()[i] << ";" << endl;
    }
    Info << "]" << endl;

    Info<< "min/max(T) = " << min(T).value() << " , " << max(T).value() << endl;
    TEqn.solve();
    Info<< "min/max(T) = " << min(T).value() << " , " << max(T).value() << endl;
}

{
    // resolve TEqn
    fvScalarMatrix TEqn
    (
        - fvm::laplacian(kappa, T)
    );

    TEqn.solve();
    Info<< "min/max(T) = " << min(T).value() << " , " << max(T).value() << endl;
}