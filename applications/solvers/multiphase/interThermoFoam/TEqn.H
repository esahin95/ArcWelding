{        
    rhoCp == rho * mixture.Cp();
    rhoPhiCp == rhoPhi * fvc::interpolate(mixture.Cp());
    
    T.storePrevIter();
    while (pimple.correctNonOrthogonal())
    {        
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
            + fvm::div(rhoPhiCp, T)
            - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
            - fvm::laplacian(mixture.kappa(), T)
            ==
            L * fvc::ddt(rho, alphaSolid)
            + fvModels.source(rhoCp, T)
        );

        TEqn.relax();
        
        fvConstraints.constrain(TEqn);
        
        TEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            Info<< "max change in T: " << gMax(mag(T.prevIter().field()-T.field())) << endl;
    
            T.relax();
        }
    }

    // solid fraction
    alphaSolid.field() = min(alphaSolidT->value(T.field()), alpha1);
    alphaSolid.correctBoundaryConditions();
}
