{    
    // derived quantities for numerical integration
    rhoCp == rho * (Cp + CpApp);
    rhoPhiCp == rhoPhi * fvc::interpolate(Cp + CpApp);
    
    T.storePrevIter();
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
            + fvm::div(rhoPhiCp, T)
            - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
            - fvm::laplacian(kappa, T)
            ==
            fvModels.source(rhoCp, T)
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
}
