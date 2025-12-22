{        
    rhoCp == rho * mixture.Cp();
    rhoPhiCp == rhoPhi * fvc::interpolate(mixture.Cp());
    
    T.storePrevIter();
    while (pimple.correctNonOrthogonal())
    {        
        Info<<"Check dimensions: "<< endl;
        Info<<"kappa: " << mixture.kappa().ref().dimensions() << endl;
        
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
            + fvm::div(rhoPhiCp, T)
            - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
            - fvm::laplacian(mixture.kappa(), T)
            ==
            fvModels.source(rhoCp, T)
        );

        TEqn.relax();
        
        fvConstraints.constrain(TEqn);
        
        TEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            Info<< "max change in T before relaxation: " << gMax(mag(T.prevIter().field()-T.field())) << endl;
    
            T.relax();
        }
    }
}