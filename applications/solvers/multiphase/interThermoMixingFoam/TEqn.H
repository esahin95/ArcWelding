{        
    rhoCp == rho * mixture.Cp();
    rhoPhiCp == rhoPhi * fvc::interpolate(mixture.Cp());

    volScalarField L(mixture.L());
    
    T.storePrevIter();
    while (pimple.correctNonOrthogonal())
    {                
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
            + fvm::div(rhoPhiCp, T)
            - fvm::SuSp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
            - fvm::laplacian(mixture.kappa(), T)
            ==
            fvModels.source(rhoCp, T)
            + fvc::ddt(rho, L) + fvc::div(rhoPhi, L, "div(rhoPhiCp,T)")
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