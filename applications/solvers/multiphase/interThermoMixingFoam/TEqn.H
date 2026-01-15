{        
    rhoCp == rho * mixture.Cp();
    rhoPhiCp == rhoPhi * fvc::interpolate(mixture.Cp());
    
    volScalarField L(mixture.L());
    while (pimple.correctNonOrthogonal())
    {                
        ddtL == fvc::ddt(rho, L);
        
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
            + fvm::div(rhoPhiCp, T)
            - fvm::SuSp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
            - fvm::laplacian(mixture.kappa(), T)
            ==
            fvModels.source(rhoCp, T)
            + fvc::ddt(rho, L) + fvc::div(rhoPhi, L, "div(rhoPhiCp,T)") // mixture.qL()
        );

        TEqn.relax();
        
        fvConstraints.constrain(TEqn);
        
        TEqn.solve();

        fvConstraints.constrain(T);

        if (pimple.finalNonOrthogonalIter())
        {    
            T.relax();
        }
    }
}