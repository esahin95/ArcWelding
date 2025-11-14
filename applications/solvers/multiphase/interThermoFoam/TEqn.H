{    
    // apparent heat capacity
    CpApp == alpha1 * L * fb * Foam::exp(-Foam::sqr(fa * (T - Tmid)));

    // derived quantities for numerical integration
    rhoCp == rho * (Cp + CpApp);
    rhoPhiCp == rhoPhi * fvc::interpolate(Cp + CpApp);
    
    // build linear equation
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
    
    T.relax();
    
    Info<< "min/max(T) = " << min(T).value() << " , " << max(T).value() << endl;
}
