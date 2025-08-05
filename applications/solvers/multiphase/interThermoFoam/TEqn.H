{    
    // derived quantities for numerical integration
    rhoCp == rho * Cp;
    rhoPhiCp == rhoPhi * fvc::interpolate(Cp);
    
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
    
    Info<< "min/max(T) = "
        << min(T).value() << " , " << max(T).value() << endl;
}
