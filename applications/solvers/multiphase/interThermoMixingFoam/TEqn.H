{        
    rhoCp == rho * mixture.Cp();
    rhoPhiCp == rhoPhi * fvc::interpolate(mixture.Cp());

    //volScalarField qL(mixture.qL());
    
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
            + mixture.qL()//+ fvc::ddt(rho, L) //+ fvc::div(rhoPhi, L, "div(rhoPhiCp,T)")
        );

        TEqn.relax();
        
        fvConstraints.constrain(TEqn);
        
        TEqn.solve();

        fvConstraints.constrain(T);

        if (pimple.finalNonOrthogonalIter())
        {
            Info<< "max change in T before relaxation: " << gMax(mag(T.prevIter().field()-T.field())) << endl;
    
            T.relax();
        }
    }
}