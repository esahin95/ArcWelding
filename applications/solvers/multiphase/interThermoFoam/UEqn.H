    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevTau(rho, U)
     ==
       phaseChange.SU(rho, rhoPhi, U)
     + fvModels.source(rho, U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::interpolate(beta*(T - T0)*rho) * snGradGh
                  + fvc::interpolate(dsigmadT * mag(fvc::grad(alpha1))) * fvc::snGrad(T)
                  - fvc::interpolate(dsigmadT * (mixture.n() & fvc::grad(T))) * fvc::snGrad(alpha1)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvConstraints.constrain(U);
    }
