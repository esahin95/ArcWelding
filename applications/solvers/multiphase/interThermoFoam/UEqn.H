MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
    + MRF.DDt(rho, U)
    + turbulence.divDevTau(rho, U)
    + fvm::Sp(alpha1*DC*Foam::sqr(alphaSolid)/(pow3(1-alphaSolid) + scalar(0.001)), U)
    ==
    fvModels.source(rho, U)
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
                mixture.surfaceTensionForce() * corrf
                - ghf*fvc::snGrad(rho)
                - fvc::interpolate(mixture.beta()*(T - T0)*rho) * (g & mesh.Sf() / mesh.magSf())
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvConstraints.constrain(U);
}