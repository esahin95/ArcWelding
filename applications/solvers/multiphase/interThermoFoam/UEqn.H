MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
    + MRF.DDt(rho, U)
    + turbulence.divDevTau(rho, U)
    + fvm::Sp(alpha1*DC*Foam::sqr(alphaSolid)/(pow3(1-alphaSolid) + scalar(0.001)), U)
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
                mixture.surfaceTensionForce() * corrf
                - ghf*fvc::snGrad(rho)
                - fvc::interpolate(beta*(T - T0)*rho) * snGradGh
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvConstraints.constrain(U);
}