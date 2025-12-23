if (adjustTimeStep)
{
    const volScalarField::Internal DiNumfv
    (
        fvc::surfaceSum
        (
            mesh.magSf()
          * fvc::interpolate(mixture.kappa())
          * mesh.surfaceInterpolation::deltaCoeffs()
        )()() / (mesh.V() * rho * mixture.Cp().ref()) * runTime.deltaT()
    );
    scalar DiNum(gMax(DiNumfv));
    Info<< "Diffusion number maximum: " << DiNum << endl;
    
    scalar deltaT = min
    (
        min
        (
            maxCo/(CoNum + small), 
            maxAlphaCo/(alphaCoNum + small)
        ),
        maxDi/(DiNum + small)
    ) * runTime.deltaTValue();
    
    deltaT = min(deltaT, fvModels.maxDeltaT());

    deltaT = min
    (
        min
        (
            deltaT, 
            runTime.deltaTValue() + 0.1*deltaT
        ),
        1.2*runTime.deltaTValue()
    );
    runTime.setDeltaT(min(deltaT, maxDeltaT));

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}