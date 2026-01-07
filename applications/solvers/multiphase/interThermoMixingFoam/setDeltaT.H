if (adjustTimeStep)
{
    scalar DiNum = 0.0;
    scalar meanDiNum = 0.0;

    {
        // Energy transport over faces
        scalarField DiNumf
        (
            fvc::surfaceSum
            (
                  mesh.magSf() 
                * fvc::interpolate(mixture.kappa()) 
                * mesh.surfaceInterpolation::deltaCoeffs()
            )().primitiveField()
        );

        // Energy saved in volume
        scalarField DiNumv
        (
            mesh.V().field() * rho.field() * mixture.Cp()().field()
        );

        DiNum = gMax(DiNumf / DiNumv) * runTime.deltaTValue();
        meanDiNum = gSum(DiNumf) / gSum(DiNumv) * runTime.deltaTValue();
    }

    Info<< "Diffusion Number mean: " << meanDiNum
        << " max: " << DiNum << endl;
    
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