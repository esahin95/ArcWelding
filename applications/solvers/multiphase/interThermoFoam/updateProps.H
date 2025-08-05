// heat capacity field
Cp == alpha1 * Cp1 + alpha2 * Cp2;

// heat conductivity
kappa == alpha1 * kappa1 + alpha2 * kappa2;

// thermal expansion coefficient
beta == alpha1 * beta1 + alpha2 * beta2;

// testing
//tmpScl == Cp / dimensionedScalar("tmpCp", dimensionSet(0, 2, -2, -1, 0), 1.0);
tmpVec == fvc::reconstruct(- fvc::interpolate(beta*(T - T0)*rho) * (g & mesh.Sf()) / dimensionedScalar("tmp", dimensionSet(1, -2, -2, 0, 0), 1.0));
tmpVec -= fvc::reconstruct(- fvc::interpolate(beta*(T - T0)*rho) * snGradGh * mesh.magSf()) / dimensionedScalar("tmp", dimensionSet(1, -2, -2, 0, 0), 1.0);
//tmpScl == fvc::snGrad(gh) / dimensionedScalar("tmp", dimensionSet(0, 1, -2, 0, 0), 1.0);
