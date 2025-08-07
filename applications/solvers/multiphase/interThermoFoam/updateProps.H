// heat capacity field
Cp == alpha1 * Cp1 + alpha2 * Cp2;

// heat conductivity
kappa == alpha1 * kappa1 + alpha2 * kappa2;

// thermal expansion coefficient
beta == alpha1 * beta1 + alpha2 * beta2;

// correction field for surface forces
damper == scalar(2.0) * fvc::interpolate(rho) / (rho1 + rho2);