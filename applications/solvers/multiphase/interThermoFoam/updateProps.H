// heat capacity field
Cp == alpha1 * Cp1 + alpha2 * Cp2;

// apparent heat capacity
CpApp == alpha1 * L * fb * Foam::exp(-Foam::sqr(fa * (T - Tmid)));

// heat conductivity
kappa == alpha1 * kappa1 + alpha2 * kappa2;

// thermal expansion coefficient
beta == alpha1 * beta1 + alpha2 * beta2;

// liquid fraction
fliq == scalar(0.5) * (Foam::erf(fa*(T - Tmid)) + scalar(1.0));