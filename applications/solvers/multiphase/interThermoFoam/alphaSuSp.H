// Phase change alpha1 source
const tmp<volScalarField::Internal>& Su = tmp<volScalarField::Internal>(nullptr);
const tmp<volScalarField::Internal>& Sp = tmp<volScalarField::Internal>(nullptr);

tmp<volScalarField> divU;

if (Sp.valid())
{
    divU =
    (
        mesh.moving()
      ? fvc::div(phiCN() + mesh.phi())
      : fvc::div(phiCN())
    );
}
