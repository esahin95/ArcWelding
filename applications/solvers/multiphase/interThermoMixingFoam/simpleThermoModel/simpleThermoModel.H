#ifndef simpleThermoModel_H
#define simpleThermoModel_H

#include "physicalProperties.H"

#include "volFields.H"

#include "Function1.H"

namespace Foam
{

//Forward declaration of classes
//class FvMesh;

class simpleThermoModel
:
    public physicalProperties
{

    const fvMesh& mesh_;
    
    dimensionedScalar Cp0_;
    dimensionedScalar kappa0_;
    dimensionedScalar beta0_;

    volScalarField Cp_;
    volScalarField kappa_;
    volScalarField beta_;

    dimensionedScalar L_;

    autoPtr<Function1<scalar>> alphaSolidPtr_;

public:

    simpleThermoModel(const fvMesh& mesh, const word& group);

    simpleThermoModel(const simpleThermoModel&) = delete;

    static autoPtr<simpleThermoModel> New
    (
        const fvMesh& mesh,
        const word& group = word::null
    );

    virtual ~simpleThermoModel()
    {}

    // member functions

    virtual bool read();

    virtual void correct()
    {};

    tmp<volScalarField> Cp() const;

    tmp<volScalarField> kappa() const;

    tmp<volScalarField> beta() const;

    tmp<volScalarField> alphaSolid() const;

    const dimensionedScalar& L() const
    {
        return L_;
    };
};

}

#endif