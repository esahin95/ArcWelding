#ifndef simpleThermoModel_H
#define simpleThermoModel_H

#include "physicalProperties.H"

namespace Foam
{

//Forward declaration of classes
class FvMesh;

class simpleThermoModel
:
    public physicalProperties
{


public:

    TypeName("simpleThermoModel");

    simpleThermoModel(const fvMesh& mesh, const word& group);

    simpleThermoModel(const simpleThermoModel&) = delete;

    virtual ~simpleThermoModel()
    {}

    // member functions

    virtual bool read();

    tmp<volScalarField> Cp() const;

    tmp<volScalarField> kappa() const;


}

}

#endif