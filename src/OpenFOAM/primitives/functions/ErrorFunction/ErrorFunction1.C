#include "ErrorFunction1.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    makeFunction1Type(ErrorFunction, scalar);
}

template<class Type>
Foam::Function1s::ErrorFunction<Type>::ErrorFunction
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction1<Type, ErrorFunction<Type>>(name),
    Tliq_(dict.lookup<scalar>("Tliq")),
    Tsol_(dict.lookup<scalar>("Tsol")),
    Tmid_(0.5 * (Tsol_ + Tliq_)),
    a_(4.0 / (max(1e-6, Tliq_ - Tsol_)))
{}

template<class Type>
Foam::Function1s::ErrorFunction<Type>::ErrorFunction(const ErrorFunction& errf)
:
    FieldFunction1<Type, ErrorFunction<Type>>(errf),
    Tliq_(errf.Tliq_),
    Tsol_(errf.Tsol_),
    Tmid_(errf.Tmid_),
    a_(errf.a_)
{}

template<class Type>
Foam::Function1s::ErrorFunction<Type>::~ErrorFunction()
{}

template<class Type>
Type Foam::Function1s::ErrorFunction<Type>::value(const scalar x) const 
{
    return 0.5 * (Foam::erf(a_ * (x - Tmid_)) + 1.0);
}

/*
template<class Type>
Type Foam::Function1s::ErrorFunction<Type>::derivative(const scalar x) const 
{
    return b_ * Foam::exp(-Foam::sqr(a_ * (x - Tmid_)));
}
*/

template<class Type>
Type Foam::Function1s::ErrorFunction<Type>::integral
(
    const scalar x1, 
    const scalar x2
) const 
{
    Type y(Zero);
    return y;
}

template<class Type>
void Foam::Function1s::ErrorFunction<Type>::write(Ostream& os) const 
{
    writeKeyword(os, "Tsol-Tliq") << Tsol_ << "-" << Tliq_ << token::END_STATEMENT << nl;
}