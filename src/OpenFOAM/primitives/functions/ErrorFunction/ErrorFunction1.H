#ifndef ErrorFunction_H
#define ErrorFunction_H

#include "Function1.H"

namespace Foam
{
namespace Function1s
{

template<class Type>
class ErrorFunction
:
    public FieldFunction1<Type, ErrorFunction<Type>>
{
    
    scalar Tliq_;
    scalar Tsol_;
    scalar Tmid_;
    scalar a_;

public:

    TypeName("errorFunction");

    ErrorFunction(const word& name, const dictionary& dict);

    ErrorFunction(const ErrorFunction& errf);

    virtual ~ErrorFunction();

    virtual Type value(const scalar x) const;

    virtual Type integral(const scalar x1, const scalar x2) const;
    
    virtual void write(Ostream& os) const;

    void operator=(const ErrorFunction<Type>&) = delete;

};

}
}

/*
#ifdef NoRepository
    #include "ErrorFunction1.C"
#endif
*/

#endif