#ifndef CMVECTORIZED_H
#define CMVECTORIZED_H


#include "cubicmatrix.h"


class CMVectorized: public CubicMatrix
{
public:
    CMVectorized();
    virtual ~CMVectorized();

    virtual const CubicMatrix& dot(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3);
    virtual const CubicMatrix& dotBeta(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3);

private:

};


class CMVectorizedFactory: public CMFactory
{
public:
    virtual ~CMVectorizedFactory() {}
    virtual CubicMatrix* create() const { return new CMVectorized(); }
};


#endif // CMVECTORIZED_H
