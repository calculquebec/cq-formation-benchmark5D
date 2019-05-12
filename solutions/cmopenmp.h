#ifndef CMOPENMP_H
#define CMOPENMP_H


#include "cubicmatrix.h"


class CMOpenMP: public CubicMatrix
{
public:
    CMOpenMP();
    virtual ~CMOpenMP();

    virtual const CubicMatrix& dot(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3);
    virtual const CubicMatrix& dotBeta(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3);

    virtual void recenter3D(const double value = 0.0);

    virtual CubicMatrix& T();

private:

};


class CMOpenMPFactory: public CMFactory
{
public:
    virtual ~CMOpenMPFactory() {}
    virtual CubicMatrix* create() const { return new CMOpenMP(); }
};


#endif // CMOPENMP_H
