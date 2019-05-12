#include "cmopenmp.h"

#include <algorithm>
#include <stdlib.h>


using namespace std;


CMOpenMP::CMOpenMP(): CubicMatrix()
{
}


CMOpenMP::~CMOpenMP()
{
}


const CubicMatrix& CMOpenMP::dot(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3)
{
    if (cm1.n() < 1 || cm2.n() < 1 || cm3.n() < 1 ||
            cm1.n() != cm2.n() || cm2.n() != cm3.n()) {
        cleanup();
        return *this;
    }

    resize(cm1.n());

    for (size_t i1 = 0; i1 < n_; i1++) {
        for (size_t j2 = 0; j2 < n_; j2++) {
            for (size_t k3 = 0; k3 < n_; k3++) {
                double sum = 0.0;

                for (size_t s = 0; s < n_; s++) {
                    size_t j1 = s, k2 = s, i3 = s;

                    for (size_t t = 0; t < n_; t++) {
                        size_t k1 = t, i2 = t, j3 = t;

                        sum += cm1.matrix()[i1][j1][k1] *
                               cm2.matrix()[i2][j2][k2] *
                               cm3.matrix()[i3][j3][k3];
                    }
                }

                matrix_[i1][j2][k3] = sum;
            }
        }
    }

    nflop_ += n_* n_ * (n_* n_ + n_* n_* n_* 2);

    return *this;
}


const CubicMatrix& CMOpenMP::dotBeta(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3)
{
    if (cm1.n() < 1 || cm2.n() < 1 || cm3.n() < 1 ||
            cm1.n() != cm2.n() || cm2.n() != cm3.n()) {
        cleanup();
        return *this;
    }

    resize(cm1.n());

    for (size_t i1 = 0; i1 < n_; i1++) {
        for (size_t j2 = 0; j2 < n_; j2++) {
            for (size_t k3 = 0; k3 < n_; k3++) {
                double sum = 0.0;

                for (size_t s = 0; s < n_; s++) {
                    size_t j1 = s, k2 = s, i3 = s;

                    for (size_t t = 0; t < n_; t++) {
                        size_t k1 = t, i2 = t, j3 = t;

                        sum += cm1.matrix()[i1][j1][k1] *
                               cm2.matrix()[i2][j2][k2] *
                               cm3.matrix()[i3][j3][k3];
                    }
                }

                matrix_[i1][j2][k3] = sum;
            }
        }
    }

    nflop_ += n_* n_ * (n_* n_ + n_* n_* n_* 2);

    return *this;
}


