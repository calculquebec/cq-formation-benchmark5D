#include "cmvectorized.h"

#include <algorithm>
#include <iostream>
#include <stdlib.h>


using namespace std;


CMVectorized::CMVectorized(): CubicMatrix()
{
}


CMVectorized::~CMVectorized()
{
}


const CubicMatrix& CMVectorized::dot(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3)
{
    if (cm1.n() < 1 || cm2.n() < 1 || cm3.n() < 1 ||
            cm1.n() != cm2.n() || cm2.n() != cm3.n()) {
        cleanup();
        return *this;
    }

    resize(cm1.n());

    cm2.T();
    cm3.T().T();

    for (size_t i1 = 0; i1 < n_; i1++) {
        for (size_t i2 = 0; i2 < n_; i2++) {
            for (size_t i3 = 0; i3 < n_; i3++) {
                double sum = 0.0;

                for (size_t j = 0; j < n_; j++) {
                    for (size_t k = 0; k < n_; k++) {
                        sum += cm1.matrix()[i1][j][k] *
                               cm2.matrix()[i2][j][k] *
                               cm3.matrix()[i3][j][k];
                    }
                }

                matrix_[i1][i2][i3] = sum;
            }
        }
    }

    nflop_ += n_* n_ * (n_* n_ + n_* n_* n_* 2);

    return *this;
}


const CubicMatrix& CMVectorized::dotBeta(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3)
{
    if (cm1.n() < 1 || cm2.n() < 1 || cm3.n() < 1 ||
            cm1.n() != cm2.n() || cm2.n() != cm3.n()) {
        cleanup();
        return *this;
    }

    resize(cm1.n());

    cm2.T();
    cm3.T().T();

    // Block loops
    for (size_t bMin1 = 0; bMin1 < n_; bMin1 += bsize1_) {
      size_t bMax1 = min(bMin1 + bsize1_, n_);
     for (size_t bMin2 = 0; bMin2 < n_; bMin2 += bsize2_) {
       size_t bMax2 = min(bMin2 + bsize2_, n_);
      for (size_t bMin3 = 0; bMin3 < n_; bMin3 += bsize3_) {
        size_t bMax3 = min(bMin3 + bsize3_, n_);

        // Initialize the current block
        for (size_t i1 = bMin1; i1 < bMax1; i1++) {
         for (size_t i2 = bMin2; i2 < bMax2; i2++) {
          for (size_t i3 = bMin3; i3 < bMax3; i3++) {
            matrix_[i1][i2][i3] = 0.0;
        }}}

        // Main dot-product loops
        for (size_t j = 0; j < n_; j++) {
          for (size_t kmin = 0; kmin < n_; kmin += ksize_) {
            size_t kmax = min(kmin + ksize_, n_);

            // Prism loops
            for (size_t pmin1 = bMin1; pmin1 < bMax1; pmin1 += psize1_) {
              size_t pmax1 = min(pmin1 + psize1_, bMax1);
             for (size_t pmin2 = bMin2; pmin2 < bMax2; pmin2 += psize2_) {
               size_t pmax2 = min(pmin2 + psize2_, bMax2);
              for (size_t pmin3 = bMin3; pmin3 < bMax3; pmin3 += psize3_) {
                size_t pmax3 = min(pmin3 + psize3_, bMax3);

                // Partial sum of current prism
                for (size_t i1 = pmin1; i1 < pmax1; i1++) {
                 for (size_t i2 = pmin2; i2 < pmax2; i2++) {
                  for (size_t i3 = pmin3; i3 < pmax3; i3++) {
                    double sum = 0.0;

                    for (size_t k = kmin; k < kmax; k++) {
                        sum += cm1.matrix()[i1][j][k] *
                               cm2.matrix()[i2][j][k] *
                               cm3.matrix()[i3][j][k];
                    }

                    matrix_[i1][i2][i3] += sum;
                }}} // End partial sum of current prism
            }}} // End prism loops
        } } // End dot-product loops
    }}} // End block loops

    nflop_ += n_* n_ * (n_* n_ + n_* n_* n_* 2);

    return *this;
}


void CMVectorized::recenter3D(const double value)
{
    reductionSum();

    const size_t bsize = 8;

    for (size_t imin = 0; imin < n_; imin += bsize) {
        size_t imax = min(imin + bsize, n_);

        for (size_t jmin = 0; jmin < n_; jmin += bsize) {
            size_t jmax = min(jmin + bsize, n_);

            for (size_t kmin = 0; kmin < n_; kmin += bsize) {
                size_t kmax = min(kmin + bsize, n_);

                for (size_t i = imin; i < imax; i++) {
                    for (size_t j = jmin; j < jmax; j++) {
                        for (size_t k = kmin; k < kmax; k++) {
                            matrix_[i][j][k] += value -
                                (   reductionAxisI_[j][k] +
                                    reductionAxisJ_[k][i] +
                                    reductionAxisK_[i][j]   ) / (3 * n_);
                        }
                    }
                }
            }
        }
    }

    nflop_ += n_ * n_ * n_ * 5;
}


void CMVectorized::reductionSum()
{
    initReduction(0.0);

    const size_t bsize = 8;

    for (size_t imin = 0; imin < n_; imin += bsize) {
        size_t imax = min(imin + bsize, n_);

        for (size_t jmin = 0; jmin < n_; jmin += bsize) {
            size_t jmax = min(jmin + bsize, n_);

            for (size_t kmin = 0; kmin < n_; kmin += bsize) {
                size_t kmax = min(kmin + bsize, n_);

                for (size_t i = imin; i < imax; i++) {
                    for (size_t j = jmin; j < jmax; j++) {
                        for (size_t k = kmin; k < kmax; k++) {
                            reductionAxisI_[j][k] += matrix_[i][j][k];
                            reductionAxisJ_[k][i] += matrix_[i][j][k];
                            reductionAxisK_[i][j] += matrix_[i][j][k];
                        }
                    }
                }
            }
        }
    }

    nflop_ += n_ * n_ * n_ * 3;
}


