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

    cm2.T();
    cm3.T().T();

#pragma omp parallel for
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


const CubicMatrix& CMOpenMP::dotBeta(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3)
{
    if (cm1.n() < 1 || cm2.n() < 1 || cm3.n() < 1 ||
            cm1.n() != cm2.n() || cm2.n() != cm3.n()) {
        cleanup();
        return *this;
    }

    resize(cm1.n());

    cm2.T();
    cm3.T().T();

    const size_t nbb1 = (n_ + bsize1_ - 1) / bsize1_;
    const size_t nbb2 = (n_ + bsize2_ - 1) / bsize2_;
    const size_t nbb3 = (n_ + bsize3_ - 1) / bsize3_;

    // Block loops
#pragma omp parallel for
    for (size_t b = 0; b < nbb1 * nbb2 * nbb3; b++) {
        size_t b1 = (b / (nbb2 * nbb3));
        size_t b2 = (b % (nbb2 * nbb3)) / nbb3;
        size_t b3 = (b % (nbb2 * nbb3)) % nbb3;

        size_t bMin1 = b1 * bsize1_;
        size_t bMin2 = b2 * bsize2_;
        size_t bMin3 = b3 * bsize3_;

        size_t bMax1 = min(bMin1 + bsize1_, n_);
        size_t bMax2 = min(bMin2 + bsize2_, n_);
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
    } // End block loops

    nflop_ += n_* n_ * (n_* n_ + n_* n_* n_* 2);

    return *this;
}


void CMOpenMP::recenter3D(const double value)
{
    const size_t bsize = 16;
    const size_t nbb = (n_ + bsize - 1) / bsize;
    const size_t nbb2 = nbb * nbb, nbb3 = nbb2 * nbb;

#pragma omp parallel
    {
#pragma omp for
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n32_; k++) {
                reductionAxisI_[j][k] = 0.0;
                reductionAxisJ_[j][k] = 0.0;
                reductionAxisK_[j][k] = 0.0;
            }
        }

#pragma omp for
        for (size_t j = 0; j < n_; j++) {
            for (size_t imin = 0; imin < n_; imin += bsize) {
                size_t imax = min(imin + bsize, n_);

                for (size_t kmin = 0; kmin < n_; kmin += bsize) {
                    size_t kmax = min(kmin + bsize, n_);

                    for (size_t i = imin; i < imax; i++) {
                        for (size_t k = kmin; k < kmax; k++) {
                            reductionAxisI_[j][k] += matrix_[i][j][k];
                        }
                    }
                }
            }
        }

#pragma omp for
        for (size_t k = 0; k < n_; k++) {
            for (size_t jmin = 0; jmin < n_; jmin += bsize) {
                size_t jmax = min(jmin + bsize, n_);

                for (size_t imin = 0; imin < n_; imin += bsize) {
                    size_t imax = min(imin + bsize, n_);

                    for (size_t j = jmin; j < jmax; j++) {
                        for (size_t i = imin; i < imax; i++) {
                            reductionAxisJ_[k][i] += matrix_[i][j][k];
                        }
                    }
                }
            }
        }

#pragma omp for
        for (size_t i = 0; i < n_; i++) {
            for (size_t kmin = 0; kmin < n_; kmin += bsize) {
                size_t kmax = min(kmin + bsize, n_);

                for (size_t jmin = 0; jmin < n_; jmin += bsize) {
                    size_t jmax = min(jmin + bsize, n_);

                    for (size_t k = kmin; k < kmax; k++) {
                        for (size_t j = jmin; j < jmax; j++) {
                            reductionAxisK_[i][j] += matrix_[i][j][k];
                        }
                    }
                }
            }
        }

#pragma omp for
        for (size_t b = 0; b < nbb3; b++) {
            size_t bi = (b / nbb2);
            size_t bj = (b % nbb2) / nbb;
            size_t bk = (b % nbb2) % nbb;

            size_t imin = bi * bsize;
            size_t jmin = bj * bsize;
            size_t kmin = bk * bsize;

            size_t imax = min(imin + bsize, n_);
            size_t jmax = min(jmin + bsize, n_);
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

    nflop_ += n_ * n_ * 3 + n_ * n_ * n_ * (3 + 5);
}


CubicMatrix& CMOpenMP::T()
{
    const size_t bsize = 8;
    const size_t nbb = (n_ + bsize - 1) / bsize;
    const size_t nbb2 = nbb * nbb;
    const size_t nbb3 = nbb * nbb * nbb;
    size_t nbBl = 0;
    size_t *bIs = new size_t[(nbb3 - nbb) / 3 + nbb];
    size_t *bJs = new size_t[(nbb3 - nbb) / 3 + nbb];
    size_t *bKs = new size_t[(nbb3 - nbb) / 3 + nbb];

    for (size_t bi = 0; bi < nbb; bi++) {
        // Main diagonal
        bIs[nbBl] = bi;
        bJs[nbBl] = bi;
        bKs[nbBl] = bi;
        nbBl++;

        for (size_t bj = bi; bj < nbb; bj++) {
            for (size_t bk = bi + 1; bk < nbb; bk++) {
                size_t m = (bi + bj + bk) % 3;

                if (m == 0) { bIs[nbBl] = bi; bJs[nbBl] = bj; bKs[nbBl] = bk; }
                if (m == 1) { bIs[nbBl] = bj; bJs[nbBl] = bk; bKs[nbBl] = bi; }
                if (m == 2) { bIs[nbBl] = bk; bJs[nbBl] = bi; bKs[nbBl] = bj; }
                nbBl++;
            }
        }
    }

#pragma omp parallel for
    for (size_t b = 0; b < nbBl; b++) {
        size_t bi = bIs[b];
        size_t bj = bJs[b];
        size_t bk = bKs[b];

        size_t imin = bi * bsize;
        size_t jmin = bj * bsize;
        size_t kmin = bk * bsize;

        size_t imax = min(imin + bsize, n_);
        size_t jmax = min(jmin + bsize, n_);
        size_t kmax = min(kmin + bsize, n_);

        if ((bi != bj) || (bj != bk)) {
            for (size_t i = imin; i < imax; i++) {
                for (size_t j = jmin; j < jmax; j++) {
                    for (size_t k = kmin; k < kmax; k++) {
                        double matrixIJK = matrix_[i][j][k];
                        matrix_[i][j][k] = matrix_[k][i][j];
                        matrix_[k][i][j] = matrix_[j][k][i];
                        matrix_[j][k][i] = matrixIJK;
                    }
                }
            }
        }
        else {
            for (size_t i = imin; i < imax; i++) {
                for (size_t j = i; j < jmax; j++) {
                    for (size_t k = i+1; k < kmax; k++) {
                        double matrixIJK = matrix_[i][j][k];
                        matrix_[i][j][k] = matrix_[k][i][j];
                        matrix_[k][i][j] = matrix_[j][k][i];
                        matrix_[j][k][i] = matrixIJK;
                    }
                }
            }
        }
    }

    delete [] bKs;
    delete [] bJs;
    delete [] bIs;

    nflop_ += (n_* n_* n_ - n_) / 3;

    return *this;
}


