#include "cubicmatrix.h"


#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>


using namespace std;


CubicMatrix::CubicMatrix():
    nflop_(0),
    n_(0), n32_(0), ksize_(0),
    bsize1_(0), bsize2_(0), bsize3_(0),
    psize1_(0), psize2_(0), psize3_(0),
    data_(NULL), vctrs_(NULL), matbuf_(NULL), matrix_(NULL),
    dataI_(NULL), reductBufI_(NULL), reductionAxisI_(NULL),
    dataJ_(NULL), reductBufJ_(NULL), reductionAxisJ_(NULL),
    dataK_(NULL), reductBufK_(NULL), reductionAxisK_(NULL)
{
}


CubicMatrix::~CubicMatrix()
{
    cleanup();
}


void CubicMatrix::cleanup()
{
    reductionAxisK_ = NULL;
    reductionAxisJ_ = NULL;
    reductionAxisI_ = NULL;

    delete [] reductBufK_; reductBufK_ = NULL;
    delete [] reductBufJ_; reductBufJ_ = NULL;
    delete [] reductBufI_; reductBufI_ = NULL;

    delete [] dataK_; dataK_ = NULL;
    delete [] dataJ_; dataJ_ = NULL;
    delete [] dataI_; dataI_ = NULL;

    delete [] matbuf_; matbuf_ = NULL; matrix_ = NULL;
    delete [] vctrs_;  vctrs_ = NULL;
    delete [] data_;   data_ = NULL;
    n32_ = 0;
    n_ = 0;
}


double CubicMatrix::dist(CubicMatrix &cm)
{
    if (n_ != cm.n_) return -1.0;

    double sumDist2 = 0.0;

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                double diff = cm.matrix_[i][j][k] - matrix_[i][j][k];
                sumDist2 += diff * diff;
            }
        }
    }

    nflop_ += n_ * n_ * n_ * 3;

    return sqrt(sumDist2);
}



const CubicMatrix& CubicMatrix::dot(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3)
{
    if (cm1.n_ < 1 || cm2.n_ < 1 || cm3.n_ < 1 ||
            cm1.n_ != cm2.n_ || cm2.n_ != cm3.n_) {
        cleanup();
        return *this;
    }

    resize(cm1.n_);

    for (size_t i1 = 0; i1 < n_; i1++) {
        for (size_t j2 = 0; j2 < n_; j2++) {
            for (size_t k3 = 0; k3 < n_; k3++) {
                double sum = 0.0;

                for (size_t s = 0; s < n_; s++) {
                    size_t j1 = s, k2 = s, i3 = s;

                    for (size_t t = 0; t < n_; t++) {
                        size_t k1 = t, i2 = t, j3 = t;

                        sum += cm1.matrix_[i1][j1][k1] *
                               cm2.matrix_[i2][j2][k2] *
                               cm3.matrix_[i3][j3][k3];
                    }
                }

                matrix_[i1][j2][k3] = sum;
            }
        }
    }

    nflop_ += n_* n_ * (n_* n_ + n_* n_* n_ * 2);

    return *this;
}


const CubicMatrix& CubicMatrix::dotBeta(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3)
{
    if (cm1.n_ < 1 || cm2.n_ < 1 || cm3.n_ < 1 ||
            cm1.n_ != cm2.n_ || cm2.n_ != cm3.n_) {
        cleanup();
        return *this;
    }

    resize(cm1.n_);

    // Block loops
    for (size_t bMin1 = 0; bMin1 < n_; bMin1 += bsize1_) {
        size_t bMax1 = min(bMin1 + bsize1_, n_);

        for (size_t bMin2 = 0; bMin2 < n_; bMin2 += bsize2_) {
            size_t bMax2 = min(bMin2 + bsize2_, n_);

            for (size_t bMin3 = 0; bMin3 < n_; bMin3 += bsize3_) {
                size_t bMax3 = min(bMin3 + bsize3_, n_);

                // Prism loops
                for (size_t pmin1 = bMin1; pmin1 < bMax1; pmin1 += psize1_) {
                    size_t pmax1 = min(pmin1 + psize1_, bMax1);

                    for (size_t pmin2 = bMin2; pmin2 < bMax2; pmin2 += psize2_) {
                        size_t pmax2 = min(pmin2 + psize2_, bMax2);

                        for (size_t pmin3 = bMin3; pmin3 < bMax3; pmin3 += psize3_) {
                            size_t pmax3 = min(pmin3 + psize3_, bMax3);

                            // Main dot-product loops
                            for (size_t i1 = pmin1; i1 < pmax1; i1++) {
                                for (size_t j2 = pmin2; j2 < pmax2; j2++) {
                                    for (size_t k3 = pmin3; k3 < pmax3; k3++) {
                                        double sum = 0.0;

                                        for (size_t s = 0; s < n_; s++) {
                                            size_t j1 = s, k2 = s, i3 = s;

                                            for (size_t t = 0; t < n_; t++) {
                                                size_t k1 = t, i2 = t, j3 = t;

                                                sum += cm1.matrix_[i1][j1][k1] *
                                                       cm2.matrix_[i2][j2][k2] *
                                                       cm3.matrix_[i3][j3][k3];
                                            }
                                        }

                                        matrix_[i1][j2][k3] = sum;
                                    }
                                }
                            } // End dot-product loops
                        }
                    }
                } // End prism loops
            }
        }
    } // End block loops

    nflop_ += n_* n_ * (n_* n_ + n_* n_* n_ * 2);

    return *this;
}


void CubicMatrix::head(const size_t n)
{
    const size_t N = min(n_, n);

    if (N < 1) {
        cout << "Null" << endl;
        return;
    }

    for (size_t i = 0; i < N; i++) {
        for (size_t k = 0; k < N; k++) {
            for (size_t j = 0; j < N - k; j++) cout << "  ";

            for (size_t j = 0; j < N; j++) {
                cout << matrix_[i][j][N-1 - k] << " ";
            }
            cout << "->" << reductionAxisJ_[N-1 - k][i] << endl;
        }

        cout << "=>";
        for (size_t j = 0; j < N; j++) cout << reductionAxisK_[i][j] << " ";
        cout << endl;

        for (size_t j = 0; j < N-1; j++) cout << "  ";
        cout << "-----" << endl;
    }

    for (size_t k = 0; k < N; k++) {
        for (size_t j = 0; j < N-1 - k; j++) cout << "  ";
        cout << "=>";

        for (size_t j = 0; j < N; j++) {
            cout << reductionAxisI_[j][N-1 - k] << " ";
        }
        cout << endl;
    }
}


void CubicMatrix::initCount()
{
    size_t ind = 0;

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                matrix_[i][j][k] = ++ind;
            }
            for (size_t k = n_; k < n32_; k++) {
                matrix_[i][j][k] = 0.0;
            }
        }
    }

    nflop_ += n_ * n_ * n_;
}


void CubicMatrix::initRandomOddInt(unsigned int seed)
{
    srand(seed);

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                int randomValue = rand() & 0xff;
                int sign = 1 - (randomValue & 1) * 2;
                int oddValue = (randomValue | 1);
                matrix_[i][j][k] = sign * oddValue;
            }
            for (size_t k = n_; k < n32_; k++) {
                matrix_[i][j][k] = 0.0;
            }
        }
    }

    nflop_ += n_ * n_ * n_;
}


void CubicMatrix::initReduction(const double value)
{
    for (size_t s = 0; s < n_; s++) {
        for (size_t t = 0; t < n32_; t++) {
            reductionAxisI_[s][t] = value;
            reductionAxisJ_[s][t] = value;
            reductionAxisK_[s][t] = value;
        }
    }

    nflop_ += n_ * n_ * 3;
}


double CubicMatrix::mean()
{
    if (n_ < 1) return NAN;

    const size_t N = n_ * n_ * n_;
    double sum = 0.0;

    // Initialize reduction
    for (size_t j = 0; j < n_; j++) {
        for (size_t k = 0; k < n_; k++) {
            reductionAxisI_[j][k] = matrix_[0][j][k];
        }
    }

    // Reduction following i-axis
    for (size_t i = 1; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                reductionAxisI_[j][k] += matrix_[i][j][k];
            }
        }
    }

    // Reduction following j-axis
    for (size_t j = 1; j < n_; j++) {
        for (size_t k = 0; k < n_; k++) {
            reductionAxisI_[0][k] += reductionAxisI_[j][k];
        }
    }

    // Reduction following k-axis
    for (size_t k = 0; k < n_; k++) {
        sum += reductionAxisI_[0][k];
    }

    nflop_ += n_ * n_ * n_;

    return sum / N;
}


void CubicMatrix::printConfig(const char * algo) const
{
    cout << " Algo: " << algo << endl;
    cout << "    N: " << n_ << endl;
    cout << "Block: " << bsize1_ << "x" << bsize2_ << "x" << bsize3_ << endl;
    cout << "Prism: " << psize1_ << "x" << psize2_ << "x" << psize3_ << endl;
    cout << "ksize: " << ksize_ << endl;
}


bool CubicMatrix::operator==(const CubicMatrix &cm) const
{
    if (n_ != cm.n_) return false;

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                if (matrix_[i][j][k] != cm.matrix_[i][j][k]) return false;
            }
        }
    }

    return true;
}


void CubicMatrix::recenter3D(const double value)
{
    reductionSum();

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                matrix_[i][j][k] += value -
                    (reductionAxisI_[j][k] + reductionAxisJ_[k][i] + reductionAxisK_[i][j]) / (3 * n_);
            }
        }
    }

    nflop_ += n_ * n_ * n_ * 5;
}


void CubicMatrix::reductionProd()
{
    initReduction(1.0);

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                reductionAxisI_[j][k] *= matrix_[i][j][k];
                reductionAxisJ_[k][i] *= matrix_[i][j][k];
                reductionAxisK_[i][j] *= matrix_[i][j][k];
            }
        }
    }

    nflop_ += n_ * n_ * n_ * 3;
}


void CubicMatrix::reductionSum()
{
    initReduction(0.0);

    for (size_t i = 0; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                reductionAxisI_[j][k] += matrix_[i][j][k];
                reductionAxisJ_[k][i] += matrix_[i][j][k];
                reductionAxisK_[i][j] += matrix_[i][j][k];
            }
        }
    }

    nflop_ += n_ * n_ * n_ * 3;
}


void CubicMatrix::resize(const size_t n)
{
    if (n_ == n) return;

    cleanup();

    if (n > 0) {

        // Main cubic matrix

        n_ = n;
        n32_ = ((n + 31) & -32);

        data_   = new double  [32 + n * n32_ * n32_];
        vctrs_  = new double* [32 + n * n32_];
        matbuf_ = new double**[32 + n];

        size_t offsetD = 0, offsetV = 0, offsetM = 0;

        while ((size_t)(&data_[offsetD]) & 0xff) offsetD++;
        while ((size_t)(&vctrs_[offsetV]) & 0xff) offsetV++;
        while ((size_t)(&matbuf_[offsetM]) & 0xff) offsetM++;

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                vctrs_[offsetV + i * n32_ + j] = &data_[offsetD + (i * n32_ + j) * n32_ + 0];
            }
            matbuf_[offsetM + i] = &vctrs_[offsetV + i * n32_ + 0];
        }
        matrix_ = &matbuf_[offsetM + 0];

        // Reduction matrices

        dataI_ = new double[32 + n * n32_];
        dataJ_ = new double[32 + n * n32_];
        dataK_ = new double[32 + n * n32_];

        reductBufI_ = new double*[32 + n];
        reductBufJ_ = new double*[32 + n];
        reductBufK_ = new double*[32 + n];

        size_t offsetDI= 0, offsetDJ = 0, offsetDK = 0;
        size_t offsetVI= 0, offsetVJ = 0, offsetVK = 0;

        while ((size_t)(&dataI_[offsetDI]) & 0xff) offsetDI++;
        while ((size_t)(&dataJ_[offsetDJ]) & 0xff) offsetDJ++;
        while ((size_t)(&dataK_[offsetDK]) & 0xff) offsetDK++;

        while ((size_t)(&reductBufI_[offsetVI]) & 0xff) offsetVI++;
        while ((size_t)(&reductBufJ_[offsetVJ]) & 0xff) offsetVJ++;
        while ((size_t)(&reductBufK_[offsetVK]) & 0xff) offsetVK++;

        for (size_t s = 0; s < n; s++) {
            reductBufI_[offsetVI + s] = &dataI_[offsetDI + s * n32_ + 0];
            reductBufJ_[offsetVJ + s] = &dataJ_[offsetDJ + s * n32_ + 0];
            reductBufK_[offsetVK + s] = &dataK_[offsetDK + s * n32_ + 0];
        }

        reductionAxisI_ = &reductBufI_[offsetVI + 0];
        reductionAxisJ_ = &reductBufJ_[offsetVJ + 0];
        reductionAxisK_ = &reductBufK_[offsetVK + 0];
    }
}


void CubicMatrix::setBlockSize(const size_t bsize1, const size_t bsize2, const size_t bsize3)
{
    bsize1_ = bsize1;
    bsize2_ = bsize2;
    bsize3_ = bsize3;
}


void CubicMatrix::setPrismSize(const size_t psize1, const size_t psize2, const size_t psize3)
{
    psize1_ = psize1;
    psize2_ = psize2;
    psize3_ = psize3;
}


CubicMatrix& CubicMatrix::T()
{
    for (size_t i = 0; i < n_; i++) {
        for (size_t j = i; j < n_; j++) {
            for (size_t k = i+1; k < n_; k++) {
                double matrixIJK = matrix_[i][j][k];
                matrix_[i][j][k] = matrix_[k][i][j];
                matrix_[k][i][j] = matrix_[j][k][i];
                matrix_[j][k][i] = matrixIJK;
            }
        }
    }

    nflop_ += (n_*n_*n_ - n_) / 3;

    return *this;
}


double CubicMatrix::varP()
{
    if (n_ < 1) return NAN;

    const size_t N = n_ * n_ * n_;
    double sum1 = 0.0;
    double sum2 = 0.0;

    // Initialize reduction
    for (size_t j = 0; j < n_; j++) {
        for (size_t k = 0; k < n_; k++) {
            reductionAxisI_[j][k] = matrix_[0][j][k];
            reductionAxisJ_[j][k] = matrix_[0][j][k] * matrix_[0][j][k];
        }
    }

    // Reduction following i-axis
    for (size_t i = 1; i < n_; i++) {
        for (size_t j = 0; j < n_; j++) {
            for (size_t k = 0; k < n_; k++) {
                reductionAxisI_[j][k] += matrix_[i][j][k];
                reductionAxisJ_[j][k] += matrix_[i][j][k] * matrix_[i][j][k];
            }
        }
    }

    // Reduction following j-axis
    for (size_t j = 1; j < n_; j++) {
        for (size_t k = 0; k < n_; k++) {
            reductionAxisI_[0][k] += reductionAxisI_[j][k];
            reductionAxisJ_[0][k] += reductionAxisJ_[j][k];
        }
    }

    // Reduction following k-axis
    for (size_t k = 0; k < n_; k++) {
        sum1 += reductionAxisI_[0][k];
        sum2 += reductionAxisJ_[0][k];
    }

    nflop_ += n_ * n_ * n_ * 3 + 4;

    double ct_mean = sum1 / N;

    return sum2 / N - ct_mean * ct_mean;
}


