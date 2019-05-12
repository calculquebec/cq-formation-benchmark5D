#ifndef CUBICMATRIX_H
#define CUBICMATRIX_H

#include <cstddef>


#define CM_MIN_N 2
#define CM_MAX_N 960


class CubicMatrix
{
public:
    CubicMatrix();
    virtual ~CubicMatrix();

    void cleanup();

    double dist(CubicMatrix &cm);

    virtual const CubicMatrix& dot(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3);
    virtual const CubicMatrix& dotBeta(CubicMatrix &cm1, CubicMatrix &cm2, CubicMatrix &cm3);

    inline size_t getNflop() const { return nflop_; }

    void head(const size_t n = 3);

    void initCount();
    void initRandomOddInt(unsigned int seed);
    virtual void initReduction(const double value = 0.0);

    double mean();

    inline double const* const* const* matrix() const { return matrix_; }
    inline size_t n() const { return n_; }

    bool operator==(const CubicMatrix &cm) const;

    virtual void recenter3D(const double value = 0.0);
    void reductionProd();
    virtual void reductionSum();

    inline void resetNflop() { nflop_ = 0; }
    void resize(const size_t n = 3);

    inline void setKSize(const size_t ksize) { ksize_ = ksize; }
    void setBlockSize(const size_t bsize1, const size_t bsize2, const size_t bsize3);
    void setPrismSize(const size_t psize1, const size_t psize2, const size_t psize3);

    virtual CubicMatrix& T();

    double varP();

protected:
    void printConfig(const char * algo) const;

protected:
    size_t    nflop_;

    size_t    n_;       // n_ <= n32_
    size_t    n32_;     // n_ rounded to the next multiple of 32
    size_t    ksize_;   // Block size in K axis

    size_t    bsize1_;  // Block size 1, bsize1_ <= n_
    size_t    bsize2_;  // Block size 2, bsize2_ <= n_
    size_t    bsize3_;  // Block size 3, bsize3_ <= n_

    size_t    psize1_;  // Prism size 1, psize1_ <= bsize1_
    size_t    psize2_;  // Prism size 2, psize2_ <= bsize2_
    size_t    psize3_;  // Prism size 3, psize3_ <= bsize3_

    double*** matrix_;

    double** reductionAxisI_;
    double** reductionAxisJ_;
    double** reductionAxisK_;

private:
    double*   data_;
    double**  vctrs_;
    double*** matbuf_;

    double* dataI_;
    double* dataJ_;
    double* dataK_;

    double** reductBufI_;
    double** reductBufJ_;
    double** reductBufK_;

};


class CMFactory
{
public:
    virtual ~CMFactory() {}
    virtual CubicMatrix* create() const { return new CubicMatrix(); }
};


#endif // CUBICMATRIX_H

