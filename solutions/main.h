#ifndef MAIN_H
#define MAIN_H


#include "cmvectorized.h"
#include "cmopenmp.h"

#include <string>


class CMArgs
{
public:
    CMArgs();
    ~CMArgs();

    bool parse(int argc, char* argv[]);
    void printArgs() const;
    void printCSV() const;

    void usage(char * cmd) const;

public:
    std::string algo;
    bool beta;

    size_t n;
    size_t ksize;
    size_t n_runs;

    size_t bsize1;  // Block size 1, bsize1 <= n
    size_t bsize2;  // Block size 2, bsize2 <= n
    size_t bsize3;  // Block size 3, bsize3 <= n

    size_t psize1;  // Prism size 1, psize1 <= bsize1
    size_t psize2;  // Prism size 2, psize2 <= bsize2
    size_t psize3;  // Prism size 3, psize3 <= bsize3

};


#endif // MAIN_H
