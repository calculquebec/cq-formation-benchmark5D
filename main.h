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
    size_t n_runs;

};


#endif // MAIN_H
