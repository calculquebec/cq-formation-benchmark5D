#include "main.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <time.h>


using namespace std;


#define NB_CM 5

#define TAG_NORESULT  0
#define TAG_NEWRUN    1
#define TAG_PASSED    7
#define TAG_FAILED   13
#define TAG_NO_RUN   99


///////////////////////////////////////
//  b5D_Wtime()
///////////////////////////////////////
double b5D_Wtime()
{
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return tp.tv_sec + tp.tv_nsec / 1E9;
}


///////////////////////////////////////
//  computeRun()
///////////////////////////////////////
size_t computeRun(CubicMatrix* cm[NB_CM], const CMArgs &cmargs, const size_t runID = 0)
{
    size_t totalFlop = 0;

    // Rotation between all 5 CubicMatrices
    for (size_t cmID = 0; cmID < NB_CM; cmID++) {
        // Use a different target CubicMatrix
        size_t mCurr = (cmID + 0) % NB_CM;
        size_t m1    = (cmID + 1) % NB_CM;
        size_t m2    = (cmID + 2) % NB_CM;
        size_t m3    = (cmID + 3) % NB_CM;
        size_t mPrev = (cmID + NB_CM - 1) % NB_CM;

        // Reset counters
        cm[m1]->resetNflop();
        cm[m2]->resetNflop();
        cm[m3]->resetNflop();
        cm[mCurr]->resetNflop();

        // Reset m1, m2 and m3 matrices
        cm[m1]->initRandomOddInt(3 * runID + 1);
        cm[m2]->initRandomOddInt(3 * runID + 2);
        cm[m3]->initRandomOddInt(3 * runID + 3);

        // 3D dot-product
        if (!cmargs.beta) {
            cm[mCurr]->dot(    *cm[m1], *cm[m2], *cm[m3]);
        } else {
            cm[mCurr]->dotBeta(*cm[m1], *cm[m2], *cm[m3]);
        }

        // Recenter values
        for (size_t iter = 0; iter < 8; iter++) {
            cm[mCurr]->recenter3D(0.0);
        }

        // If not first iteration
        if (mCurr > 0) {
            // cm[mCurr]->recenter3D(0.0); // Uncomment to insert errors
            double error = cm[mCurr]->dist(*cm[mPrev]);

            cout << "Dist-" << runID + 1 << "."
                 << mPrev + 1 << "-" << mCurr + 1 << ","
                 << error << ","
                 << setprecision(12) << cm[mCurr]->varP() << endl;

            if (error < 0.0 || 0.0 < error) return 0;
        }

        totalFlop += cm[m1]->getNflop();
        totalFlop += cm[m2]->getNflop();
        totalFlop += cm[m3]->getNflop();
        totalFlop += cm[mCurr]->getNflop();
    }

    return totalFlop;
}


///////////////////////////////////////
//  main()
///////////////////////////////////////
int main(int argc, char* argv[])
{
    CMArgs cmargs;

    if (!cmargs.parse(argc, argv)) {
        exit(1);
    }

    CMFactory* cmf = NULL;

    if      (cmargs.algo == "vect") { cmf = new CMVectorizedFactory(); }
    else if (cmargs.algo == "opmp") { cmf = new CMOpenMPFactory(); }
    else /**/ { cmargs.algo = "base"; cmf = new CMFactory(); }

    CubicMatrix* cm[NB_CM];

    for (size_t cmID = 0; cmID < NB_CM; cmID++) {
        cm[cmID] = cmf->create();

        cm[cmID]->resize(cmargs.n);
    }

    bool passed = true;
    size_t totalFlop = 0;
    double timeStart = b5D_Wtime();

    cout << "DistID,Distance,VarP" << endl;

    for (size_t runID = 0; runID < cmargs.n_runs; runID++) {
        size_t moreFlop = computeRun(cm, cmargs, runID);

        if (moreFlop == 0) {
            passed = false;
            break;
        }

        totalFlop += moreFlop;
    }

    double seconds = b5D_Wtime() -  timeStart;
    double teraFlops = totalFlop / (1000.0 * 1000.0 * 1000.0 * 1000.0) / seconds;

    cout << "     Final Nb of FLOP: " << totalFlop << endl;
    cout << "  Elapsed time (sec.): " << seconds << endl;
    cout << "         Global speed: " << teraFlops << " TFlop/s" << endl;

    cout << "Result=" << (passed ? "passed" : "failed") << "," << teraFlops << ",";
    cmargs.printCSV();
    cout << endl;

    for (size_t cmID = 0; cmID < NB_CM; cmID++) {
        delete cm[cmID];
    }

    return 0;
}


///////////////////////////////////////
//  class CMArgs
///////////////////////////////////////
CMArgs::CMArgs():
    algo(""), beta(false),
    n(0), n_runs(1)
{
}


CMArgs::~CMArgs()
{
}


bool CMArgs::parse(int argc, char* argv[])
{
    int i = 0;

    while (++i < argc) {
        if (strcmp(argv[i], "-h") == 0) { usage(argv[0]); return false; }

        if (strcmp(argv[i], "-a") == 0 && (i + 1) < argc) algo = argv[++i];
        if (strcmp(argv[i], "-b") == 0) beta = true;

        if (strcmp(argv[i], "-n") == 0 && (i + 1) < argc) n = atoi(argv[++i]);
        if (strcmp(argv[i], "-r") == 0 && (i + 1) < argc) n_runs = atoi(argv[++i]);
    }

    if (algo.length() == 0) {
        cerr << "Error: please specify an AlgoName ('base', 'vect' or 'opmp')" << endl;
        return false;
    }

    if (n < CM_MIN_N || CM_MAX_N < n) {
        cerr << "Error: N is mandatory and must be inside range ["
             << CM_MIN_N << ", " << CM_MAX_N << "]" << endl;
        return false;
    }

    if (n_runs < 1) {
        cerr << "Error: Number of runs must be at least 1" << endl;
        return false;
    }

    return true;
}


void CMArgs::printArgs() const
{
    cout << " Algo: " << algo << (beta ? "(beta)":"") << endl;
    cout << "    N: " << n << endl;

    cout << "nruns: " << n_runs << endl;
}


void CMArgs::printCSV() const
{
    cout << algo << (beta ? "(beta)":"") << "," << n << ",";
    cout << n_runs;
}


void CMArgs::usage(char *cmd) const
{
    cout << "Usage: " << cmd << " ";
    cout << "-a AlgoName [-b] -n N [-b1 size1 -b2 size2 -b3 size3] \\" << endl;
    cout << "\t[-p1 size4 -p2 size5 -p3 size6] [-k ksize] [-r R]" << endl;
    cout << endl;

    cout << "-a AlgoName # Either 'base', 'vect' or 'opmp'" << endl;
    cout << "-b          # Active beta code, or optimized version" << endl;
    cout << "-n N        # Cubic matrix size, min=" << CM_MIN_N << ", max=" << CM_MAX_N << endl;
    cout << "-r R        # Number of different runs (default = 1)" << endl;
}


