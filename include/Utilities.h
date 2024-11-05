#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cmath>
#include <vector>
#include <limits>

using namespace std;

static constexpr struct PhysicalConstants{
    double hbarSqrOver2m = 3.80998212e-2;
    double kConstant = 1./hbarSqrOver2m;
    double imageChargeConstant = .359991137;
    double BoltzmannConstant = 8.617333262e-5;
    double SommerfeldConstant = 1.618311e-4;
    double exponentLimit = - 0.5 * log(numeric_limits<double>::epsilon());

} CONSTANTS;

struct Utilities{
    static vector<double> linspace(double start, double end, int n);

    static double fermiDiracFunction(double energy, double kT);

    static double logFermiDiracFunction(double energy, double kT);
};

#endif //UTILITIES_H_
