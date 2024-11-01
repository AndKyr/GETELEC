#include "Utilities.h"

vector<double> Utilities::linspace(double start, double end, int n) {
    vector<double> result(n);
    if (n == 1) {
        result[0] = start;
        return result;
    }

    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

double Utilities::fermiDiracFunction(double energy, double kT){
    if (energy > CONSTANTS.exponentLimit * kT)
        return exp(-energy / kT);
    else if(energy < - CONSTANTS.exponentLimit * kT)
        return 1. - exp(energy / kT);
    else
        return 1. / (1. + exp(energy / kT));
}

double Utilities::logFermiDiracFunction(double energy, double kT){
    if (energy > CONSTANTS.exponentLimit * kT)
        return exp(-energy / kT);
    else if(energy < - CONSTANTS.exponentLimit * kT)
        return -energy / kT + exp(energy / kT);
    else
        return log(1. + exp(-energy / kT));
}