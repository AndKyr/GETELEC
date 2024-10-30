#include "transmissionCalculator.h"
#include <vector>
#include <chrono>

std::vector<double> linspace(double start, double end, int n) {
    std::vector<double> result(n);
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


int main(){

    ModifiedSNBarrier tunnelFunction;
    TransmissionCalculator calculator =  TransmissionCalculator(&tunnelFunction, 3);

    int Nruns = 1.e4;
    vector<double> energies = linspace(1., -10., Nruns);
    
    auto start = std::chrono::high_resolution_clock::now();

    for (auto& energy : energies){
        tunnelFunction.setEnergy(energy);
        calculator.updateKappaAtLimits();
        calculator.solveDifferentialSystem();
        double D = calculator.transmissionCoefficient();
        // cout << " D = " << D << endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time per D calculation: " << duration.count() / (double) Nruns << " us" << std::endl;

}