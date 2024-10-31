#include "transmissionCalculator.h"
#include <vector>
#include <iostream>
#include <chrono>



int main(){

    ModifiedSNBarrier tunnelFunction;
    TransmissionSolver calculator =  TransmissionSolver(&tunnelFunction);

    int Nruns = 100;
    vector<double> energies = Utilities::linspace(1., -10., Nruns);
    
    auto start = std::chrono::high_resolution_clock::now();

    for (auto& energy : energies){
        tunnelFunction.setEnergy(energy);
        calculator.updateKappaAtLimits();
        calculator.solveNoSave();
        double D = calculator.transmissionCoefficient();
        cout << " D = " << D << endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time per D calculation: " << duration.count() / (double) Nruns << " us" << std::endl;

}