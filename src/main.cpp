#include "TransmissionSolver.h"
#include "BandEmitter.h"
#include <vector>
#include <iostream>
#include <chrono>


int main(){

    ModifiedSNBarrier tunnelFunction;
    TransmissionSolver calculator =  TransmissionSolver(&tunnelFunction);

    BandEmitter emitter = BandEmitter(&calculator);

    auto start = std::chrono::high_resolution_clock::now();


    emitter.calculateCurrentDensityAndSpectra();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time for current density: " << duration.count() << " us. Number of Calls: " << calculator.getNumberOfCalls() << std::endl;

    emitter.writeSolution("spectra.dat");

}