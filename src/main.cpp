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

    TransmissionInterpolator interpolator = TransmissionInterpolator(calculator, 4.5, 0.025, -5.4, 3.5, 4);

    interpolator.refineToTolerance(20);
    interpolator.writeSplineNodes();
    ofstream outFile("interpolatedTransmission.dat", ios::out);        

    for (double x = -5.4; x < 3.5; x+=0.001){
        double D = calculator.calculateTransmissionCoefficientForEnergy(x - 4.5);
        outFile << x << " " << D << " " << interpolator.evaluate(x) << " " << abs(D - interpolator.evaluate(x)) << endl;
    }
}