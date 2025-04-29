#include "TransmissionSolver.h"
#include "BandEmitter.h"
#include <vector>
#include <iostream>
#include <chrono>

using namespace std;
using namespace getelec;

int main(){

    ModifiedSNBarrier tunnelFunction;


    tunnelFunction.setField(3.);

    TransmissionSolver calculator =  TransmissionSolver(&tunnelFunction);

    BandEmitter emitter = BandEmitter(calculator);


    emitter.integrateTotalEnergyDistributionODEAndSaveSpectra();
    cout << emitter.getCurrentDensityODE() << endl;

    auto fields = Utilities::linspace(1., 10, 64);
    auto kTs = Utilities::linspace(.01, .3, 64);
    auto start = std::chrono::high_resolution_clock::now();
    for (double kT : kTs){
        emitter.setParameters(4.5, kT);
        emitter.integrateTotalEnergyDistributionODEAndSaveSpectra();
        cout << kT  << " " << emitter.getCurrentDensityODE() << endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time for current density: " << duration.count() / fields.size() << " us. Number of Calls: " << calculator.getNumberOfCalls() << std::endl;

    emitter.writeSolution("spectra.dat");

    emitter.writePlottingData();

    return 0;


    // TransmissionInterpolator interpolator = TransmissionInterpolator(calculator, 4.5, 0.025, 1.e-12, 1.e-4);
    // interpolator.initialize(-5.4, 3.5, 4);
    // interpolator.refineToTolerance();
    // interpolator.writeSplineNodes();
    // ofstream outFile("interpolatedTransmission.dat", ios::out);        
    // outFile << " E D_calc D_interp error" << endl;
    // for (double x = -5.4; x < 3.5; x+=0.001){
    //     double D = calculator.calculateTransmissionCoefficientForEnergy(x - 4.5);
    //     outFile << x << " " << D << " " << interpolator.evaluate(x) << " " << abs(D - interpolator.evaluate(x)) << endl;
    // }
}
