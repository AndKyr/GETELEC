#include "transmissionCalculator.h"
#include <vector>
int main(){

    ModifiedSNBarrier tunnelFunction;
    TransmissionCalculator calculator =  TransmissionCalculator(&tunnelFunction);

    calculator.solveDifferentialSystem();
    calculator.writeSolution();
}