#include <gtest/gtest.h>
#include "Utilities.h"
#include "ODESolver.h"
#include "TransmissionSolver.h"
#include "TransmissionInterpolator.h"
#include "TunnelingFunction.h"

#include <random>
#include <iostream>
// Test for TransmissionSolver::setEnergy
TEST(TransmissionSolverTest, DefaultValueTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(8.0);
    double transmission = solver.calculateTransmissionCoefficientForEnergy(-4.5);
    EXPECT_NEAR(transmission, 0.00066697781489034457, 1.e-10);
}

// Test for TransmissionInterpolator::calculateYforX
TEST(TransmissionInterpolatorTest, evaluationTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(12.);
    double testEnergy = Utilities::getUniformRandomDouble(-7., 0.);
    double calculatedValue = solver.calculateTransmissionCoefficientForEnergy(-4.5 + testEnergy);

    TransmissionInterpolator interpolator(solver);
    interpolator.initialize(-7., 0.0, 8);
    interpolator.refineToTolerance();
    double interpolatedValue = interpolator.evaluate(testEnergy);
    
    EXPECT_NEAR(interpolatedValue, calculatedValue, interpolator.calculateTolerance(testEnergy, log(calculatedValue))); // Transmission coefficient is expected to be less than 1
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
