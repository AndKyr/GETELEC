#include <gtest/gtest.h>
#include "Utilities.h"
#include "ODESolver.h"
#include "TransmissionSolver.h"
#include "TransmissionInterpolator.h"
#include "TunnelingFunction.h"
#include "BandEmitter.h"

#include <random>
#include <iostream>
// Test for TransmissionSolver:: check that the transmission coefficient for the default barrier is calculated correctly
TEST(TransmissionSolverTest, DefaultValueTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(8.0);
    double transmission = solver.calculateTransmissionCoefficientForEnergy(-4.5);
    EXPECT_NEAR(transmission, 0.00066697781489034457, 1.e-10);
}

// Test for TransmissionInterpolator:: check that the inteprolated and calculated values are close
TEST(TransmissionInterpolatorTest, evaluationTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(12.);
    mt19937 generator(1987);

    double testEnergy = Utilities::getUniformRandomDouble(-7., 0., generator);
    double calculatedValue = solver.calculateTransmissionCoefficientForEnergy(-4.5 + testEnergy);

    TransmissionInterpolator interpolator(solver);
    interpolator.initialize(-7., 0.0, 8);
    interpolator.refineToTolerance();
    double interpolatedValue = interpolator.evaluate(testEnergy);
    
    EXPECT_NEAR(interpolatedValue, calculatedValue, interpolator.calculateTolerance(testEnergy, log(calculatedValue))); // Transmission coefficient is expected to be less than 1
}

TEST(BandEmitterTest, DefaultValueTest){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    BandEmitter emitter(solver);
    emitter.calculateCurrentDensityAndSpectra();
    double currentDensity = emitter.getCurrentDensity();
    double nottinghamHeat = emitter.getNottinghamHeat();
    EXPECT_NEAR(3.676829392936353e-09, currentDensity, 1.e-12);
    EXPECT_NEAR(-7.5900179516309179e-10, nottinghamHeat, 1.e-13);
}

TEST(BandEmitterTest, CurrentDensityMethodComparison){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    BandEmitter emitter(solver);
    mt19937 generator(1987);
    emitter.setGenerator(&generator);
    barrier.setGenerator(&generator);

    for (int i = 0; i < 100; i++){
        barrier.setRandomParameters();
        emitter.setRandomParameters();
        emitter.updateBarrier();
        emitter.calculateCurrentDensityAndSpectra();
        double currentDensity = emitter.getCurrentDensity();
        double currentDensity2 = emitter.calcualteCurrentDensity();
        EXPECT_NEAR(currentDensity, currentDensity2, 100*emitter.getToleranceForValue(currentDensity));
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(catch_exceptions) = false;

    return RUN_ALL_TESTS();
}
