#include <gtest/gtest.h>
#include "Utilities.h"
#include "ODESolver.h"
#include "TransmissionSolver.h"
#include "TransmissionInterpolator.h"
#include "TunnelingFunction.h"
#include "BandEmitter.h"
#include "Getelec.h"

#include <random>
#include <iostream>
// Test for TransmissionSolver:: check that the transmission coefficient for the default barrier is calculated correctly
TEST(TransmissionSolverTest, DefaultValueTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(8.0);
    double transmission = solver.calculateTransmissionCoefficientForEnergy(-4.5);
    EXPECT_NEAR(transmission, 0.00066697796753640074, 1.e-10);
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
    EXPECT_NEAR(3.6745963706309422e-09, currentDensity, emitter.getToleranceForValue(currentDensity));
    EXPECT_NEAR(-7.5900179516309179e-10, nottinghamHeat, emitter.getToleranceForValue(nottinghamHeat));
}

TEST(BandEmitterTest, CurrentDensityMethodComparison){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    BandEmitter emitter(solver);
    mt19937 generator(1987);
    emitter.setGenerator(&generator);
    barrier.setGenerator(&generator);

    for (int i = 0; i < 64; i++){
        barrier.setRandomParameters();
        emitter.setRandomParameters();
        emitter.updateSolverAndInterpolator();
        emitter.calculateCurrentDensityAndNottingham();
        double currentDensity = emitter.getCurrentDensity();
        double currentDensity2 = emitter.calcualteCurrentDensity();
        EXPECT_NEAR(currentDensity, currentDensity2, 100*emitter.getToleranceForValue(currentDensity));
    }
}

TEST(ConfigTest, ConfigFileReadTest){

    Config config;
    config.bandEmitterParams.absoluteTolerance *= 1.e-1;
    config.bandEmitterParams.relativeTolerance *= 1.e-1;
    config.bandEmitterParams.maxSteps += 1;
    config.print_all_params("tempConfig.txt");
    Config config2("tempConfig.txt");
    EXPECT_EQ(config.bandEmitterParams.absoluteTolerance, config2.bandEmitterParams.absoluteTolerance);
    EXPECT_EQ(config.bandEmitterParams.relativeTolerance, config2.bandEmitterParams.relativeTolerance);
    EXPECT_EQ(config.bandEmitterParams.maxSteps, config2.bandEmitterParams.maxSteps);
    remove("tempConfig.txt");
}


TEST(GetelecObjectTest, RunParalleltest){
    Getelec getelec;
    auto fields = Utilities::linspace(2., 10., 64);
    getelec.setField(fields);
    EXPECT_NO_THROW(getelec.run());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(catch_exceptions) = false;

    return RUN_ALL_TESTS();
}
