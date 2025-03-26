#include <gtest/gtest.h>
#include "Utilities.h"
#include "ODESolver.h"
#include "TransmissionSolver.h"
#include "TransmissionInterpolator.h"
#include "TunnelingFunction.h"
#include "BandEmitter.h"
#include "Getelec.h"
#include "BSpline.h"

#include <random>
#include <iostream>

namespace getelec{
// Test for TransmissionSolver:: check that the transmission coefficient for the default barrier is calculated correctly
TEST(TransmissionSolverTest, DefaultValueTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(8.0);
    double transmission = solver.calculateTransmissionProbability(-4.5);
    EXPECT_NEAR(transmission, 0.00066697796753639933, 1.e-10);
}


TEST(TransmissionSolverTest, DerivativeTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    solver.setXlimits(10.0);

    double E = -4.5;
    double dE = 1.e-2;
    double deltaKappaSquared = CONSTANTS.kConstant * dE;

    solver.setEnergyAndInitialValues(-4.5);
    solver.solve(true);
    vector<double> y = solver.getSolution();
    solver.writeSolution("odeSolution.dat");

    solver.setEnergyAndInitialValues(-4.5 + dE);
    solver.solve(true);
    vector<double> y2 = solver.getSolution();
    solver.writeSolution("odeSolution2.dat");

    for (int i = 0; i < 3; i++){
        double derivative = y[i+3];
        double derivativeApprox = (y2[i] - y[i]) / deltaKappaSquared;
        EXPECT_NEAR(derivative, derivativeApprox, 1.e-2);
    }
}


// TEST(BSplineTest, ValueTest){
//     vector<double> x = Utilities::linspace(0., 2*M_PI, 6);
//     vector<vector<double>> valuesAndDerivatives;
    
//     int nDerivs = 1;
//     for (int i = 0; i <= nDerivs; i++){
//         auto val = vector<double>(6, 0.);
//         for (int j = 0; j < 6; j++){
//             val[j] = sin(x[j] + i * 0.5 * M_PI);
//         }
//         valuesAndDerivatives.push_back(val);
//     }
//     auto spline = GeneralizedHermiteSpline(x, valuesAndDerivatives);

//     auto cubicSpline = CubicHermiteSpline(x, valuesAndDerivatives[0], valuesAndDerivatives[1]);

//     auto testPoints = Utilities::linspace(0., 2*M_PI, 64);
//     for (auto x : testPoints){
//         double calculatedValue = sin(x);
//         double interpolatedValue = spline.evaluate(x);
//         double cubicSplineInterpolatedValue = cubicSpline.evaluate(x);

//         EXPECT_NEAR(calculatedValue, interpolatedValue, 1.e-3);
//         EXPECT_NEAR(calculatedValue, cubicSplineInterpolatedValue, 1.e-3);
//     }

// }

TEST(TransmissionSolverTest, SplineTest){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    solver.setXlimits(10.0);

    solver.setSolutionSplines(-6., 0., 8);
    solver.writeSplineSolution("splineSolution.dat");
}


// Test for TransmissionInterpolator:: check that the inteprolated and calculated values are close
TEST(TransmissionInterpolatorTest, evaluationTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(12.);

    const double bandDepth = 7.5;
    const double workFunction = 4.5;

    TransmissionInterpolator interpolator(solver, workFunction, 0.025, bandDepth);
    interpolator.initialize(-bandDepth+0.1, 1.0, 4);
    interpolator.refineToTolerance();


    mt19937 generator(1987);

    for (int i = 0; i < 64; i++){
        double testEnergy = Utilities::getUniformRandomDouble(-bandDepth+0.1, 1., generator);
        double waveVector = sqrt(testEnergy + bandDepth) * CONSTANTS.sqrt2mOverHbar;
        double calculatedValue = solver.calculateTransmissionProbability(-4.5 + testEnergy, waveVector);
        double interpolatedValue = interpolator.evaluate(testEnergy);

        double error = interpolator.calculateError(testEnergy, log(calculatedValue));
        double tolerance = interpolator.calculateTolerance(testEnergy, log(calculatedValue));
        EXPECT_LE(error, 10 * tolerance);
        
    }
}

TEST(BandEmitterTest, DefaultValueTest){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    BandEmitter emitter(solver);
    emitter.calculateCurrentDensityAndSpectra();
    double currentDensity = emitter.getCurrentDensity();
    double nottinghamHeat = emitter.getNottinghamHeat();
    EXPECT_NEAR(4.2191382392219547e-09, currentDensity, emitter.getToleranceForValue(currentDensity));
    EXPECT_NEAR(-8.7275419520125675e-10, nottinghamHeat, emitter.getToleranceForValue(nottinghamHeat));
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
    ifstream file("tempConfig.txt");
    string line;
    while (getline(file, line)) cout << line << endl;
    file.close();
    remove("tempConfig.txt");
}


TEST(GetelecObjectTest, RunParalleltest){
    Getelec getelec;
    auto fields = Utilities::linspace(2., 10., 64);
    getelec.setField(fields);
    EXPECT_NO_THROW(getelec.run());
}

TEST(GeneralXCFunctionTest, ValueTest){
    ModifiedSNBarrierWithDftXC barrier;
    barrier.setRadius(1.e3);
    auto xValues = Utilities::linspace(-0.1484608302148185, 1., 8);
    vector<double> expectedBarrierValues = {-27.33779495 ,  -4.500531602652986,  -2.555235757718333,
        -2.755385558646906,  -3.246461076744525,  -3.892712462231659,
        -4.606645547865526,  -5.354816780773875};
    vector<double> barrierValues;
    for (auto x : xValues){
        barrierValues.push_back(barrier.potentialFunction(x));
    }
    for (size_t i = 0; i < xValues.size(); i++){
        EXPECT_NEAR(barrierValues[i], expectedBarrierValues[i], 1.e-5);
    }
}

TEST(GeneralXCFunctionTest, DerivativeTest){
    ModifiedSNBarrierWithDftXC barrier;
    
    auto xValues = Utilities::linspace(-0.1484608302148185, 1., 64);
    double dx = 1.e-5;
    for (auto&x : xValues){
        double derivative = barrier.potentialFunctionDerivative(x);
        double derivativeApprox = (barrier.potentialFunction(x + dx) - barrier.potentialFunction(x - dx)) / (2*dx);
        EXPECT_NEAR(derivative, derivativeApprox, 1.e-5);
    }

}


} // namespace getelec

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(catch_exceptions) = false;

    return RUN_ALL_TESTS();
}