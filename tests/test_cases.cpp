#include <gtest/gtest.h>
#include "Utilities.h"
#include "ODESolver.h"
#include "TransmissionSolver.h"
#include "TransmissionInterpolator.h"
#include "TunnelingFunction.h"
#include "BandEmitter.h"
#include "Getelec.h"
#include "BSpline.h"
#include "TransmissionSplines.h"

#include <random>
#include <iostream>

namespace getelec{
/**
 * @brief Test for TransmissionSolver:: check that the transmission coefficient for the default barrier is the one expected
 */
TEST(TransmissionSolverTest, DefaultValueTest) {
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier);
    solver.setXlimits(8.0);
    double transmission = solver.calculateTransmissionProbability(-4.5);
    EXPECT_NEAR(transmission, 0.00066697796753639933, 1.e-9);
    EXPECT_NEAR(transmission, 0.00066697796753639933, 1.e-9);
}

/**
 * @brief Test for TransmissionSolver:: check that the derivative of the transmission solution is correct
 * @details This test checks the accuracy of the derivative of the transmission solution calculated by the TransmissionSolver class.
 *        It compares the derivative of the solution with an approximate value calculated using finite differences.
 */
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

    solver.setEnergyAndInitialValues(-4.5 + dE);
    solver.solve(true);
    vector<double> y2 = solver.getSolution();

    for (int i = 0; i < 3; i++){
        double derivative = y[i+3];
        double derivativeApprox = (y2[i] - y[i]) / deltaKappaSquared;
        EXPECT_NEAR(derivative, derivativeApprox, 1.e-2);
    }
}

/**
 * @brief Test for TransmissionSolver:: check that the top barrier is found correctly
 * @details This test checks the accuracy of the top barrier found by the TransmissionSolver class.
 *        It compares the calculated top barrier with the value found by densely sampling the barrier
 */
TEST(TransmissionSolverTest, TopBarrierFinderTest){
    ModifiedSNBarrier barrier;
    barrier.setBarrierTopFinder(true);

    TransmissionSolver solver(&barrier);
    solver.setXlimits(10.0);
    solver.setEnergyAndInitialValues(-4.5);
    solver.solveNoSave();

    double topBarrier = -numeric_limits<double>::infinity();
    vector<double> xPoints = Utilities::linspace(solver.getXFinal(), solver.getXInitial(), 256);
    for (auto x : xPoints){
        if (barrier.potentialFunction(x) > topBarrier)
            topBarrier = barrier.potentialFunction(x);
    }

    EXPECT_NEAR(topBarrier, barrier.getBarrierTop(), 1.e-2);
}

/**
 * @brief Test for TransmissionSpline:: check that the interpolated and calculated values are close
 * @details This test checks the accuracy of the TransmissionSpline class by comparing the interpolated values
 *          with the calculated values from the TransmissionSolver class. It uses a modified SN barrier and
 *          a TransmissionSolver to calculate the transmission probability and emission current estimates.
 *          The test samples the transmission coefficient at various energy levels and checks if the
 *          interpolated values are within a specified tolerance of the calculated values.
 * @note The test uses a relative tolerance of 1.e-2 and an absolute tolerance of 1.e-13 for the comparison.
 * @note The test also writes the spline solution and nodes to files for further analysis.
 */
TEST(TransmissionSplineTest, SplineTest){
    double atol = 1.e-13;
    double rtol = 1.e-2;
    double workFunction = 4.5;
    double kT = 0.025;
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    TransmissionSpline interpolator(solver, workFunction, kT, atol, rtol);
    solver.setXlimits(12.0);

    interpolator.smartInitialSampling();
    interpolator.refineSamplingToTolerance();

    interpolator.writeSplineSolution("splineSolution.dat", 256);
    interpolator.writeSplineNodes("splineNodes.dat");

    auto testEnergies = Utilities::linspace(interpolator.getMinimumSampleEnergy(), interpolator.getMaximumSampleEnergy(), 64);
    

    double k = sqrt(7.) * CONSTANTS.sqrt2mOverHbar;
    for (auto energy : testEnergies){
        solver.calculateTransmissionProbability(energy, k);
        double solverCurrentEstimate = solver.getEmissionEstimate(k, workFunction, kT);
        double interpolatedCurrentEstimate = interpolator.normalEnergyDistributionEstimate(energy, k);
        double tolerance = atol + interpolator.getmaximumCurrentEstimate() * rtol;
        EXPECT_NEAR(solverCurrentEstimate, interpolatedCurrentEstimate, tolerance);
    }
}

/**
 * @brief Test for TransmissionSpline:: check that the maximum emission current estimate is found correctly
 * @details This test checks the accuracy of the maximum emission current estimate found by the TransmissionSpline class.
 */
TEST(TransmissionSplineTest, maxEmissionFindTest){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    TransmissionSpline interpolator(solver);
    solver.setXlimits(12.0);

    double tolerance = 0.01;
    int maxIterations = 15;
    interpolator.smartInitialSampling();
    int iterations;
    EXPECT_NO_THROW(iterations = interpolator.findMaximumCurrentEstimate(maxIterations, tolerance));

    auto energies = Utilities::linspace(interpolator.getMinimumSampleEnergy(), interpolator.getMaximumSampleEnergy(), 256);
    
    double maxEmission = 0.;

    for (auto energy : energies){
        double emissionEstimate = interpolator.normalEnergyDistributionEstimate(energy);
        if (emissionEstimate > maxEmission) 
            maxEmission = emissionEstimate;
    }
    EXPECT_NEAR(maxEmission, interpolator.getmaximumCurrentEstimate(), 2*tolerance);

    cout << "maximum found within " << iterations << " iterations" << endl;
}

/**
 * @brief Test for BandEmitter:: check that the current density and Nottingham heat are calculated correctly
 * @details This test checks the accuracy of the current density and Nottingham heat calculated by the BandEmitter class.
 *         It uses a modified SN barrier and a TransmissionSolver to calculate the current density and Nottingham heat.
 *        The test compares the calculated values with known expected values to ensure correctness.
 */
TEST(BandEmitterTest, DefaultValueTest){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    BandEmitter emitter(solver);
    emitter.integrateTotalEnergyDistributionODEAndSaveSpectra();
    double currentDensity = emitter.getCurrentDensityODE();
    double nottinghamHeat = emitter.getNottinghamHeatODE();
    EXPECT_NEAR(4.1646907052681002e-09, currentDensity, emitter.getToleranceForValue(currentDensity));
    EXPECT_NEAR(-8.660864615724323e-10, nottinghamHeat, emitter.getToleranceForValue(nottinghamHeat));
}

/**
 * @brief Calculates the total energy distribution with two methods and checks that they match
 * @details This test calculates the total energy distribution by two independent methods: ordinary differential equation (fast but specific for effectiveMass = 1)
 *  and integration of the double integral over parallel energies (slower but more general). Then it demands that they match to the tolerance of the band emitter.
 *  The test is run for 64 different random parameters of the barrier, and the results are compared.
 */
TEST(BandEmitterTest, totalEnergyDistributionMethodComparison){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    BandEmitter emitter(solver);
    mt19937 generator(1987);
    emitter.setGenerator(&generator);
    barrier.setGenerator(&generator);

    for (int i = 0; i < 64; i++){
        barrier.setRandomParameters();
        emitter.setParameters();
        emitter.integrateTotalEnergyDistributionODEAndSaveSpectra();
        emitter.writePlottingData("bandEmitterPlotting.dat");
        auto [totalEnergies, totalEnergyDistributions, dummyDerivatives] = emitter.getSpectra();

        for (size_t j = 0; j < totalEnergies.size(); j++){
            double TEDFromIntegral = emitter.totalEnergyDistributionIntegrateParallel(totalEnergies[j]);
            EXPECT_NEAR(totalEnergyDistributions[j], TEDFromIntegral , 5 * emitter.getToleranceForValue(totalEnergyDistributions[j]));
        }
    }
}


/**
 * @brief Test for BandEmitter:: check that the current density and Nottingham heat are calculated correctly
 * @details This test checks the accuracy of the current density and Nottingham heat calculated by the BandEmitter class.
 *         It uses a modified SN barrier and a TransmissionSolver to calculate the current density and Nottingham heat.
 *        The test compares the calculated values with two independent methods (TED vs NED integration) and confirms that they match.
 */
TEST(BandEmitterTest, CurrentDensityMethodComparison){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    BandEmitter emitter(solver);
    mt19937 generator(1987);
    emitter.setGenerator(&generator);
    barrier.setGenerator(&generator);

    for (int i = 0; i < 64; i++){
        barrier.setRandomParameters();
        emitter.setParameters();
        emitter.integrateTotalEnergyDistributionODE();
        double currentDensity = emitter.getCurrentDensityODE();
        double nottingham = emitter.getNottinghamHeatODE();
        double nottingham2 = emitter.nottinghamIntegrateTotalPrallel();

        double currentDensity2 = emitter.currentDensityIntegrateNormal();
        double currentDensity3 = emitter.currentDensityIntegrateTotalParallel();
        double currentDensity4 = emitter.currentDensityIntegrateParallelTotal();
        double currentDensity5 = emitter.currentDensityIntegrateNormalParallel();
        EXPECT_NEAR(currentDensity, currentDensity2, 10*emitter.getToleranceForValue(currentDensity));
        EXPECT_NEAR(currentDensity, currentDensity3, 10*emitter.getToleranceForValue(currentDensity));
        EXPECT_NEAR(currentDensity, currentDensity4, 10*emitter.getToleranceForValue(currentDensity));
        EXPECT_NEAR(currentDensity, currentDensity5, 10*emitter.getToleranceForValue(currentDensity));
        EXPECT_NEAR(nottingham, nottingham2, 100*emitter.getToleranceForValue(nottingham));
    }
}

/**
 * @brief Test for BandEmitter:: check that the normal energy distribution is calculated correctly
 * @details This test checks the accuracy of the normal energy distributtion calculated by the BandEmitter class.
 *         It uses a modified SN barrier and a TransmissionSolver to calculate the current density and Nottingham heat.
 *        The test compares the calculated values with two independent methods (single integral vs double integral) and confirms that they match.
 */
TEST(BandEmitterTest, normalEnergyDistributionMethodComparison){
    ModifiedSNBarrier barrier;
    TransmissionSolver solver(&barrier, Config().transmissionSolverParams, 10., 1);
    BandEmitter emitter(solver);
    mt19937 generator(1987);
    emitter.setGenerator(&generator);
    barrier.setGenerator(&generator);

    for (int i = 0; i < 64; i++){
        barrier.setRandomParameters();
        emitter.setParameters();
        double currentDensity = emitter.currentDensityIntegrateNormalParallel(true);
        auto [normalEnergies, normalEnergyDistribution] = emitter.getNormalEnergyDistribution();

        for (size_t j = 0; j < normalEnergies.size(); j++){
            double NEDSimple = emitter.normalEnergyDistributionForEnergy(normalEnergies[j]) * emitter.getkT() * CONSTANTS.SommerfeldConstant;
            EXPECT_NEAR(normalEnergyDistribution[j], NEDSimple , currentDensity / (emitter.getXFinal() - emitter.getXInitial()));
        }
    }
}

/**
 * @brief Test for reading the configuration file
 * @details This test checks the ability to read a configuration file and verify that the parameters are set correctly.
 *         It creates a temporary configuration file, modifies some parameters, and then reads the file to check if the values match.
 *         The test also prints the contents of the configuration file to the console for debugging purposes.
 *         Finally, it removes the temporary file after the test is complete.
 */
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


/**
 * @brief Test for Getelec:: check that the Getelec class can be constructed and run without errors
 */
TEST(GetelecTest, runRandomCasesTest){
    tbb::global_control tbbGlobalControl(tbb::global_control::max_allowed_parallelism, 1);

    mt19937 generator(1987);

    Getelec getelec = Getelec("GetelecConfig.txt", "modifiedSN", &generator);
    // getelec.setGenerator(&generator);
    getelec.setRandomParameters(64);
    EXPECT_NO_THROW(getelec.run(CalculationFlags::All));
}


/**
 * @brief Test for Getelec:: check that the Getelec class can be constructed and run without errors
 */
TEST(GetelecObjectTest, RunParalleltest){
    Getelec getelec;
    auto fields = Utilities::linspace(2., 10., 64);
    getelec.setField(fields);
    EXPECT_NO_THROW(getelec.run());
}



/**
 * @brief Test for ModifiedSNBarrierWithDftXC:: check that the potential function and its derivative are calculated correctly
 */
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

/**
 * @brief Test for ModifiedSNBarrierWithDftXC:: check that the derivative of the potential function is calculated correctly
 * @details This test checks the accuracy of the derivative of the potential function calculated by the ModifiedSNBarrierWithDftXC class.
 *         It compares the calculated derivative with an approximate value calculated using finite differences.
 */
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