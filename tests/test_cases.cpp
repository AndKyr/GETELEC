#include <gtest/gtest.h>
#include "Utilities.h"
#include "ODESolver.h"
#include "TransmissionSolver.h"
#include "TransmissionInterpolator.h"
#include "TunnelingFunction.h"

// Test for Utilities::linspace
TEST(UtilitiesTest, Linspace) {
    auto result = Utilities::linspace(0.0, 10.0, 5);
    std::vector<double> expected = {0.0, 2.5, 5.0, 7.5, 10.0};
    EXPECT_EQ(result.size(), expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
        EXPECT_NEAR(result[i], expected[i], 1e-6);
    }
}

// Test for Utilities::fermiDiracFunction
TEST(UtilitiesTest, FermiDiracFunction) {
    double energy = 0.5;
    double kT = 0.025;
    double result = Utilities::fermiDiracFunction(energy, kT);
    EXPECT_GT(result, 0.0);
    EXPECT_LT(result, 1.0);
}

// Test for FunctionInterpolator
TEST(FunctionInterpolatorTest, Evaluate) {
    FunctionInterpolator interpolator;
    interpolator.initialize(0.0, 10.0, 5);
    double value = interpolator.evaluate(5.0);
    EXPECT_NEAR(value, std::exp(5.0), 1e-6);
}

// Test for TransmissionSolver::setEnergy
TEST(TransmissionSolverTest, SetEnergy) {
    TunnelingFunction tunnelingFunc;
    TransmissionSolver solver(&tunnelingFunc);
    solver.setEnergy(1.0);
    EXPECT_NEAR(tunnelingFunc.getEnergy(), 1.0, 1e-6);
}

// Test for TransmissionInterpolator::calculateYforX
TEST(TransmissionInterpolatorTest, CalculateYforX) {
    TunnelingFunction tunnelingFunc;
    TransmissionSolver solver(&tunnelingFunc);
    TransmissionInterpolator interpolator(solver);
    double value = interpolator.calculateYforX(0.5);
    EXPECT_LT(value, 0.0); // Transmission coefficient is expected to be less than 1
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
