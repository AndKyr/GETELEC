#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cmath>
#include <vector>
#include <limits>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <random>
#include <cassert>
#include <algorithm>

namespace getelec{

using namespace std;

/**
 * @struct PhysicalConstants
 * @brief Defines physical constants used throughout the project.
 * 
 * These constants include fundamental physical values such as Planck's constant,
 * the Boltzmann constant, and others used in various calculations.
 */
static constexpr struct PhysicalConstants {
    double hbarSqrOver2m = 3.80998212e-2; /**< Reduced Planck's constant squared over 2m (eV·nm²). */
    double kConstant = 1. / hbarSqrOver2m; /**< Constant connecting wavevector and energy */
    double sqrt2mOverHbar = sqrt(kConstant); /**<  sqrt(2m)/hbar. */
    double imageChargeConstant = 0.359991137; /**< Constant for image charge effects (eV·nm). */
    double BoltzmannConstant = 8.617333262e-5; /**< Boltzmann constant in eV/K. */
    double SommerfeldConstant = 1.618311e-4; /**< Sommerfeld constant for current density [A / nm^2 / (eV)^2 ]. */
    double exponentLimit = -0.5 * log(numeric_limits<double>::epsilon()); /**< Exponent limit to prevent underflow. */
    double electronCharge = 1.602176634e-19; /**< Electron charge in Cb. */
    
} CONSTANTS;

/**
 * @struct Utilities
 * @brief Provides a collection of static utility functions for numerical and physical calculations.
 * 
 * This struct includes functions for generating linearly spaced values, 
 * as well as computing the Fermi-Dirac distribution and its logarithmic counterpart.
 */
struct Utilities {
    /**
     * @brief Generates a vector of linearly spaced values.
     * 
     * @param start The starting value of the range.
     * @param end The ending value of the range.
     * @param n The number of points to generate.
     * @return A vector containing `n` linearly spaced values between `start` and `end`.
     * 
     * @note If `n == 1`, the result will contain only the `start` value.
     */
    static vector<double> linspace(double start, double end, int n);

    /**
     * @brief Computes the Fermi-Dirac distribution for a given energy.
     * 
     * @param energy The energy level in eV.
     * @param kT The thermal energy (kT) in eV.
     * @return The value of the Fermi-Dirac function at the specified energy level.
     * 
     * @details This function accounts for extreme energy cases to prevent underflow:
     * - For very high energy values, it approximates the exponential behavior.
     * - For very low energy values, it approximates to 1 minus the exponential.
     */
    static double fermiDiracFunction(double energy, double kT);

    /**
     * @brief Computes the logarithmic Fermi-Dirac distribution.
     * 
     * @param energy The energy level in eV.
     * @param kT The thermal energy (kT) in eV.
     * @return The logarithmic value of the Fermi-Dirac function at the specified energy level.
     * 
     * @details This function handles extreme cases:
     * - For very high energy values, it directly computes the exponential.
     * - For very low energy values, it calculates `-energy / kT` with adjustments.
     * - For intermediate values, it computes the logarithmic representation of the distribution.
     */
    static double logFermiDiracFunction(double energy, double kT);

    /**
     * @brief Generates a random double uniformly distributed in the range [a, b].
     * 
     * @param a The lower bound of the range.
     * @param b The upper bound of the range.
     * @return A random double value in the range [a, b].
     * 
     * @details This function uses the Mersenne Twister engine for random number generation.
     */
    static double getUniformRandomDouble(double a, double b, mt19937& generator) {
        std::uniform_real_distribution<> dis(a, b); // Distribution in range [a, b]
        return dis(generator);
    }

    /**
     * @brief Generates a vector of uniformly distributed random numbers
     * @param lowerLimit The bottom limit of the distribution
     * @param upperLimit The upper limit of the distribution
     * @param size The size of the generated vector
     * @param generator The RNG object to be used.
     */
    static vector<double> getUniformRandomDoubles(double lowerLimit, double upperLimit, int size, mt19937& generator){
        vector<double> result(size);
        std::uniform_real_distribution<> dis(lowerLimit, upperLimit); // Distribution in range [a, b]
        for (int i = 0; i < size; ++i) 
            result[i] = dis(generator);
        
        return result;
    }

    /**
     * @brief Calcualtes and returns the vector of indices that sort the input vector x
     * @param x The input vector
     * @return The sorting indices vector
     */
    static vector<size_t> sortingIndices(const vector<double>& x){
        vector<size_t> sortIndices(x.size()); // Fill indices with 0, 1, ..., n-1
        iota(sortIndices.begin(), sortIndices.end(), 0);
        auto comparisonLambda = [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; };
        sort(sortIndices.begin(), sortIndices.end(), comparisonLambda);
        return sortIndices;
    }

    /**
     * @brief Sort two vectors of data based on the first.
     * @param x The first vector (absiccae) to be sorted
     * @param y The second vector to be sorted based on x
     */
    static void sortCoordinates(vector<double>& x, vector<double>& y){
        assert(x.size() == y.size() && "x and y must have the same size");
        
        vector<size_t> sortIndices = sortingIndices(x); //get the indices that sort x

        //copy the vectors to avoid overwriting
        auto xCopy = x;
        auto yCopy = y;

        //transform x and y by replacing with saved copies of the indices.
        transform(sortIndices.begin(), sortIndices.end(), x.begin(), [&xCopy](size_t sortedIndex) { return xCopy[sortedIndex]; });
        transform(sortIndices.begin(), sortIndices.end(), y.begin(), [&yCopy](size_t sortedIndex) { return yCopy[sortedIndex]; });
    }
};


/**
 * @class FunctionInterpolator
 * @brief Provides functionality for interpolating a given function.
 * 
 * This class uses GSL's interpolation framework to evaluate functions at arbitrary points,
 * refine sampling to meet tolerance, and output spline nodes.
 */
class FunctionInterpolator {
protected:
    double absoluteTolerance = 1.e-12; /**< Absolute tolerance for the interpolation process. */
    double relativeTolerance = 1.e-4; /**< Relative tolerance for the interpolation process. */
    gsl_interp_accel* accelerator = NULL; /**< GSL interpolation accelerator for improved performance. */
    gsl_spline* spline = NULL; /**< GSL spline object for interpolation. */

    /**
     * @struct SplineElement
     * @brief Represents an element in the spline sampling list.
     */
    struct SplineElement {
        double x; /**< The x-coordinate of the sample. */
        double y; /**< The y-coordinate of the sample. */
        bool bisect; /**< Flag indicating if this element requires refinement. */

        /**
         * @brief Constructs a SplineElement with given x, y, and bisect values.
         * @param xx The x-coordinate.
         * @param yy The y-coordinate.
         * @param b Whether the element requires refinement.
         */
        SplineElement(double xx, double yy, bool b) : x(xx), y(yy), bisect(b) {}
    };

    list<SplineElement> samplingList; /**< List of sampling points used for spline interpolation. */

    /**
     * @brief Calculates the y-value for a given x-coordinate.
     * @param x The x-coordinate.
     * @return The y-value calculated for the given x.
     * 
     * @note The default implementation is `exp(x)`. Override for custom behavior.
     */
    virtual double calculateYforX(double x) { return exp(x); }

    /**
     * @brief Calculates the error between a calculated y-value and the spline estimate.
     * @param x The x-coordinate.
     * @param yCalculated The calculated y-value.
     * @return The error between the calculated and estimated y-values.
     */
    virtual double calculateError(double x, double yCalculated){
        return abs(yCalculated - gsl_spline_eval(spline, x, accelerator));
    }

    /**
     * @brief Calculates the tolerance for a given x and y-value.
     * @param x The x-coordinate.
     * @param yValue The y-value.
     * @return The calculated tolerance.
     */
    virtual double calculateTolerance(double x, double yValue){ return absoluteTolerance + relativeTolerance * abs(yValue) + x * 0.;}

    int updateSpline(); /**< Updates the spline object based on the current sampling list. */
    int refineSampling(); /**< Refines the sampling list to meet the specified tolerance. */

public:
    /**
     * @brief Constructs a FunctionInterpolator with specified tolerances.
     * @param aTol Absolute tolerance for interpolation.
     * @param rTol Relative tolerance for interpolation.
     */
    FunctionInterpolator(double aTol = 1.e-12, double rTol = 1.e-5)
        : absoluteTolerance(aTol), relativeTolerance(rTol) {}

    /**
     * @brief Destructor that frees GSL resources.
     */
    ~FunctionInterpolator() {
        gsl_spline_free(spline);
        gsl_interp_accel_free(accelerator);
    }

    /**
     * @brief Evaluates the spline at a given x-coordinate.
     * @param x The x-coordinate.
     * @return The interpolated y-value.
     */
    virtual double evaluate(double x) {
        return gsl_spline_eval(spline, x, accelerator);
    }

    /**
     * @brief Initializes the spline with initial sampling points.
     * @param xInit The initial x-coordinate.
     * @param xFinal The final x-coordinate.
     * @param numberOfInitialElements The number of initial sampling points.
     */
    void initialize(double xInit, double xFinal, int numberOfInitialElements);

    /**
     * @brief Refines the sampling list to meet the specified tolerance.
     * @param maxRefiningSteps The maximum number of refining steps allowed.
     * @return The number of refining steps actually taken
     */
    int refineToTolerance(int maxRefiningSteps = 10);

    /**
     * @brief Writes the spline nodes to a file.
     * @param filename The name of the output file.
     */
    void writeSplineNodes(string filename = "SplineNodes.dat");

    /**
     * @brief Getter for the absolute tolerance.
     */
    double getAbsoluteTolerance() const { return absoluteTolerance; }
    
    /**
     * @brief Getter for the relative tolerance.
     */
    double getRelativeTolerance() const { return relativeTolerance; }

    /**
     * @brief Gets the size of the spline sampling list.
     * @return The number of elements in the sampling list.
     */
    size_t size() const { return samplingList.size(); }
};

} // namespace getelec

#endif // UTILITIES_H_
