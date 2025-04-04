#ifndef BSPLINE_H_
#define BSPLINE_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf.h>

#include <typeinfo>
#include <vector>
#include <stdexcept>
#include <cassert>

using namespace std;

namespace getelec{

/**
 * @class BSpline
 * @brief Class to handle B-spline interpolation (wrapping to native GSL code)
 * @note This class should not be used directly, but rather through derived classes
 * @note The class (and its derived ones) uses raw C-style arrays for the input data, so it is not safe to use with vectors or other containers, nor with dynamically allocated arrays, as the class does not take ownership of the data
 * @note The class is not thread-safe
 * @note The class is not copyable
 * @note The class is not movable
 */
class BSpline{
public:

    /**
     * @brief Construct a new BSpline object
     */
    BSpline() = default;

    /**
     * @brief Set the positions for the spline (control points)
     * @param x Vector of positions
     */
    void setPositions(const vector<double>& x);

    /**
     * @brief Set the positions for the spline (control points) and the values of the function at the points
     * @param x Vector of positions
     * @param y Vector of values
     */
    void setPositionsAndValues(const vector<double>& x, const vector<double>& y);

    /**
     * @brief Set the positions for the spline (control points) and the values of the function at the points
     * @param x Vector of positions
     * @param y Vector of values
     */
    void setPositionsAndValues(gsl_vector* x, gsl_vector* y);

    /**
     * @brief Delete the object and deallocate gsl objects
     */
    ~BSpline() {
        if (workSpace) {
            gsl_bspline_free(workSpace);
            workSpace = nullptr;
        }
        if (coefficients) {
            gsl_vector_free(coefficients);
            coefficients = nullptr;
        }
        
        if (positions) {
            gsl_vector_free(positions);
            positions = nullptr;
        }

        if (values) {
            gsl_vector_free(values);
            values = nullptr;
        }

        for (auto& coeffSet : coefficientSets){
            gsl_vector_free(coeffSet);
        }
    }

    /**
     * @brief Evaluate the spline at a given point
     * @param x Point at which to evaluate the spline
     * @return Value of the spline at the given point
     */
    double evaluate(double x) const;

    /**
     * @brief Evaluate multiple functions at a given point
     * @param x Point at which to evaluate the functions
     * @return Vector of values of the functions at the given point
     */
    vector<double> evaluateMultiple(double x) const;


    double evaluateDerivative(double x, size_t nDerivative = 1) const;

    /**
     * @brief Evaluate the derivatives of multiple splines (on the same workspace and absiccae) at a given point
     * @param x Point at which to evaluate the derivative
     * @param nDerivative Order of the derivative to evaluate
     * @return Value of the derivative at the given point
     */
    vector<double> evaluateDerivativeMultiple(double x, size_t nDerivative = 1) const ;

    /**
     * @brief Evaluate the integral of the spline between two points
     * @param a Lower limit of the integral
     * @param b Upper limit of the integral
     * @return Value of the integral between the two points
     */
    double evaluateIntegral(double a, double b) const;

    /**
     * @brief Evaluate the integral of multiple splines between two points
     * @param a Lower limit of the integral
     * @param b Upper limit of the integral
     * @return Vector of values of the integrals between the two points
     * @note The function assumes that the splines are on the same workspace and have the same abscissa
     */
    vector<double> evaluateIntegralMultiple(double a, double b) const;

protected:
    gsl_bspline_workspace* workSpace = nullptr;  // Workspace for spline evaluation
    gsl_vector* coefficients = nullptr;
    vector<gsl_vector*> coefficientSets;  
    gsl_vector* positions = nullptr;
    gsl_vector* values = nullptr;
};

/**
 * @class CubicHermiteSpline
 * @brief Class to handle cubic Hermite spline interpolation
 * @note The class is not thread-safe
 * @note The class is not copyable
 * @note The class is not movable
 */
class CubicHermiteSpline : public BSpline{
public:
    CubicHermiteSpline() = default;

    /**
     * @brief Initialize a new CubicHermiteSpline object from positions, values and derivatives
     * @param x Vector of positions
     * @param y Vector of values
     * @param dy_dx Vector of derivatives
     */
    void initialize(const vector<double>& x, const vector<double>& y, const vector<double>& dy_dx);

    /**
     * @brief Initialize a new CubicHermiteSpline object from positions, values and derivatives
     * @param x Vector of positions
     * @param y Vector of values
     * @param dy_dx Vector of derivatives
     * @note The input vectors must have the same size
     * @note The class takes no ownership of the input vectors. They need to be deallocated by the user 
     */
    void initialize(gsl_vector* x, gsl_vector* y, gsl_vector* dy_dx);

    /**
     * @brief Initialize with multiple functions from positions, values and derivatives
     * @param x Vector of positions
     * @param ySets Vector of vectors of values
     * @param dy_dxSets Vector of vectors of derivatives 
     */
    void initializeMultiple(const vector<double>& x, const vector<vector<double>>& ySets, const vector<vector<double>>& dy_dxSets);

    /**
     * @brief Delete the object and deallocate gsl objects
     * @note The destructor is virtual to allow for proper cleanup of derived classes
     */
    ~CubicHermiteSpline() {
        if (derivs){
            gsl_vector_free(derivs);
            derivs = nullptr;
        }
    }

private:
    gsl_vector* derivs = nullptr;

    void reallocateWorkSpace();

    void reallocateCoefficients();

    void allocateCoefficients();

    void reallocate();
    
    void intitialize();
};

}


// class GeneralizedHermiteSpline : public BSpline{

// public:
//     GeneralizedHermiteSpline(vector<double>& x, vector<vector<double>>& valuesAndDerivatives_) : BSpline(x, valuesAndDerivatives_[0]) {
//         size_t nDeriv = valuesAndDerivatives_.size() - 1;
//         size_t nPoints = x.size();
//         valuesAndDerivatives = gsl_matrix_alloc(nPoints, nDeriv+1);
//         for (size_t i = 0; i < nPoints; i++) {
//             for (size_t j = 0; j < nDeriv; j++) {
//                 gsl_matrix_set(valuesAndDerivatives, i, j, valuesAndDerivatives_[j][i]);
//             }
//         }
//         initialize();
//     }

//     GeneralizedHermiteSpline(gsl_vector* x, gsl_matrix* valuesAndDerivatives_, size_t nDeriv) : 
//         BSpline(x, nullptr) 
//     {
//         valuesAndDerivatives = valuesAndDerivatives_;
//         initialize();
//     }

//     ~GeneralizedHermiteSpline() {
//         gsl_matrix_free(valuesAndDerivatives);
//     }

// private:
//     gsl_matrix* valuesAndDerivatives;

//     void initialize(){
//         size_t nDeriv = valuesAndDerivatives->size2 - 1;
//         size_t nPoints = valuesAndDerivatives->size1;
//         coefficients = gsl_vector_alloc((nDeriv+1) * nPoints);
//         workSpace = gsl_bspline_alloc_ncontrol(2*nDeriv + 2, (nDeriv + 1) * nPoints);
//         gsl_bspline_init_hermite(nDeriv, positions, workSpace);
//         gsl_bspline_interp_hermite_general();
//     }


//     /* --------------------------------------------------------------------------
//     Implements generalized Hermite interpolation with B-splines, based on
//     M. S. Mummy, "Hermite Interpolation with B-splines," CAGD 6 (1989), 177–179.

//     We assume:
//         - x:  n interpolation nodes (gsl_vector of length n)
//         - valuesAndDerivatives:  n x (m+1) matrix of derivative data:
//             valuesAndDerivatives(i,0) = f(x_i),
//             valuesAndDerivatives(i,1) = f'(x_i),
//             ...
//             valuesAndDerivatives(i,m) = f^(m)(x_i).
//         - m:  highest derivative order to match
//         - coefficients:  output B-spline coefficients (length (m+1)*n)
//         - w:  gsl_bspline_workspace with order = 2*m + 2,
//             knots set by gsl_bspline_init_hermite(m, x, w)

//     This code constructs Hermite B-spline coefficients via:
//         A_i^r = sum_{s=0..m} [ choose(2m+1 - s, m-s) *
//                                 ( (x_{i+1} - x_i)^s / s! ) *
//                                 valuesAndDerivatives(i, m-s) ],
//     storing them in coefficients at index i*(m+1) + r.

//     For i = n-1, we take x_{i+1} = x_i (end-knot repetition).
//     -------------------------------------------------------------------------- */
//     int gsl_bspline_interp_hermite_general(){
//         const size_t n = positions->size; //number of interpolation nodes
//         const size_t m = valuesAndDerivatives->size2; //number of values to be matched at each node


//         /* Build the coefficients for i=0..(n-1). At the first and last, x[-1] = x[0] and x[n] = x[n-1] */
//         for (size_t i = 0; i < n; i++) {
//             const double xi   = gsl_vector_get(positions, i); //get the x value at the i-th node
//             const double xip1 = (i < n-1) ? gsl_vector_get(positions, i + 1) : xi;  /* repeated end knot for i=n-1 */
            
//             const double xim1 = i > 0 ? gsl_vector_get(positions, i - 1) : xi;  /* repeated end knot for i=n-2 */
//             const double deltaPlus = xip1 - xi;
//             const double deltaMinus = xi - xim1;
        
//             /* For each r in 0..m, store A_i^r. */
//             for (int j = 0; j < m; j++) {
//                 double sum_r = 0.0;
            
//                 /* Sum over s=0..m in Mummy's eq. (3.1). */
//                 for (int s = 0; s < m; s++) {
//                     /* binomial coefficient: choose(2m+1 - s, m-s). */
//                     const double binomialCoefficient = gsl_sf_choose(m, s) / gsl_sf_choose(2*m + 1, s);
            
//                     /* factor = (delta^s / s!) */
//                     double ds_term = pow(deltaPlus, s) * pow(deltaMinus, m-s-1) / gsl_sf_fact(s);
            
//                     /* derivative data valuesAndDerivatives(i, m-s) = f^(m-s)(x_i). */
//                     const double f_deriv = gsl_matrix_get(valuesAndDerivatives, i, s);
            
//                     sum_r += binomialCoefficient * ds_term * f_deriv;
//                 }
            
//                 /* Store in coefficients at index [i*(m+1) + r]. */
//                 gsl_vector_set(coefficients, i * m + j, sum_r);
//             }
//         }
        
//         return GSL_SUCCESS;
//     }
// };

#endif //BSPLINE_H_