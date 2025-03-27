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

using namespace std;

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
    void setPositions(vector<double>& x){
        if (positions)
            gsl_vector_free(positions);
        positions = gsl_vector_alloc(x.size());
        for (size_t i = 0; i < x.size(); i++) 
            gsl_vector_set(positions, i, x[i]); 
    }

    /**
     * @brief Set the positions for the spline (control points) and the values of the function at the points
     * @param x Vector of positions
     * @param y Vector of values
     */
    void setPositionsAndValues(vector<double>& x, vector<double>& y){

        // Ensure the input vectors have the same size
        assert(x.size() == y.size() && "Input vectors must have the same size");

        //reallocate the positions and values
        if (positions)
            gsl_vector_free(positions);
        positions = gsl_vector_alloc(x.size());
        for (size_t i = 0; i < x.size(); i++) 
            gsl_vector_set(positions, i, x[i]); 

        if (values)
            gsl_vector_free(values);
        values = gsl_vector_alloc(y.size());
        for (size_t i = 0; i < y.size(); i++) 
            gsl_vector_set(values, i, y[i]); 

    }

    void setPositionsAndValues(gsl_vector* x, gsl_vector* y){
        // Ensure the input vectors have the same size
        assert(x->size == y->size && "Input vectors must have the same size");
        positions = x;
        values = y;
    }

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
    double evaluate(double x) const {
        assert(workSpace && coefficients && "The spline has not been initialized");
        double result;
        gsl_bspline_calc(x, coefficients, &result, workSpace);
        return result;
    }

    /**
     * @brief Evaluate multiple functions at a given point
     * @param x Point at which to evaluate the functions
     * @return Vector of values of the functions at the given point
     */
    vector<double> evaluateMultiple(double x) const{
        vector<double> result(coefficientSets.size());
        for (size_t i = 0; i < coefficientSets.size(); i++){
            assert(workSpace && coefficientSets[i] && "The spline has not been initialized");
            gsl_bspline_calc(x, coefficientSets[i], &result[i], workSpace);
        }
        return result;
    }


    double evaluateDerivative(double x, size_t nDerivative = 1) const {
        assert(workSpace && coefficients && "The spline has not been initialized");
        double result;
        gsl_bspline_calc_deriv(x, coefficients, nDerivative, &result, workSpace);
        return result;
    }

    vector<double> evaluateDerivativeMultiple(double x, size_t nDerivative = 1) const {
        vector<double> result(coefficientSets.size());
        for (size_t i = 0; i < coefficientSets.size(); i++){
            assert(workSpace && coefficientSets[i] && "The spline has not been initialized");
            gsl_bspline_calc_deriv(x, coefficientSets[i], nDerivative, &result[i], workSpace);
        }
        return result;
    }

    double evaluateIntegral(double a, double b){
        assert(workSpace != nullptr && coefficients != nullptr);

        double result;
        gsl_bspline_calc_integ(a, b, coefficients, &result, workSpace);
        return result;
    }

    vector<double> evaluateIntegralMultiple(double a, double b){
        vector<double> result(coefficientSets.size());
        for (size_t i = 0; i < coefficientSets.size(); i++){
            assert(workSpace && coefficientSets[i] && "The spline has not been initialized");
            gsl_bspline_calc_integ(a, b, coefficientSets[i], &result[i], workSpace);
        }
        return result;
    }


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
    void initialize(vector<double>& x, vector<double>& y, vector<double>& dy_dx) {
            // Ensure the input vectors have the same size
        setPositionsAndValues(x, y);
        
        if (derivs)
            gsl_vector_free(derivs);
        derivs = gsl_vector_alloc(dy_dx.size());
        for (size_t i = 0; i < dy_dx.size(); i++) 
            gsl_vector_set(derivs, i, dy_dx[i]); 
        
        reallocate();
        intitialize();
    }

    /**
     * @brief Initialize a new CubicHermiteSpline object from positions, values and derivatives
     * @param x Vector of positions
     * @param y Vector of values
     * @param dy_dx Vector of derivatives
     * @note The input vectors must have the same size
     * @note The class takes no ownership of the input vectors. They need to be deallocated by the user 
     */
    void initialize(gsl_vector* x, gsl_vector* y, gsl_vector* dy_dx) {
        // Ensure the input vectors have the same size
        if (x->size != y->size || x->size != dy_dx->size)
            throw invalid_argument("Input vectors must have the same size");
        
        positions = x;
        values = y;
        derivs = dy_dx;
        reallocate();
        intitialize();

        // Since the object has no ownership of the data, avoid freeing the input vectors
        positions = nullptr;
        values = nullptr;
        derivs = nullptr;
    }

    /**
     * @brief Initialize with multiple functions from positions, values and derivatives
     * @param x Vector of positions
     * @param ySets Vector of vectors of values
     * @param dy_dxSets Vector of vectors of derivatives 
     */
    void initializeMultiple(vector<double>& x, vector<vector<double>>& ySets, vector<vector<double>>& dy_dxSets){
        setPositions(x);
        assert(ySets.size() == dy_dxSets.size() && "The number of sets of values and derivatives must be the same");
        
        size_t nPoints = x.size();
        values = gsl_vector_alloc(nPoints);
        derivs = gsl_vector_alloc(nPoints);
        reallocate();

        for (int i = 0; i < ySets.size(); i++){
            assert(ySets[i].size() == nPoints && dy_dxSets[i].size() == nPoints && "The number of points in the sets must be the same");
            for (size_t j = 0; j < nPoints; j++){
                gsl_vector_set(values, j, ySets[i][j]);
                gsl_vector_set(derivs, j, dy_dxSets[i][j]);
            }
            allocateCoefficients();
            intitialize();
            coefficientSets.push_back(coefficients);
        }
        coefficients = nullptr;
    }



    ~CubicHermiteSpline() {
        if (derivs){
            gsl_vector_free(derivs);
            derivs = nullptr;
        }

    }

private:
    gsl_vector* derivs = nullptr;

    void reallocateWorkSpace(){
        size_t nControlPoints = 2 * positions->size;
        if (workSpace)
            gsl_bspline_free(workSpace);
        workSpace = gsl_bspline_alloc_ncontrol(4, nControlPoints);
    }

    void reallocateCoefficients(){
        size_t nControlPoints = 2 * positions->size;
        if (coefficients)
            gsl_vector_free(coefficients);
        coefficients = gsl_vector_alloc(nControlPoints);  
    }

    void allocateCoefficients(){
        size_t nControlPoints = 2 * positions->size;
        coefficients = gsl_vector_alloc(nControlPoints);
    }

    void reallocate(){
        reallocateWorkSpace();
        reallocateCoefficients();
    }
    
    void intitialize(){
        gsl_bspline_init_hermite(1, positions, workSpace);
        gsl_bspline_interp_chermite(positions, values, derivs, coefficients, workSpace);
    }
};


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