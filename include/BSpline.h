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

class BSpline{
public:
    BSpline() = default;
    void init(vector<double>& x, vector<double>& y){

        // Ensure the input vectors have the same size
        if (x.size() != y.size())
            throw invalid_argument("Input vectors must have the same size");

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

    void init(gsl_vector* x, gsl_vector* y){

        // Ensure the input vectors have the same size
        if (x->size != y->size)
            throw invalid_argument("Input vectors must have the same size");

        positions = x;
        values = y;
    }

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
    }

    // Function to evaluate the spline at a given point
    double evaluate(double x) const {
        assert(workSpace != nullptr && coefficients != nullptr);
        double result;
        gsl_bspline_calc(x, coefficients, &result, workSpace);
        return result;
    }

    // Function to evaluate the derivative of the spline at a given point
    double evaluateDerivative(double x, size_t nDerivative = 1) const {
        assert(workSpace != nullptr && coefficients != nullptr);

        double result;
        gsl_bspline_calc_deriv(x, coefficients, nDerivative, &result, workSpace);
        return result;
    }

    double evaluateIntegral(double a, double b){
        assert(workSpace != nullptr && coefficients != nullptr);

        double result;
        gsl_bspline_calc_integ(a, b, coefficients, &result, workSpace);
        return result;
    }

protected:
    gsl_bspline_workspace* workSpace = nullptr;  // Workspace for spline evaluation
    gsl_vector* coefficients = nullptr;
    gsl_vector* positions = nullptr;
    gsl_vector* values = nullptr;
};

class CubicHermiteSpline : public BSpline{
public:
    CubicHermiteSpline() = default;

    void init(vector<double>& x, vector<double>& y, vector<double>& dy_dx) {
            // Ensure the input vectors have the same size
        if (x.size() != y.size() || x.size() != dy_dx.size())
            throw invalid_argument("Input vectors must have the same size");
        
        BSpline::init(x, y);
        if (derivs)
            gsl_vector_free(derivs);
        derivs = gsl_vector_alloc(dy_dx.size());
        for (size_t i = 0; i < dy_dx.size(); i++) 
            gsl_vector_set(derivs, i, dy_dx[i]); 
        
        intitialize();
    }

    void init(gsl_vector* x, gsl_vector* y, gsl_vector* dy_dx) {
        // Ensure the input vectors have the same size
        if (x->size != y->size || x->size != dy_dx->size)
            throw invalid_argument("Input vectors must have the same size");
        
        positions = x;
        values = y;
        derivs = dy_dx;
        intitialize();
    }

    ~CubicHermiteSpline() {
        if (derivs){
            gsl_vector_free(derivs);
            derivs = nullptr;
        }

    }

private:
    gsl_vector* derivs = nullptr;
    
    void intitialize(){
        size_t nControlPoints = 2 * positions->size;
        if (workSpace)
            gsl_bspline_free(workSpace);
        if (coefficients)
            gsl_vector_free(coefficients);
        workSpace = gsl_bspline_alloc_ncontrol(4, nControlPoints);
        coefficients = gsl_vector_alloc(nControlPoints);
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