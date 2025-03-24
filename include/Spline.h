#ifndef SPLINE_H_
#define SPLINE_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_bspline.h>

#include <typeinfo>
#include <vector>
#include <stdexcept>

using namespace std;

class BSpline{
public:
    BSpline(vector<double>& x, vector<double>& y){

        // Ensure the input vectors have the same size
        if (x.size() != y.size() || x.size())
            throw invalid_argument("Input vectors must have the same size");

        positions = gsl_vector_alloc(x.size());
        for (size_t i = 0; i < x.size(); i++) 
            gsl_vector_set(positions, i, x[i]); 

        values = gsl_vector_alloc(y.size());
        for (size_t i = 0; i < y.size(); i++) 
            gsl_vector_set(values, i, y[i]); 

    }

    BSpline(gsl_vector* x, gsl_vector* y){

        // Ensure the input vectors have the same size
        if (x->size != y->size)
            throw invalid_argument("Input vectors must have the same size");

        positions = x;
        values = y;
    }

    ~BSpline() {
        gsl_bspline_free(workSpace);
        gsl_vector_free(coefficients);
        gsl_vector_free(positions);
    }

    // Function to evaluate the spline at a given point
    double evaluate(double x) const {
        double result;
        gsl_bspline_calc(x, coefficients, &result, workSpace);
        return result;
    }

    // Function to evaluate the derivative of the spline at a given point
    double evaluateDerivative(double x) const {
        double result;
        gsl_bspline_calc_deriv(x, coefficients, 1, &result, workSpace);
        return result;
    }

    double evaluateIntegral(double a, double b){
        double result;
        gsl_bspline_calc_integ(a, b, coefficients, &result, workSpace);
        return result;
    }

protected:
    gsl_bspline_workspace* workSpace;  // Workspace for spline evaluation
    gsl_vector* coefficients;
    gsl_vector* positions;
    gsl_vector* values;
};

class HermiteSpline : public BSpline{
public:
    HermiteSpline(vector<double>& x, vector<double>& y, vector<double>& dy_dx) : BSpline(x, y) {
            // Ensure the input vectors have the same size
        if (x.size() != y.size() || x.size() != dy_dx.size())
            throw invalid_argument("Input vectors must have the same size");

        gsl_vector* derivs = gsl_vector_alloc(y.size());
        for (size_t i = 0; i < dy_dx.size(); i++) gsl_vector_set(derivs, i, dy_dx[i]); 

        gsl_bspline_init_hermite(1, positions, workSpace);
        gsl_bspline_interp_chermite(positions, values, derivs, coefficients, workSpace);
    }

    HermiteSpline(gsl_vector* x, gsl_vector* y, gsl_vector* dy_dx) : BSpline(x, y) {
        // Ensure the input vectors have the same size
        if (x->size != y->size || x->size != dy_dx->size)
            throw invalid_argument("Input vectors must have the same size");

        gsl_bspline_init_hermite(1, positions, workSpace);
        gsl_bspline_interp_chermite(positions, values, dy_dx, coefficients, workSpace);
    }
};

#endif //SPLINE_H_