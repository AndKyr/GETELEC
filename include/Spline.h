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

class HermiteSpline{
public:
    HermiteSpline(vector<double>& x, vector<double>& y, vector<double>& dy_dx){

        // Ensure the input vectors have the same size
        if (x.size() != y.size() || x.size() != dy_dx.size())
            throw invalid_argument("Input vectors must have the same size");

        positions = gsl_vector_alloc(x.size());
        for (size_t i = 0; i < x.size(); i++) gsl_vector_set(positions, i, x[i]); 

        gsl_vector* values = gsl_vector_alloc(y.size());
        for (size_t i = 0; i < y.size(); i++) gsl_vector_set(values, i, y[i]); 

        gsl_vector* derivs = gsl_vector_alloc(y.size());
        for (size_t i = 0; i < dy_dx.size(); i++) gsl_vector_set(derivs, i, dy_dx[i]); 

        gsl_bspline_init_hermite(1, positions, workSpace);
        gsl_bspline_interp_chermite(positions, values, derivs, coefficients, workSpace);

    }

    ~HermiteSpline() {
        gsl_bspline_free(workSpace);
        gsl_vector_free(coefficients);
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

private:
    gsl_bspline_workspace* workSpace;  // Workspace for spline evaluation
    gsl_vector* coefficients;
    gsl_vector* positions;
};

#endif //SPLINE_H_