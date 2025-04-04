#include "BSpline.h"

namespace getelec{


void BSpline::setPositions(const vector<double>& x){
    if (positions)
        gsl_vector_free(positions);
    positions = gsl_vector_alloc(x.size());
    for (size_t i = 0; i < x.size(); i++) 
        gsl_vector_set(positions, i, x[i]); 
}

void BSpline::setPositionsAndValues(const vector<double>& x, const vector<double>& y){

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

void BSpline::setPositionsAndValues(gsl_vector* x, gsl_vector* y){
    // Ensure the input vectors have the same size
    assert(x->size == y->size && "Input vectors must have the same size");
    positions = x;
    values = y;
}

double BSpline::evaluate(double x) const {
    assert(workSpace && coefficients && "The spline has not been initialized");
    double result;
    gsl_bspline_calc(x, coefficients, &result, workSpace);
    return result;
}


vector<double> BSpline::evaluateMultiple(double x) const{
    vector<double> result(coefficientSets.size());
    for (size_t i = 0; i < coefficientSets.size(); i++){
        assert(workSpace && coefficientSets[i] && "The spline has not been initialized");
        gsl_bspline_calc(x, coefficientSets[i], &result[i], workSpace);
    }
    return result;
}


double BSpline::evaluateDerivative(double x, size_t nDerivative) const {
    assert(workSpace && coefficients && "The spline has not been initialized");
    double result;
    gsl_bspline_calc_deriv(x, coefficients, nDerivative, &result, workSpace);
    return result;
}

vector<double> BSpline::evaluateDerivativeMultiple(double x, size_t nDerivative) const {
    vector<double> result(coefficientSets.size());
    for (size_t i = 0; i < coefficientSets.size(); i++){
        assert(workSpace && coefficientSets[i] && "The spline has not been initialized");
        gsl_bspline_calc_deriv(x, coefficientSets[i], nDerivative, &result[i], workSpace);
    }
    return result;
}

double BSpline::evaluateIntegral(double a, double b) const{
    assert(workSpace != nullptr && coefficients != nullptr);

    double result;
    gsl_bspline_calc_integ(a, b, coefficients, &result, workSpace);
    return result;
}

vector<double> BSpline::evaluateIntegralMultiple(double a, double b) const{
    vector<double> result(coefficientSets.size());
    for (size_t i = 0; i < coefficientSets.size(); i++){
        assert(workSpace && coefficientSets[i] && "The spline has not been initialized");
        gsl_bspline_calc_integ(a, b, coefficientSets[i], &result[i], workSpace);
    }
    return result;
}

void CubicHermiteSpline::initialize(const vector<double>& x, const vector<double>& y, const vector<double>& dy_dx) {
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

void CubicHermiteSpline::initialize(gsl_vector* x, gsl_vector* y, gsl_vector* dy_dx) {
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


void CubicHermiteSpline::initializeMultiple(const vector<double>& x, const vector<vector<double>>& ySets, const vector<vector<double>>& dy_dxSets){
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

void CubicHermiteSpline::reallocateWorkSpace(){
    size_t nControlPoints = 2 * positions->size;
    if (workSpace)
        gsl_bspline_free(workSpace);
    workSpace = gsl_bspline_alloc_ncontrol(4, nControlPoints);
}

void CubicHermiteSpline::reallocateCoefficients(){
    size_t nControlPoints = 2 * positions->size;
    if (coefficients)
        gsl_vector_free(coefficients);

    for (auto coeffs : coefficientSets)
        if (coeffs) gsl_vector_free(coeffs);

    coefficientSets.clear();
    coefficients = gsl_vector_alloc(nControlPoints);  
}

void CubicHermiteSpline::allocateCoefficients(){
    size_t nControlPoints = 2 * positions->size;
    coefficients = gsl_vector_alloc(nControlPoints);
}

void CubicHermiteSpline::reallocate(){
    reallocateWorkSpace();
    reallocateCoefficients();
}

void CubicHermiteSpline::intitialize(){
    gsl_bspline_init_hermite(1, positions, workSpace);
    gsl_bspline_interp_chermite(positions, values, derivs, coefficients, workSpace);
}



}