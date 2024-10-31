#ifndef TRANSMISSIONCALCULATOR_H_
#define TRANSMISSIONCALCULATOR_H_

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

#include <vector>
#include <typeinfo>
#include <fstream>
#include <string>

using namespace std;

static constexpr struct PhysicalConstants{
    double hbarSqrOver2m = 3.80998212e-2;
    double kConstant = 1./hbarSqrOver2m;
    double imageChargeConstant = .359991137;

} CONSTANTS;

class TunnelingFunctionBase{
private:
    double energy = 0.;

public:
    TunnelingFunctionBase(){}
    TunnelingFunctionBase(double E) : energy(E){}
    ~TunnelingFunctionBase(){}

    virtual double potentialFunction(double x){
        return 0. * x;
    }

    virtual double potentialFunctionDerivative(double x){
        return 0. * x;
    }

    double kappaSquared(double x){
        return CONSTANTS.kConstant * (energy - this->potentialFunction(x));
    }

    void setEnergy(double E){energy = E;}

    double getEnergy(){return energy;}
};

class ModifiedSNBarrier : public TunnelingFunctionBase{
private:
    double radius = 1.e3;
    double field = 5.;
    double gamma = 10.;

    double imagePotential(double z){
        return CONSTANTS.imageChargeConstant * ( 1. / z - .5 / radius);
    }

    double imagePotentialDerivative(double z){
        return -CONSTANTS.imageChargeConstant * (1./(z * z));
    }

    double electrostaticPotential(double z){
        return field * (radius * (gamma - 1.) * z + z*z) / (gamma * z + radius * (gamma - 1.));
    }

    double electrostaticPotentialDerivative(double z){
        return field * ((gamma-1)*(gamma-1) * radius*radius + 2*(gamma-1)*radius*z + gamma*z*z) / 
        (((gamma - 1)*radius + gamma*z) * ((gamma - 1)*radius + gamma*z));
    }

public:
    ModifiedSNBarrier(){}
    ModifiedSNBarrier(double f, double r, double g) : field(f), radius(r), gamma(g){}
    ~ModifiedSNBarrier(){}

    void setBarrierParameters(double f, double r, double g){
        field = f;
        radius = r;
        gamma = g;
    } 

    void setField(double f){
        field = f;
    }

    void setRadius(double R){
        radius = R;
    }

    void setGamma(double g){
        gamma = g;
    }

    double potentialFunction(double z) override{
        return -imagePotential(z) - electrostaticPotential(z);
    }

    double potentialFunctionDerivative(double z) override{
        return -imagePotentialDerivative(z) - electrostaticPotentialDerivative(z);
    }
    
};

/*
class HermiteSpline{
public:
    HermiteSpline(vector<double>& x, vector<double>& y, vector<double>& dy_dx) : positions(vector2gsl(x)){
        // Ensure the input vectors have the same size
        if (x.size() != y.size() || x.size() != dy_dx.size())
            throw std::invalid_argument("Input vectors must have the same size");

        gsl_bspline_init_hermite(1, positions, workSpace);
        gsl_bspline_interp_chermite(vector2gsl(x), vector2gsl(y), vector2gsl(dy_dx), coefficients, workSpace);

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
*/

class ODESolver{
protected:
    int systemDimension;
    double xInitial;
    double xFinal;
    double relativeTolerance;
    double absoluteTolerance;
    int maxAllowedSteps;
    int stepsExpectedForInitialStep = 64;
    double initialStep;
    
    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rk4;
    gsl_odeiv2_step *step;
    gsl_odeiv2_control *controller;
    gsl_odeiv2_evolve *evolver;
    gsl_odeiv2_system sys;
    vector<double> initialValues;
    vector<vector<double>> savedSolution;
    vector<double>xSaved;
    void* parameters;

    int (*differentialSystem)(double, const double*, double* , void*);
    int (*differentialSystemJacobian)(double, const double*, double*, double*, void*);

    vector<double> solutionVector = vector<double>(systemDimension, 0.0);

public:

    ODESolver( 
                vector<double>initialValues,
                int (*differentialSystem)(double, const double*, double*, void*),
                int systemDimension = 3, 
                vector<double> xLims = {0., 1.},
                double relativeTolerance = 1.e-4,
                double absoluteTolerance = 1.e-4,
                const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd,
                int maxSteps = 4096,
                int stepExpectedForInitialStep = 64,
                int (*differentialSystemJacobian)(double, const double*, double*, double*, void*) = NULL,
                void* params = NULL
            );

    ~ODESolver(){
        gsl_odeiv2_evolve_free(evolver);
        gsl_odeiv2_control_free(controller);
        gsl_odeiv2_step_free(step);
    }

    void setTolerances(double absoluteTolerance = 1.e-5, double relativeTolerance = 1.e-5){
        relativeTolerance = relativeTolerance;
        absoluteTolerance = absoluteTolerance;
        controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    }

    void setInitialValues(vector<double> initValues){initialValues = initValues;}

    void reinitialize(){
        solutionVector = initialValues;
        gsl_odeiv2_step_reset(step);
        gsl_odeiv2_evolve_reset(evolver);
        savedSolution.clear();
        xSaved.clear();
    }

    int solve(bool saveSolution = false);

    int solveNoSave();

    void writeSolution(string filename = "odeSolution.dat");

};


class TransmissionSolver : public ODESolver{
private:
    TunnelingFunctionBase* tunnelingFunction;
     
    double kappaInitial;
    double kappaFinal;

    static int tunnelingDifferentialSystem(double x, const double y[], double f[], void *params){
        TunnelingFunctionBase* barrier = (TunnelingFunctionBase*) params;
        f[0] =  -barrier->kappaSquared(x) - y[0]*y[0] + y[1]*y[1];
        f[1] = - 2. * y[0] * y[1];
        f[2] = y[0];
        return GSL_SUCCESS;
    }

    static int tunnelingSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params){
        TunnelingFunctionBase* barrier = (TunnelingFunctionBase*) params;
        gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
        gsl_matrix *matrix = &dfdy_mat.matrix;

        gsl_matrix_set(matrix, 0, 0, -2*y[0]); gsl_matrix_set(matrix, 0, 1, 2*y[1]);  gsl_matrix_set(matrix, 0, 2, 0.);
        gsl_matrix_set(matrix, 1, 0, -2*y[1]); gsl_matrix_set(matrix, 1, 1, -2*y[0]); gsl_matrix_set(matrix, 1, 2, 0.);
        gsl_matrix_set(matrix, 2, 0, 1.);      gsl_matrix_set(matrix, 2, 1, 0.);      gsl_matrix_set(matrix, 2, 2, 0.);

        dfdt[0] = CONSTANTS.kConstant * barrier->potentialFunctionDerivative(x);
        dfdt[1] = 0.0; dfdt[2] = 0.0;
        return GSL_SUCCESS;
    }

public:   

    TransmissionSolver(
            TunnelingFunctionBase* tunnelFunctionPtr, 
            double relativeTolerance = 1.e-4,
            double absoluteTolerance = 1.e-4,
            const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd,
            int maxSteps = 4096,
            int stepExpectedForInitialStep = 64
        );


    void updateKappaAtLimits(){
        kappaInitial = sqrt(tunnelingFunction->kappaSquared(xInitial));
        kappaFinal = sqrt(tunnelingFunction->kappaSquared(xFinal));
        initialValues = {0., kappaInitial, 0.};
    }

    double transmissionCoefficient() const{
        double CplusCoefficientSquared = 0.25 * exp(2. * solutionVector[2]) * (1. + 
            (solutionVector[0]*solutionVector[0] + solutionVector[1]*solutionVector[1]) / (kappaFinal*kappaFinal) +
                2. * solutionVector[1] / kappaFinal);

        return kappaInitial / (kappaFinal * CplusCoefficientSquared);
    }
};

#endif /* TRANSMISSIONCALCULATOR */
