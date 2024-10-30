#ifndef TRANSMISSIONCALCULATOR_H_
#define TRANSMISSIONCALCULATOR_H_

#include <math.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_interp.h>

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


int differentialSystem2D(double x, const double y[], double f[], void *params);

int differentialSystem3D(double x, const double y[], double f[], void *params);

int differentialSystem4D(double x, const double y[], double f[], void *params);


int differentialSystemJacobian2D(double x, const double y[], double *dfdy, double dfdt[], void *params);

int differentialSystemJacobian3D(double x, const double y[], double *dfdy, double dfdt[], void *params);

int differentialSystemJacobian4D(double x, const double y[], double *dfdy, double dfdt[], void *params);

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

class TransmissionCalculator{
private:

    vector<double> xLimits = {0.03599847, 2.00400712}; //limits within which the potential function is defined
    vector<double> kAtLimits = {1., 1.}; //Values of k vector at the potential edges

    int systemDimension = 3;
    double relativeTolerance = 1.e-3;
    double absoluteTolerance = 1.e-1;
    double initialStep;

    TunnelingFunctionBase* tunnelingFunction;

    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rk4;
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(stepType, systemDimension);
    gsl_odeiv2_control *controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    gsl_odeiv2_evolve *evolver = gsl_odeiv2_evolve_alloc(systemDimension);
    gsl_odeiv2_system sys;

    vector<double> solutionVector = vector<double>(systemDimension, 0.0);

public:
    
    void updateKappaAtLimits(){
        kAtLimits[0] = sqrt(tunnelingFunction->kappaSquared(xLimits[0]));
        kAtLimits[1] = sqrt(tunnelingFunction->kappaSquared(xLimits[1]));
    }

    void setTolerances(double absoluteTolerance = 1.e-5, double relativeTolerance = 1.e-5){
        relativeTolerance = relativeTolerance;
        absoluteTolerance = absoluteTolerance;
        controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    }


    TransmissionCalculator(TunnelingFunctionBase* tunnelFunctionPtr, 
                            int systemDimension = 3, 
                            vector<double> xLims = {2.00400712, 0.03599847},
                            double relativeTolerance = 1.e-4,
                            double absoluteTolerance = 1.e-4,
                            const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd
                        );


    ~TransmissionCalculator(){
        gsl_odeiv2_evolve_free(evolver);
        gsl_odeiv2_control_free(controller);
        gsl_odeiv2_step_free(step);
    }


    int solveDifferentialSystem();

    double transmissionCoefficient() const;

};  



#endif /* TRANSMISSIONCALCULATOR */
