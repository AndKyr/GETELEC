#ifndef TRANSMISSIONCALCULATOR_H_
#define TRANSMISSIONCALCULATOR_H_

#include <math.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <vector>

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
        return CONSTANTS.kConstant * (energy - potentialFunction(x));
    }

    void setEnergy(double E){energy = E;}

    double getEnergy(){return energy;}
};

class ModifiedSNBarrier : public TunnelingFunctionBase{
private:
    double radius = 1.e5;
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
        return field * ((gamma*gamma - 2*gamma + 1) * radius*radius + (2*gamma - 2) * radius*z + gamma*z*z) / 
        (((gamma - 1)*radius + gamma*z, 2) * ((gamma - 1)*radius + gamma*z, 2));
    }

public:
    ModifiedSNBarrier(){}
    ModifiedSNBarrier(double f, double r, double g) : field(f), radius(r), gamma(g){}
    ~ModifiedSNBarrier();

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

    double potentialFunction(double z){
        return -imagePotential(z) - electrostaticPotential(z);
    }

    double potentialFunctionDerivative(double z){
        return -imagePotentialDerivative(z) - electrostaticPotentialDerivative(z);
    }
    
};

int differentialSystem(double x, const double y[], double f[], void *params){

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;
    
    f[0] =  -tunnelingFunction->kappaSquared(x) - y[0]*y[0] + y[1]*y[1];
    f[1] = - 2. * y[0] * y[1];
    return GSL_SUCCESS;
}

int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params){

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *matrix = &dfdy_mat.matrix;
    gsl_matrix_set(matrix, 0, 0, -2 * y[0]);
    gsl_matrix_set(matrix, 0, 1, 2 * y[1]);
    gsl_matrix_set(matrix, 1, 0, -2.0*y[1]);
    gsl_matrix_set(matrix, 1, 1, - 2 * y[0]);
    dfdt[0] = CONSTANTS.kConstant * tunnelingFunction->potentialFunctionDerivative(x);
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

class TransmissionCalculator{
private:

    vector<double> xLimits = {1.e-5, 10.}; //limits within which the potential function is defined
    double kAtLimits[2]; //Values of k vector at the potential edges

    TunnelingFunctionBase tunnelingFunction = ModifiedSNBarrier();

    vector<double> solutionRealPart;
    vector<double> solutionImaginaryPart;
    int maxIntegrationSteps = 1024;

public:
    TransmissionCalculator();
    TransmissionCalculator(int max_steps, vector<double> xLims, const TunnelingFunctionBase& tunnelFunction) : maxIntegrationSteps(max_steps), xLimits(xLims), tunnelingFunction(tunnelFunction){}


    void calculateEdgeWavevectors(){
        kAtLimits[0] = sqrt(tunnelingFunction.kappaSquared(xLimits[0]));
        kAtLimits[1] = sqrt(tunnelingFunction.kappaSquared(xLimits[1]));
    }


    int solveDifferentialSystem(){
        const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rkf45;

        gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(stepType, 2);
        gsl_odeiv2_control *controller = gsl_odeiv2_control_y_new(1e-6, 0.0);
        gsl_odeiv2_evolve *evolver = gsl_odeiv2_evolve_alloc(2);

        gsl_odeiv2_system sys = {differentialSystem, differentialSystemJacobian, 2, NULL};

        double x = xLimits[1];
        double dx = -1e-2;
        double solutionVector[2] = {1.0, 0.0};

        for(unsigned i = 0; i < maxIntegrationSteps; i++){
            int status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xLimits[0], &dx, solutionVector);

            if (status != GSL_SUCCESS || x <= xLimits[0])
                break;

            printf ("%g %g %g %g %g %g\n", x, dx, solutionVector[0], solutionVector[1], evolver->dydt_out[0], evolver->dydt_out[1]);
        }

        gsl_odeiv2_evolve_free(evolver);
        gsl_odeiv2_control_free(controller);
        gsl_odeiv2_step_free(step);
        return 0;
    }

};  


int main(){

    TransmissionCalculator calculator;
    calculator.solveDifferentialSystem();
    return 0;

}

#endif /* TRANSMISSIONCALCULATOR */
