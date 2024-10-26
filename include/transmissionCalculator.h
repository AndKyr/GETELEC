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


gsl_vector* vector2gsl(vector<double>& vec);


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


class TransmissionCalculator{
private:

    vector<double> xLimits = {1.e-1, 4.}; //limits within which the potential function is defined
    // double kAtLimits[2]; //Values of k vector at the potential edges

    const int blockSize = 128; ///< Length of the vector block that the system solution will be saved.

    TunnelingFunctionBase* tunnelingFunction;

    struct Solution{
        size_t length = 0;
        vector<double> position;
        vector<double> realPart;
        vector<double> imaginaryPart;  
        vector<double> realPartDerivative;
        vector<double> imaginaryPartDerivative;
        void resize(size_t N);
        void write(string filename = "differentialSystemSolution.dat");
    } solution;


    int maxIntegrationSteps = 1024;

public:
    TransmissionCalculator(){}
    TransmissionCalculator(TunnelingFunctionBase* tunnelFunctionPtr) : tunnelingFunction(tunnelFunctionPtr){}
    ~TransmissionCalculator(){}

    int solveDifferentialSystem();

    int getSolutionSpline();

    void writeSolution(string filename = "differentialSystemSolution.dat");

};  


int differentialSystem(double x, const double y[], double f[], void *params);

int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

#endif /* TRANSMISSIONCALCULATOR */
