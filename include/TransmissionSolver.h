#ifndef TRANSMISSIONSOLVER_H_
#define TRANSMISSIONSOLVER_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>

#include "ODESolver.h"
#include "TunnelingFunction.h"
#include "Utilities.h"

using namespace std;


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
            int stepExpectedForInitialStep = 64,
            double maxPotentialDepth = 10.
        );

    void setXlimits(double maxPotentialDepth){
        xInitial = tunnelingFunction->findRightXLimit(maxPotentialDepth);
        xFinal = tunnelingFunction->findLeftXLimit(maxPotentialDepth);
    }

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

    void printXLimits(){ cout << "xInitial = " << xInitial << " xFinal = " << xFinal << endl;}
};

#endif /* TRANSMISSIONSOLVER_H_ */
