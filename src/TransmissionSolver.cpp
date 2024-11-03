#include "TransmissionSolver.h"

TransmissionSolver::TransmissionSolver(TunnelingFunction* tunnelFunctionPtr, double rtol, double atol, const gsl_odeiv2_step_type* stepType,
                         int maxSteps, int stepExpectedForInitialStep, double maxPotentialDepth
                        ) : ODESolver(vector<double>(3, 0.0), tunnelingDifferentialSystem, 3, {2.00400712, 0.03599847},
                                        rtol, atol, stepType, maxSteps, stepExpectedForInitialStep, tunnelingSystemJacobian, tunnelFunctionPtr),
                            tunnelingFunction(tunnelFunctionPtr)
                            

{
    setXlimits(maxPotentialDepth);
    updateKappaAtLimits();
}

int TransmissionSolver::tunnelingDifferentialSystem(double x, const double y[], double f[], void *params){
        TunnelingFunction* barrier = (TunnelingFunction*) params;
        f[0] =  -barrier->kappaSquared(x) - y[0]*y[0] + y[1]*y[1];
        f[1] = - 2. * y[0] * y[1];
        f[2] = y[0];
        return GSL_SUCCESS;
}

int TransmissionSolver::tunnelingSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params){
        TunnelingFunction* barrier = (TunnelingFunction*) params;
        gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
        gsl_matrix *matrix = &dfdy_mat.matrix;

        gsl_matrix_set(matrix, 0, 0, -2*y[0]); gsl_matrix_set(matrix, 0, 1, 2*y[1]);  gsl_matrix_set(matrix, 0, 2, 0.);
        gsl_matrix_set(matrix, 1, 0, -2*y[1]); gsl_matrix_set(matrix, 1, 1, -2*y[0]); gsl_matrix_set(matrix, 1, 2, 0.);
        gsl_matrix_set(matrix, 2, 0, 1.);      gsl_matrix_set(matrix, 2, 1, 0.);      gsl_matrix_set(matrix, 2, 2, 0.);

        dfdt[0] = CONSTANTS.kConstant * barrier->potentialFunctionDerivative(x);
        dfdt[1] = 0.0; dfdt[2] = 0.0;
        return GSL_SUCCESS;
}

double TransmissionSolver::transmissionCoefficient() const{
    double CplusCoefficientSquared = 0.25 * exp(2. * solutionVector[2]) * (1. + 
        (solutionVector[0]*solutionVector[0] + solutionVector[1]*solutionVector[1]) / (kappaFinal*kappaFinal) +
            2. * solutionVector[1] / kappaFinal);

    return kappaInitial / (kappaFinal * CplusCoefficientSquared);
}