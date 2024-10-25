#include "transmissionCalculator.h"

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


int TransmissionCalculator::solveDifferentialSystem(){
    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(stepType, 2);
    gsl_odeiv2_control *controller = gsl_odeiv2_control_y_new(1e-4, 0.0);
    gsl_odeiv2_evolve *evolver = gsl_odeiv2_evolve_alloc(2);

    gsl_odeiv2_system sys = {differentialSystem, differentialSystemJacobian, 2, tunnelingFunction};

    double x = xLimits[1];
    double dx = -1e-1;
    double solutionVector[2] = {0, sqrt(tunnelingFunction->kappaSquared(xLimits[1]))};

    unsigned iStep = 0;
    for (unsigned i = 0; i < maxIntegrationSteps / blockSize; i++){
        //resize the vectors by a block
        solution.resize(solution.length + blockSize);

        for(unsigned j = 0; j < blockSize; j++){
            int status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xLimits[0], &dx, solutionVector);

            solution.realPart[iStep] = solutionVector[0];
            solution.imaginaryPart[iStep] = solutionVector[1];
            solution.realPartDerivative[iStep] = evolver->dydt_out[0];
            solution.imaginaryPartDerivative[iStep] = evolver->dydt_out[1];
            iStep++;
            
            if (status != GSL_SUCCESS || x <= xLimits[0]){
                solution.resize(iStep);
                return status;
            }

        }
    }

    return GSL_CONTINUE;
}

void TransmissionCalculator::Solution::resize(size_t N){
    position.resize(N);
    realPart.resize(N);
    imaginaryPart.resize(N);
    realPartDerivative.resize(N);
    imaginaryPartDerivative.resize(N);
    length = N;
}

void TransmissionCalculator::Solution::write(string fileName){
    ofstream outFile(fileName, ios::out);

    for (size_t i = 0; i < length; i++)
        outFile << i << " " << position[i] << " " << realPart[i] << " " << imaginaryPart[i] << " " << 
                        realPartDerivative[i] << " " << imaginaryPartDerivative[i] << endl;
}
