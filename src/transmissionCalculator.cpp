#include "transmissionCalculator.h"


gsl_vector* vector2gsl(vector<double>& vec){
    gsl_vector_view gsl_vec_view = gsl_vector_view_array(vec.data(), vec.size());
    return &gsl_vec_view.vector;
}

int differentialSystem(double x, const double y[], double f[], void *params){

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;
    
    f[0] =  -tunnelingFunction->kappaSquared(x) - y[0]*y[0] + y[1]*y[1];
    f[1] = - 2. * y[0] * y[1];
    return GSL_SUCCESS;
}

int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params){
    ///< TODO: It seems that this jacobian is incorrect. Test and fix it.

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

    gsl_odeiv2_system sys = {differentialSystem, NULL, 2, tunnelingFunction};

    double x = xLimits[1];
    double dx = -1e-1;
    double solutionVector[2] = {0, sqrt(tunnelingFunction->kappaSquared(xLimits[1]))};
    double solutionDerivativeVector[2];

    solution.resize(blockSize); // allocate the first block

    //save the initial values
    solution.position[0] = x;
    solution.realPart[0] = solutionVector[0];
    solution.imaginaryPart[0] = solutionVector[1];

    //calculate the derivative at the initial points
    differentialSystem(x, solutionVector, solutionDerivativeVector, tunnelingFunction);

    //save the derivatives
    solution.realPartDerivative[0] = solutionDerivativeVector[0];
    solution.imaginaryPartDerivative[0] = solutionDerivativeVector[1];
    
    for (unsigned iStep = 1; iStep < maxIntegrationSteps; iStep++){ //loop over blocks
        int status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xLimits[0], &dx, solutionVector);

        solution.position[iStep] = x;
        solution.realPart[iStep] = solutionVector[0];
        solution.imaginaryPart[iStep] = solutionVector[1];
        solution.realPartDerivative[iStep] = evolver->dydt_out[0];
        solution.imaginaryPartDerivative[iStep] = evolver->dydt_out[1];
        
        if (status != GSL_SUCCESS || x <= xLimits[0]){
            solution.resize(iStep);
            return status;
        }

        if (iStep >= solution.length -1)
            solution.resize(solution.length + blockSize); //resize the vectors by another block
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


void TransmissionCalculator::writeSolution(string fileName){
    ofstream outFile(fileName, ios::out);

    outFile << "i    x    Re{s'}    Im{s'}    Re{s''}    Im{s''}   kappa^2" << endl;

    for (size_t i = 0; i < solution.length; i++)
        outFile << i << " " << solution.position[i] << " " << solution.realPart[i] << " " 
                << solution.imaginaryPart[i] << " " << solution.realPartDerivative[i] << " " << 
                    solution.imaginaryPartDerivative[i] << " " << tunnelingFunction->kappaSquared(solution.position[i]) << endl;
}