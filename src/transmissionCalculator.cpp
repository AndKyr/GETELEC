#include "transmissionCalculator.h"


// gsl_vector* vector2gsl(vector<double>& vec){
//     gsl_vector_view gsl_vec_view = gsl_vector_view_array(vec.data(), vec.size());
//     return &gsl_vec_view.vector;
// }

int differentialSystem2D(double x, const double y[], double f[], void *params){

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;
    
    f[0] =  -tunnelingFunction->kappaSquared(x) - y[0]*y[0] + y[1]*y[1];
    f[1] = - 2. * y[0] * y[1];
    return GSL_SUCCESS;
}

int differentialSystem3D(double x, const double y[], double f[], void *params){

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;
    
    f[0] =  -tunnelingFunction->kappaSquared(x) - y[0]*y[0] + y[1]*y[1];
    f[1] = - 2. * y[0] * y[1];
    f[2] = y[0];
    return GSL_SUCCESS;
}

int differentialSystem4D(double x, const double y[], double f[], void *params){

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;
    
    f[0] =  -tunnelingFunction->kappaSquared(x) - y[0]*y[0] + y[1]*y[1];
    f[1] = - 2. * y[0] * y[1];
    f[2] = y[0];
    f[3] = y[1];
    return GSL_SUCCESS;
}

int differentialSystemJacobian2D(double x, const double y[], double *dfdy, double dfdt[], void *params){
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

int differentialSystemJacobian3D(double x, const double y[], double *dfdy, double dfdt[], void *params){
    ///< TODO: It seems that this jacobian is incorrect. Test and fix it.

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
    gsl_matrix *matrix = &dfdy_mat.matrix;
    gsl_matrix_set(matrix, 0, 0, -2*y[0]); gsl_matrix_set(matrix, 0, 1, 2*y[1]);  gsl_matrix_set(matrix, 0, 2, 0.);
    gsl_matrix_set(matrix, 1, 0, -2*y[1]); gsl_matrix_set(matrix, 1, 1, -2*y[0]); gsl_matrix_set(matrix, 1, 2, 0.);
    gsl_matrix_set(matrix, 2, 0, 1.);      gsl_matrix_set(matrix, 2, 1, 0.);      gsl_matrix_set(matrix, 2, 2, 0.);

    dfdt[0] = CONSTANTS.kConstant * tunnelingFunction->potentialFunctionDerivative(x);
    dfdt[1] = 0.0; dfdt[2] = 0.0;
    return GSL_SUCCESS;
}

int differentialSystemJacobian4D(double x, const double y[], double *dfdy, double dfdt[], void *params){
    ///< TODO: It seems that this jacobian is incorrect. Test and fix it.

    TunnelingFunctionBase* tunnelingFunction = (TunnelingFunctionBase*) params;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
    gsl_matrix *matrix = &dfdy_mat.matrix;
    gsl_matrix_set(matrix, 0, 0, -2*y[0]); gsl_matrix_set(matrix, 0, 1, 2*y[1]);  gsl_matrix_set(matrix, 0, 2, 0.); gsl_matrix_set(matrix, 0, 3, 0.);
    gsl_matrix_set(matrix, 1, 0, -2*y[1]); gsl_matrix_set(matrix, 1, 1, -2*y[0]); gsl_matrix_set(matrix, 1, 2, 0.); gsl_matrix_set(matrix, 1, 3, 0.);
    gsl_matrix_set(matrix, 2, 0, 1.);      gsl_matrix_set(matrix, 2, 1, 0.);      gsl_matrix_set(matrix, 2, 2, 0.); gsl_matrix_set(matrix, 2, 3, 0.);
    gsl_matrix_set(matrix, 3, 0, 0.);      gsl_matrix_set(matrix, 3, 1, 1.);      gsl_matrix_set(matrix, 3, 2, 0.); gsl_matrix_set(matrix, 3, 3, 0.);


    dfdt[0] = CONSTANTS.kConstant * tunnelingFunction->potentialFunctionDerivative(x);
    dfdt[1] = 0.0; dfdt[2] = 0.0; dfdt[3] = 0.0;

    return GSL_SUCCESS;
}


int ODESolver::solve(bool saveSolution){

    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();

    for (size_t i; i < maxAllowedSteps; i++){ //loop over blocks
        if (saveSolution){
            savedSolution.push_back(solutionVector);
            xSaved.push_back(x);
        } 

        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());

        if (x == xFinal || status != GSL_SUCCESS)    
            return status;         
    }
}



TransmissionCalculator::TransmissionCalculator( 
                            TunnelingFunctionBase* tunnelFunctionPtr, int systemDimension, 
                            vector<double> xLims, double relativeTolerance, 
                            double absoluteTolerance, const gsl_odeiv2_step_type* stepType
                        ) : tunnelingFunction(tunnelFunctionPtr),
                            systemDimension(systemDimension),
                            xLimits(xLims),
                            relativeTolerance(relativeTolerance),
                            absoluteTolerance(absoluteTolerance),
                            controller(gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance)),
                            step(gsl_odeiv2_step_alloc(stepType, systemDimension)),
                            evolver(gsl_odeiv2_evolve_alloc(systemDimension)),
                            stepType(stepType),
                            initialStep((xLimits[1] - xLimits[0]) * 1.e-3)
{
    if (systemDimension == 2)
        sys = {differentialSystem2D, differentialSystemJacobian2D, 2, tunnelingFunction};
    else if (systemDimension == 3)
        sys = {differentialSystem3D, differentialSystemJacobian3D, 3, tunnelingFunction};
    else if (systemDimension == 4)
        sys = {differentialSystem4D, differentialSystemJacobian4D, 4, tunnelingFunction};
    else
        throw invalid_argument("systemDimension must be 2, 3, or 4");

    updateKappaAtLimits();
}


int TransmissionCalculator::solveDifferentialSystem(){

    double x = xLimits[0];
    double dx = initialStep;
    int status;
    solutionVector = vector<double>(systemDimension, 0.0);
    solutionVector[1] = kAtLimits[0];
    

    #ifndef RELEASE
        ofstream outFile("differentialSystemSolution.dat", ios::out);        
    #endif

    while(x > xLimits[1]){ //loop over blocks
        #ifndef RELEASE
            outFile << x << " ";
            for (auto &y : solutionVector) outFile << y << " ";
            outFile << tunnelingFunction->kappaSquared(x) << endl;       
        #endif

        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xLimits[1], &dx, solutionVector.data());             
    }

    #ifndef RELEASE
        outFile << x << " ";
        for (auto &y : solutionVector) outFile << y << " ";
            outFile << tunnelingFunction->kappaSquared(x) << endl;       
    #endif

    gsl_odeiv2_step_reset(step);
    gsl_odeiv2_evolve_reset(evolver);

    return status;
}


double TransmissionCalculator::transmissionCoefficient() const{

    double CplusCoefficientSquared = 0.25 * exp(2. * solutionVector[2]) * (1. + 
        (solutionVector[0]*solutionVector[0] + solutionVector[1]*solutionVector[1]) / (kAtLimits[1]*kAtLimits[1]) +
            2. * solutionVector[1] / kAtLimits[1]);

    return kAtLimits[0] / (kAtLimits[1] * CplusCoefficientSquared);
}