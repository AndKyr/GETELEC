#include "transmissionCalculator.h"


ODESolver::ODESolver( 
        vector<double>initialValues,
        int (*differentialSystem)(double, const double*, double* , void*),
        int systemDimension, 
        vector<double> xLims,
        double rtol,
        double atol,
        const gsl_odeiv2_step_type* stepType,
        int maxSteps,
        int stepExpectedForInitialStep,
        int (*differentialSystemJacobian)(double, const double*, double*, double*, void*),
        void* params
                ) : systemDimension(systemDimension),
                    xInitial(xLims[0]),
                    xFinal(xLims[1]),
                    stepType(stepType),
                    maxAllowedSteps(maxSteps),
                    stepsExpectedForInitialStep(stepExpectedForInitialStep),
                    initialValues(initialValues),
                    differentialSystem(differentialSystem),
                    differentialSystemJacobian(differentialSystemJacobian),
                    parameters(params),
                    step(gsl_odeiv2_step_alloc(stepType, systemDimension)),
                    evolver(gsl_odeiv2_evolve_alloc(systemDimension))
{
    initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;
    sys = {differentialSystem, differentialSystemJacobian, (unsigned long) systemDimension, parameters};
    setTolerances(atol, rtol);
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
    return GSL_CONTINUE;
}

int ODESolver::solveNoSave(){

    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();

    for (size_t i; i < maxAllowedSteps; i++){ //loop over blocks
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());
        if (x == xFinal || status != GSL_SUCCESS)    
            return status;         
    }
    return GSL_CONTINUE;

}

void ODESolver::writeSolution(string filename){
    ofstream outFile(filename, ios::out);        
    for (size_t i = 0; i < xSaved.size(); i++){
        outFile << xSaved[i] << " ";
        for (auto &y : savedSolution[i]) outFile << y << " ";
        outFile << endl; 
    }
}

TransmissionSolver::TransmissionSolver(TunnelingFunctionBase* tunnelFunctionPtr, double rtol, double atol, const gsl_odeiv2_step_type* stepType,
                         int maxSteps, int stepExpectedForInitialStep
                        ) : ODESolver(vector<double>(3, 0.0), tunnelingDifferentialSystem, 3, {2.00400712, 0.03599847},
                                        rtol, atol, stepType, maxSteps, stepExpectedForInitialStep, tunnelingSystemJacobian, tunnelFunctionPtr),
                            tunnelingFunction(tunnelFunctionPtr)
                            

{
    updateKappaAtLimits();
}