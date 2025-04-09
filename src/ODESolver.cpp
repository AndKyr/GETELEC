#include "ODESolver.h"
#include <gsl/gsl_math.h>
#include <iostream>


namespace getelec{

int ODESolver::solve(bool saveSolution){

    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();

    if (saveSolution){
        savedSolution.push_back(solutionVector);
        xSaved.push_back(x);
    }
    
    for (size_t i = 0; i < maxAllowedSteps; i++){ //loop over blocks

        dx = GSL_SIGN(dx) *  min(abs(maxStepSize), abs(dx));
        status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());

        if (saveSolution){
            savedSolution.push_back(solutionVector);
            xSaved.push_back(x);
        } 

        if (x == xFinal || status != GSL_SUCCESS)    
            return status;         
    }
    return GSL_CONTINUE;
}


const gsl_odeiv2_step_type* ODESolver::getStepTypeFromString(const string& stepTypeName){
    needsJacobian = false;
    if (stepTypeName == "rk2")
        return gsl_odeiv2_step_rk2;
    else if (stepTypeName == "rkf45")
        return gsl_odeiv2_step_rk4;
    else if (stepTypeName == "rkf45")
        return gsl_odeiv2_step_rkf45;
    else if (stepTypeName == "rkck")
        return gsl_odeiv2_step_rkck;
    else if (stepTypeName == "rk8pd")
        return gsl_odeiv2_step_rk8pd;
    
    needsJacobian = true;
    if (stepTypeName == "rk2imp")
        return gsl_odeiv2_step_rk2imp;
    else if (stepTypeName == "rk4imp")
        return gsl_odeiv2_step_rk4imp;
    else if (stepTypeName == "bsimp")
        return gsl_odeiv2_step_bsimp;
    else if (stepTypeName == "rk1imp")
        return gsl_odeiv2_step_rk1imp;
    else if (stepTypeName == "msadams")
        return gsl_odeiv2_step_msadams;
    else if (stepTypeName == "msbdf")
        return gsl_odeiv2_step_msbdf;   
    else {
        cout << "WARNING : stpeTypeName " << stepTypeName << "not valid. Reverting to default rkck" << endl;
        needsJacobian = false;
        return gsl_odeiv2_step_rkck;
    }
}

int ODESolver::solveNoSave(){

    double x = xInitial;
    double dx = initialStep;
    int status;
    reinitialize();

    for (size_t i = 0; i < maxAllowedSteps; i++){ //loop over blocks
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
        for (auto &y : savedSolution[i]) 
            outFile << y << " ";
        outFile << endl; 
    }
}

} // namespace getelec