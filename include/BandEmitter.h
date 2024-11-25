#ifndef BANDEMITTER_H_
#define BANDEMITTER_H_

#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include "TransmissionSolver.h"

// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_odeiv2.h>

// #include <vector>
// #include <typeinfo>
// #include <fstream>
// #include <string>


using namespace std;



class BandEmitter : public ODESolver{
private:


    /**
     * @brief A struct that keeps all the important parameters of the emitter. It is necessary so that a void pointer is passed to the
     * differentialSystem function as params, to avoid passing "this" object itself (ugly self-reference)
     */
    double workFunction = 4.5;
    double kT = 0.025;
    double effectiveMass = 1.;
    double bandDepth = 10.; // the depth of the band as measured from the Fermi level
    vector<double> savedEnergies;
    vector<double> savedLogD;
    TransmissionSolver& transmissionSolver;
    TransmissionInterpolator interpolator;

    double calculateIntegrand(double energy){
        
        double result = interpolator.evaluate(energy);
        if (effectiveMass != 1.){
            double aBarX = - bandDepth + (1. - effectiveMass) * (energy + bandDepth);
            result -= (1. - effectiveMass) * interpolator.evaluate(aBarX);
        }
        return result;
    }

    static int differentialSystem(double energy, const double y[], double f[], void *params){
        BandEmitter* sysParams = (BandEmitter*) params; //cast the void pointer as SystemParams

        //calculate transmission coefficient
        double D = sysParams->calculateIntegrand(energy);
        //derivatives
        f[0] = Utilities::fermiDiracFunction(energy, sysParams->kT) * (D - y[0] * exp(energy/sysParams->kT) / sysParams->kT);
        f[1] = y[0];
        f[2] = energy * y[0];
        return GSL_SUCCESS;
    }

    static int differentialSystemLog(double energy, const double y[], double f[], void *params){
        BandEmitter* sysParams = (BandEmitter*) params; //cast the void pointer as SystemParams
        double D;
        //calculate transmission coefficient
        // TODO: this works only for effectiveMass = 1. FIX IT
        
        f[0] = Utilities::fermiDiracFunction(energy, sysParams->kT) * (D * exp(-y[0]) - exp(energy/sysParams->kT) / sysParams->kT);
        f[1] = exp(y[0] - y[1]);
        return GSL_SUCCESS;
    }

    static double normalEnergyDistribution(double energy, void* params){
        BandEmitter* sysParams = (BandEmitter*) params; //cast the void pointer as SystemParams
        double result = log(sysParams->calculateIntegrand(energy));// * Utilities::logFermiDiracFunction(energy, sysParams->kT);
        sysParams->savedEnergies.push_back(energy);
        sysParams->savedLogD.push_back(result);
        return result;
    }

    gsl_function integrationFunction = {&normalEnergyDistribution, this};

    gsl_integration_workspace* integrationWorkspace = NULL;

    static int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

public:

    void setParameters(double workFunction_ = 4.5, double kT_ = 0.025, double effectiveMass_ = 1., double bandDepth_ = 10.){
        if (workFunction != workFunction_ || bandDepth_ != bandDepth){
            workFunction = workFunction_;
            bandDepth = bandDepth_;
            transmissionSolver.setXlimits(workFunction + bandDepth + 2.);
        }
        
        xInitial = -bandDepth;

        effectiveMass = effectiveMass_;
        kT = kT_;
        xFinal = workFunction + 10. * kT;
        initialStep = (xFinal - xInitial) / stepsExpectedForInitialStep;
        setInitialValues({0., 0., 0.});
    }
    

    BandEmitter(TransmissionSolver& solver,
                double workFun = 4.5,
                double kT_ = .025,
                double effMass = 1.,
                double bandDepth_ = 7.,
                double rtol = 1.e-4,
                double atol = 1.e-12,
                int maxSteps = 4096,
                int minSteps = 64,
                int stepExpectedForInitialStep = 256
                )   :   ODESolver(vector<double>(3, 0.0), differentialSystem, 3, {0., 1.}, rtol, atol, 
                            gsl_odeiv2_step_rkck, maxSteps, minSteps, stepExpectedForInitialStep, NULL, this), 
                            transmissionSolver(solver), 
                            interpolator(solver, workFun, kT_, atol, rtol)
    {
        setParameters(workFun, kT_, effMass, bandDepth_);
        updateBarrier();
    }

    void updateBarrier(){
        transmissionSolver.setXlimits(workFunction + bandDepth + 1.);
        interpolator.initialize(-bandDepth, 10 * kT + workFunction, 8);
        interpolator.refineToTolerance();
    }

    int calculateCurrentDensityAndSpectra(double convergenceTolerance = 1.e-7){
        double x = xInitial;
        double dx = initialStep;
        int status;
        reinitialize();
        vector<double> previousSolution;
        for (size_t i = 0; i < maxAllowedSteps; i++){ //loop over blocks
            previousSolution = solutionVector;
            savedSolution.push_back(solutionVector);
            xSaved.push_back(x);
            dx = GSL_SIGN(dx) *  min(abs(maxStepSize), abs(dx));
            status = gsl_odeiv2_evolve_apply(evolver, controller, step, &sys, &x, xFinal, &dx, solutionVector.data());
            bool hasConverged = abs(previousSolution[1] - solutionVector[1]) / previousSolution[1]  < convergenceTolerance;
            if (x == xFinal || status != GSL_SUCCESS || hasConverged)    
                return status;         
        }
        return GSL_CONTINUE; 
    }

    double calcualteCurrentDensity(){
        double result, error;
        if (!integrationWorkspace)
            integrationWorkspace = gsl_integration_workspace_alloc(maxAllowedSteps);

        gsl_integration_qag(&integrationFunction, xInitial, xFinal, 0., relativeTolerance, maxAllowedSteps, GSL_INTEG_GAUSS41, integrationWorkspace, &result, &error);
        return result;
    }

    void writeSavedLogD(string filename = "energyTransmission.dat"){
            ofstream outFile(filename, ios::out);        
            for (size_t i = 0; i < savedEnergies.size(); i++)
                outFile << savedEnergies[i] << " " << savedLogD[i] << endl;
    }

    void writePlottingData(string filename = "bandEmitterPlotting.dat"){
        ofstream outFile(filename, ios::out);        
        outFile << " E D_calc D_interp error NED lFD(E) tolerance" << endl;
        for (double x = xInitial; x < xFinal; x += 0.001){
            double D = transmissionSolver.calculateTransmissionCoefficientForEnergy(x - workFunction);
            double lFD = Utilities::logFermiDiracFunction(x, kT);
            double err = interpolator.calculateError(x, log(D));
            double tol = interpolator.calculateTolerance(x, log(D));
            outFile << x << " " << D << " " << interpolator.evaluate(x) << " " << err << " " << lFD * D << " " << lFD << " " << tol <<  endl;
        }
        interpolator.writeSplineNodes();
    }

    double getCurrentDensity(){
        return solutionVector[1] * CONSTANTS.SommerfeldConstant;
    }
};  



#endif /* BANDEMITTER_H_ */
