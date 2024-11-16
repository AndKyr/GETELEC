#ifndef TRANSMISSIONSOLVER_H_
#define TRANSMISSIONSOLVER_H_




#include "ODESolver.h"
#include "TunnelingFunction.h"
#include "Utilities.h"

using namespace std;


class TransmissionSolver : public ODESolver{
private:
    TunnelingFunction* tunnelingFunction;
     
    double kappaInitial;
    double kappaFinal;

    static int tunnelingDifferentialSystem(double x, const double y[], double f[], void *params);

    static int tunnelingSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

    int numberOfCalls = 0;

public:   

    TransmissionSolver(
            TunnelingFunction* tunnelFunctionPtr, 
            double relativeTolerance = 1.e-5,
            double absoluteTolerance = 1.e-5,
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

    void setEnergy(double E){
        tunnelingFunction->setEnergy(E);
        updateKappaAtLimits();
    }

    void resetNumberOfCalls(){numberOfCalls = 0;}

    double transmissionCoefficient() const;

    double calculateTransmissionCoefficientForEnergy(double energy){
        setEnergy(energy);
        solveNoSave();
        numberOfCalls++;
        return transmissionCoefficient();
    }

    int getNumberOfCalls(){return numberOfCalls;}

    void printXLimits(){ cout << "xInitial = " << xInitial << " xFinal = " << xFinal << endl;}
};


class TransmissionInterpolator : public FunctionInterpolator{
private: 
    TransmissionSolver& solver;
    double kT;
    double workFunction;

    double calculateYforX(double energy) override{
        return log(solver.calculateTransmissionCoefficientForEnergy(-workFunction + energy));
    }

    double calculateError(double energy, double yCalculated) override{
        return Utilities::fermiDiracFunction(energy, kT) * abs(exp(yCalculated) - exp(gsl_spline_eval(spline, energy, accelerator)));
    }

    double calculateTolerance(double energy, double transmission) override {
        double maxEmissionEstimate = 0.;
        for (auto& v : samplingList){
            double emissionEstimate = exp(v.y) * Utilities::fermiDiracFunction(v.x, kT);
            if (emissionEstimate > maxEmissionEstimate)
                maxEmissionEstimate = emissionEstimate;
        }

        double emissionEstimate = exp(transmission) * Utilities::fermiDiracFunction(energy, kT);

        if (energy <= 0)
            return absoluteTolerance + (emissionEstimate + maxEmissionEstimate) * relativeTolerance;
        else
            return absoluteTolerance + emissionEstimate * relativeTolerance;
    }

public:

    /**The interpolator lives in the E_F = 0 convention */
    TransmissionInterpolator(TransmissionSolver& solver_, double workFunction_ = 4.5, double kT_ = .025, 
                                double aTol = 1.e-12, double rTol = 1.e-5) 
                            : FunctionInterpolator(aTol, rTol), solver(solver_),
                                workFunction(workFunction_), kT(kT_)
    {}
    

    void setParameters(double kT_, double W){
        kT = kT_;
        workFunction = W;
    }

    double evaluate(double x) override{
        return exp(gsl_spline_eval(spline, x, accelerator));
    }
};


#endif /* TRANSMISSIONSOLVER_H_ */
