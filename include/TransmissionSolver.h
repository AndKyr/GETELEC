#ifndef TRANSMISSIONSOLVER_H_
#define TRANSMISSIONSOLVER_H_



#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <list>

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


class FunctionInterpolator{
private:
    // TransmissionSolver& solver;
    double absoluteTolerance = 1.e-12;
    double relativeTolerance = 1.e-5;


    gsl_interp_accel *accelerator;
    gsl_spline *spline;

    struct SplineElement{
        double x;
        double y;
        bool bisect;

        SplineElement(double xx, double yy, bool b) : x(xx), y(yy), bisect(b){}
    };

    list<SplineElement> samplingList;

    virtual double calculateYforX(double x){return x;}

    int updateSpline(){
        gsl_spline_free(spline);
        gsl_interp_accel_free(accelerator);

        accelerator = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_cspline, samplingList.size());

        vector<double> x(samplingList.size());
        vector<double> y(samplingList.size());

        int i = 0;
        for (auto& val : samplingList){
            x[i] = val.x;
            y[i++] = val.y;
        }

        return gsl_spline_init(spline, x.data(), y.data(), x.size());
    }

    int initialize(double xInit, double xFinal, int numberOfInitialElements){
        samplingList.clear();
        gsl_spline_free(spline);
        gsl_interp_accel_free(accelerator);

        vector<double> xPoints = Utilities::linspace(xInit, xFinal, numberOfInitialElements);

        for (const auto& x : xPoints){
            samplingList.push_back(SplineElement(x, calculateYforX(x), true));
        }

        samplingList.front().bisect = false;
    }   

    int refineSampling(double atol, double rtol){
        for(auto it = samplingList.begin(); it != samplingList.end(); it++){
            if (it->bisect){
                double xNew = .5*(it->x + prev(it)->x);
                double yNew = calculateYforX(xNew);
                double error = abs(yNew - gsl_spline_eval(spline, xNew, accelerator));
                if (error < atol + rtol * abs(yNew)){
                    samplingList.emplace(it, SplineElement(xNew, yNew, false));
                    it->bisect = false;
                } else
                    samplingList.emplace(it, SplineElement(xNew, yNew, true));
            }
        }
        return 0;
    }

public:
    FunctionInterpolator();
    FunctionInterpolator(double xa = 0, double xb = 1, int N = 8, 
                        double aTol = 1.e-12, double rTol = 1.e-5, 
                        ) : 
        absoluteTolerance(aTol), relativeTolerance(rTol)
    {
        initialize(xa, xb, N);
    }

    ~FunctionInterpolator(){
        gsl_spline_free(spline);
        gsl_interp_accel_free(accelerator);
    }
}

#endif /* TRANSMISSIONSOLVER_H_ */
