#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cmath>
#include <vector>
#include <limits>
#include <list>
#include <iostream>
#include <fstream>

#include <string>


#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace std;

static constexpr struct PhysicalConstants{
    double hbarSqrOver2m = 3.80998212e-2;
    double kConstant = 1./hbarSqrOver2m;
    double imageChargeConstant = .359991137;
    double BoltzmannConstant = 8.617333262e-5;
    double SommerfeldConstant = 1.618311e-4;
    double exponentLimit = - 0.5 * log(numeric_limits<double>::epsilon());

} CONSTANTS;

struct Utilities{
    static vector<double> linspace(double start, double end, int n);

    static double fermiDiracFunction(double energy, double kT);

    static double logFermiDiracFunction(double energy, double kT);
};


class FunctionInterpolator{
protected:
    // TransmissionSolver& solver;
    double absoluteTolerance = 1.e-12;
    double relativeTolerance = 1.e-4;

    gsl_interp_accel *accelerator = NULL;
    gsl_spline *spline = NULL;

    struct SplineElement{
        double x;
        double y;
        bool bisect;
        SplineElement(double xx, double yy, bool b) : x(xx), y(yy), bisect(b){}
    };

    list<SplineElement> samplingList;

    virtual double calculateYforX(double x){return exp(x);}

    virtual double calculateError(double x, double yCalculated){
        return abs(yCalculated - gsl_spline_eval(spline, x, accelerator));
    }

    virtual double calculateTolerance(double x, double yValue){
        return absoluteTolerance + relativeTolerance * abs(yValue);
    }

    int updateSpline();

    int refineSampling();

public:
    FunctionInterpolator(double aTol = 1.e-12, double rTol = 1.e-5
                        ) : 
        absoluteTolerance(aTol), relativeTolerance(rTol)
    {}

    ~FunctionInterpolator(){
        gsl_spline_free(spline);
        gsl_interp_accel_free(accelerator);
    }

    virtual double evaluate(double x){
        return gsl_spline_eval(spline, x, accelerator);
    }

    void initialize(double xInit, double xFinal, int numberOfInitialElements);


    void refineToTolerance(int maxRefiningSteps = 10){
        for (int i = 0; i < maxRefiningSteps; i++)
            if(refineSampling())
                updateSpline();
            else
                return;
    }

    void writeSplineNodes(string filename = "SplineNodes.dat"){
        ofstream outFile(filename, ios::out);        
        for (SplineElement& element : samplingList)
            outFile << element.x << " " << element.y << endl;
    }
};

#endif //UTILITIES_H_
