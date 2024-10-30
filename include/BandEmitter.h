#ifndef BANDEMITTER_H_
#define BANDEMITTER_H_

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

#include <vector>
#include <typeinfo>
#include <fstream>
#include <string>

#include "transmissionCalculator.h"

using namespace std;



int differentialSystem(double x, const double y[], double f[], void *params);

class BandEmitter{
private:

    TunnelingFunctionBase* barrier;
    TransmissionCalculator& transmissionCalculator;

    int systemDimension = 2;
    double relativeTolerance = 1.e-3;
    double absoluteTolerance = 1.e-1;
    double initialStep;

    TunnelingFunctionBase* tunnelingFunction;

    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rk4;
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(stepType, systemDimension);
    gsl_odeiv2_control *controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    gsl_odeiv2_evolve *evolver = gsl_odeiv2_evolve_alloc(systemDimension);
    gsl_odeiv2_system sys;

    vector<double> solutionVector = vector<double>(systemDimension, 0.0);

public:
    

    void setTolerances(double absoluteTolerance = 1.e-5, double relativeTolerance = 1.e-5){
        relativeTolerance = relativeTolerance;
        absoluteTolerance = absoluteTolerance;
        controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    }


    BandEmitter(TunnelingFunctionBase* tunnelFunctionPtr, 
                            int systemDimension = 2, 
                            double relativeTolerance = 1.e-4,
                            double absoluteTolerance = 1.e-4,
                            const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd
                        );


    ~BandEmitter(){
        gsl_odeiv2_evolve_free(evolver);
        gsl_odeiv2_control_free(controller);
        gsl_odeiv2_step_free(step);
    }


    int solveDifferentialSystem();

};  



#endif /* BANDEMITTER_H_ */
