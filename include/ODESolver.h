#ifndef ODESOLVER_H_
#define ODESOLVER_H_

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

#include <vector>
#include <typeinfo>
#include <fstream>
#include <string>


using namespace std;

class ODESolver{
protected:
    int systemDimension;
    double xInitial;
    double xFinal;
    double relativeTolerance;
    double absoluteTolerance;
    int maxAllowedSteps;
    int stepsExpectedForInitialStep = 64;
    double initialStep;
    
    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rk4;
    gsl_odeiv2_step *step;
    gsl_odeiv2_control *controller;
    gsl_odeiv2_evolve *evolver;
    gsl_odeiv2_system sys;
    vector<double> initialValues;
    vector<vector<double>> savedSolution;
    vector<double>xSaved;
    void* parameters;

    int (*differentialSystem)(double, const double*, double* , void*);
    int (*differentialSystemJacobian)(double, const double*, double*, double*, void*);

    vector<double> solutionVector = vector<double>(systemDimension, 0.0);

public:

    ODESolver( 
                vector<double>initialValues,
                int (*differentialSystem)(double, const double*, double*, void*),
                int systemDimension = 3, 
                vector<double> xLims = {0., 1.},
                double relativeTolerance = 1.e-4,
                double absoluteTolerance = 1.e-4,
                const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd,
                int maxSteps = 4096,
                int stepExpectedForInitialStep = 64,
                int (*differentialSystemJacobian)(double, const double*, double*, double*, void*) = NULL,
                void* params = NULL
            );

    ~ODESolver(){
        gsl_odeiv2_evolve_free(evolver);
        gsl_odeiv2_control_free(controller);
        gsl_odeiv2_step_free(step);
    }

    void setTolerances(double absoluteTolerance = 1.e-5, double relativeTolerance = 1.e-5){
        relativeTolerance = relativeTolerance;
        absoluteTolerance = absoluteTolerance;
        controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    }

    void setInitialValues(vector<double> initValues){initialValues = initValues;}

    void reinitialize(){
        solutionVector = initialValues;
        gsl_odeiv2_step_reset(step);
        gsl_odeiv2_evolve_reset(evolver);
        savedSolution.clear();
        xSaved.clear();
    }

    int solve(bool saveSolution = false);

    int solveNoSave();

    void writeSolution(string filename = "odeSolution.dat");

};


#endif /* ODESOLVER_H_ */
