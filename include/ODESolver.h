#ifndef ODESOLVER_H_
#define ODESOLVER_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

#include <cmath>
#include <vector>
#include <typeinfo>
#include <fstream>
#include <string>

using namespace std;

/**
 * @class ODESolver
 * @brief A class for solving ordinary differential equations (ODEs) using the GSL library.
 * 
 * This class provides an interface for solving ODE systems using the GNU Scientific Library (GSL).
 * It allows setting up the ODE system, initializing solver parameters, and performing the integration.
 */
class ODESolver {
protected:
    int systemDimension;                 /**< Dimension of the ODE system (number of equations). */
    double xInitial;                     /**< Initial value of the independent variable (e.g., time or spatial coordinate). */
    double xFinal;                       /**< Final value of the independent variable. */
    double maxStepSize;                  /**< Maximum allowed step size during integration. */
    double relativeTolerance = 1.e-4;    /**< Relative tolerance for the ODE solver. */
    double absoluteTolerance = 1.e-12;   /**< Absolute tolerance for the ODE solver. */
    int maxAllowedSteps = 4096;          /**< Maximum number of allowed steps in the solver. */
    int minAllowedSteps = 32;            /**< Minimum number of allowed steps in the solver. */
    int stepsExpectedForInitialStep = 64;/**< Expected number of steps for initial step size estimation. */
    double initialStep;                  /**< Initial step size for the solver. */
    
    const gsl_odeiv2_step_type *stepType = gsl_odeiv2_step_rk4; /**< Type of the GSL ODE stepper to use. */
    gsl_odeiv2_step *step;               /**< GSL ODE stepper object. */
    gsl_odeiv2_control *controller;      /**< GSL ODE controller for adaptive step sizing. */
    gsl_odeiv2_evolve *evolver;          /**< GSL ODE evolution object for performing the integration. */
    gsl_odeiv2_system sys;               /**< GSL ODE system structure containing function pointers and parameters. */
    vector<double> initialValues;        /**< Initial values of the dependent variables (state vector). */
    vector<vector<double>> savedSolution;/**< Storage for the solution at each step if saving is enabled. */
    vector<double> xSaved;               /**< Storage for the independent variable at each step if saving is enabled. */
    void* parameters;                    /**< Pointer to additional parameters required by the ODE system functions. */

    /**
     * @brief Function pointer to the ODE system function.
     * 
     * Should be of the form:
     * int func(double x, const double y[], double f[], void *params);
     * where y[] is the state vector at x, f[] is the derivative vector (dy/dx), and params is a pointer to any additional parameters.
     */
    int (*differentialSystem)(double, const double*, double*, void*);

    /**
     * @brief Function pointer to the Jacobian of the ODE system (optional).
     * 
     * Should be of the form:
     * int jac(double x, const double y[], double *dfdy, double dfdt[], void *params);
     * where dfdy is the Jacobian matrix, dfdt is the derivative of f with respect to x, and params is a pointer to any additional parameters.
     */
    int (*differentialSystemJacobian)(double, const double*, double*, double*, void*);

    vector<double> solutionVector = vector<double>(systemDimension, 0.0); /**< Current state vector during integration. */

public:

    /**
     * @brief Constructor for the ODESolver class.
     * 
     * Initializes the ODE solver with the given parameters.
     * 
     * @param initialValues Initial values of the dependent variables (state vector).
     * @param differentialSystem Pointer to the ODE system function.
     * @param systemDimension Dimension of the ODE system.
     * @param xLims Vector containing the initial and final values of the independent variable.
     * @param relativeTolerance Relative tolerance for the ODE solver.
     * @param absoluteTolerance Absolute tolerance for the ODE solver.
     * @param stepType Type of the GSL ODE stepper to use.
     * @param maxSteps Maximum number of allowed steps in the solver.
     * @param minSteps Minimum number of allowed steps in the solver.
     * @param stepExpectedForInitialStep Expected number of steps for initial step size estimation.
     * @param differentialSystemJacobian Pointer to the Jacobian function (optional).
     * @param params Pointer to additional parameters required by the ODE system functions.
     */
    ODESolver( 
                vector<double> initialValues,
                int (*differentialSystem)(double, const double*, double*, void*),
                int systemDimension = 3, 
                vector<double> xLims = {0., 1.},
                double relativeTolerance = 1.e-4,
                double absoluteTolerance = 1.e-4,
                const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd,
                int maxSteps = 4096,
                int minSteps = 32,
                int stepExpectedForInitialStep = 64,
                int (*differentialSystemJacobian)(double, const double*, double*, double*, void*) = NULL,
                void* params = NULL
            );

    /**
     * @brief Destructor for the ODESolver class.
     * 
     * Frees the GSL objects used in the solver.
     */
    ~ODESolver(){
        gsl_odeiv2_evolve_free(evolver);
        gsl_odeiv2_control_free(controller);
        gsl_odeiv2_step_free(step);
    }

    /**
     * @brief Sets the tolerances for the ODE solver.
     * 
     * @param aTol Absolute tolerance.
     * @param rTol Relative tolerance.
     */
    void setTolerances(double aTol = 1.e-5, double rTol = 1.e-5){
        relativeTolerance = rTol;
        absoluteTolerance = aTol;
        controller = gsl_odeiv2_control_y_new(absoluteTolerance, relativeTolerance);
    }

    /**
     * @brief Sets the initial values of the dependent variables.
     * 
     * @param initValues Vector of initial values.
     */
    void setInitialValues(vector<double> initValues){ initialValues = initValues; }

    /**
     * @brief Reinitializes the solver to start a new integration.
     * 
     * Resets the solution vector, GSL stepper, and clears saved solutions.
     */
    void reinitialize(){
        solutionVector = initialValues;
        gsl_odeiv2_step_reset(step);
        gsl_odeiv2_evolve_reset(evolver);
        savedSolution.clear();
        xSaved.clear();
    }

    /**
     * @brief Solves the ODE system.
     * 
     * @param saveSolution If true, saves the solution at each step.
     * @return GSL_SUCCESS on success, or an error code.
     */
    int solve(bool saveSolution = false);

    /**
     * @brief Solves the ODE system without saving the solution at each step.
     * 
     * @return GSL_SUCCESS on success, or an error code.
     */
    int solveNoSave();

    /**
     * @brief Writes the saved solution to a file.
     * 
     * @param filename Name of the output file.
     */
    void writeSolution(string filename = "odeSolution.dat");

};

#endif /* ODESOLVER_H_ */
