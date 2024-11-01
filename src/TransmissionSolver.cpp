#include "TransmissionSolver.h"


TransmissionSolver::TransmissionSolver(TunnelingFunctionBase* tunnelFunctionPtr, double rtol, double atol, const gsl_odeiv2_step_type* stepType,
                         int maxSteps, int stepExpectedForInitialStep, double maxPotentialDepth
                        ) : ODESolver(vector<double>(3, 0.0), tunnelingDifferentialSystem, 3, {2.00400712, 0.03599847},
                                        rtol, atol, stepType, maxSteps, stepExpectedForInitialStep, tunnelingSystemJacobian, tunnelFunctionPtr),
                            tunnelingFunction(tunnelFunctionPtr)
                            

{
    setXlimits(maxPotentialDepth);
    updateKappaAtLimits();
}