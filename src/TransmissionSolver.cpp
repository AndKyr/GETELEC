#include "TransmissionSolver.h"
#include <cassert>
#include <sstream>
#include <gsl/gsl_complex_math.h>

namespace getelec{

int TransmissionSolver::tunnelingDifferentialSystem(double x, const double y[], double f[], void *params){
        TunnelingFunction* barrier = (TunnelingFunction*) params;
        //y[0] = Re{s'}, y[1] = Im{s'}, y[2] = Re{s}, y[3] = Im{s} (no need to follow Im{s} because it is just a phase factor)
        f[0] =  -barrier->kappaSquared(x) - y[0]*y[0] + y[1]*y[1]; //Re{s''} = -k^2 - Re{s'}^2 + Im{s'}^2
        f[1] = - 2. * y[0] * y[1]; //Im{s''} = -2 Re{s'} Im{s'}
        f[2] = y[0];//Re{s'} = Re{s'}
        // f[3] = y[1];//Im{s'} = Im{s'} . This is not necessary because Im{s} just gives a phase factor that doesn't matter in the end
        return GSL_SUCCESS;
}

int TransmissionSolver::tunnelingSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params){
        TunnelingFunction* barrier = (TunnelingFunction*) params;
        gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
        gsl_matrix *matrix = &dfdy_mat.matrix;

        gsl_matrix_set(matrix, 0, 0, -2*y[0]); gsl_matrix_set(matrix, 0, 1, 2*y[1]);  gsl_matrix_set(matrix, 0, 2, 0.);
        gsl_matrix_set(matrix, 1, 0, -2*y[1]); gsl_matrix_set(matrix, 1, 1, -2*y[0]); gsl_matrix_set(matrix, 1, 2, 0.);
        gsl_matrix_set(matrix, 2, 0, 1.);      gsl_matrix_set(matrix, 2, 1, 0.);      gsl_matrix_set(matrix, 2, 2, 0.);

        dfdt[0] = CONSTANTS.kConstant * barrier->potentialFunctionDerivative(x);
        dfdt[1] = 0.0; dfdt[2] = 0.0;
        return GSL_SUCCESS;
}

void TransmissionSolver::updateKappaInitial(){
    double kappaSquaredInitial = barrier->kappaSquared(xInitial);

    if (kappaSquaredInitial <= 0.) 
        throw std::runtime_error("The tunneling energy is lower than the edge potential values. The integration interval must extend beyond the classically forbidden region.");
    
    kappaInitial = sqrt(kappaSquaredInitial);
}

double TransmissionSolver::calculateKappaFinal() const{
    double kappaSquaredFinal = barrier->kappaSquared(xFinal);

    if (kappaSquaredFinal <= 0.) 
        throw std::runtime_error("The tunneling energy is lower than the left edge potential values. The integration interval must extend beyond the classically forbidden region.");
    
    return sqrt(kappaSquaredFinal);
}

double TransmissionSolver::calculateTransmissionProbability(double energy, double waveVector){
    numberOfCalls++;

    setEnergy(energy);

    //the tunneling region is really big, there is no point to calculate
    if (xInitial > 11.)
        return 1.e-50;

    initialValues = {0, 1, 0};
    solveNoSave();

    double k = waveVector;
    if (k <= 0.)
        k = calculateKappaFinal();
    
    double out = transmissionProbabilityforWaveVector(k);
    
    // Assert with inline error message construction
    assert(isfinite(out) && out >= 0. && ([](double val) {
        stringstream ss;
        ss << "transmission coefficient appears to not be finite. Current value: " << val;
        return ss.str().c_str();
    }(out)));

    return max(out, 1.e-50);
}

gsl_complex TransmissionSolver::transmissionCoefficientForWaveVector(double k) const{
    gsl_complex incidentWaveCoeff;
    double kappaInitialDerivaitve = 0.;// = 0.5 * barrier->kappaSquaredDerivative(xInitial) / kappaInitial;
    // GSL_SET_COMPLEX(&incidentWaveCoeff, 
    //     0.5*(k + solutionVector[1] * kappaInitial) / k / sqrt(kappaInitial), 
    //     0.25*(solutionVector[1] * kappaInitialDerivaitve - 2.* solutionVector[0]*kappaInitial ) / k / pow(kappaInitial,1.5)
    // );

    GSL_SET_COMPLEX(&incidentWaveCoeff, k + solutionVector[1] * kappaInitial, -solutionVector[0]);
    
    return gsl_complex_mul_real(gsl_complex_inverse(incidentWaveCoeff), 2. * exp(-solutionVector[2]) * k * sqrt(kappaInitial));
}

double TransmissionSolver::transmissionProbabilityforWaveVector(double k) const{
    // gsl_complex transmissionCoeff = transmissionCoefficientForWaveVector(k);
    // return gsl_complex_abs2(transmissionCoeff) / k;
    double factor = k + solutionVector[1] * kappaInitial;
    factor = solutionVector[0] * solutionVector[0] + factor * factor;
    return 4 * exp(-2*solutionVector[2]) * k * kappaInitial / factor;
}

void TransmissionSolver::calculateFundamentalMatrix(){
    initialValues = {0, 1, 0};
    solveNoSave();
    double prefactor = exp(solutionVector[2]);
    fundamentalMatrix[0] = prefactor;   
    fundamentalMatrix[1] = 0.;
    fundamentalMatrix[2] = prefactor * solutionVector[0];
    fundamentalMatrix[3] = prefactor * solutionVector[1];
}

} // namespace getelec