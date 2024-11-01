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


#include "TransmissionSolver.h"

using namespace std;



class BandEmitter : public ODESolver{
private:


    struct SystemParameters{
        double workFunction = 4.5;
        double kT = 0.025;
        double effectiveMass = 1.;
        TransmissionSolver* transmissionCalculator;
    } systemParams;


    double bandDepth = 10.;

    static int differentialSystem(double energy, const double y[], double f[], void *params){
        SystemParameters* sysParams = (SystemParameters*) params;

        sysParams->transmissionCalculator->setEnergy(-sysParams->workFunction + energy);
        sysParams->transmissionCalculator->solveNoSave();
        double D = sysParams->transmissionCalculator->transmissionCoefficient();

        f[0] = Utilities::fermiDiracFunction(energy, sysParams->kT) * (D - y[0] * exp(energy/sysParams->kT) / sysParams->kT);
        f[1] = y[0];
        return GSL_SUCCESS;

    }

    static int differentialSystemJacobian(double x, const double y[], double *dfdy, double dfdt[], void *params);

public:
    



    BandEmitter(TunnelingFunctionBase* tunnelFunctionPtr, 
                            int systemDimension = 2, 
                            double relativeTolerance = 1.e-4,
                            double absoluteTolerance = 1.e-4,
                            const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd
                        );

};  



#endif /* BANDEMITTER_H_ */
