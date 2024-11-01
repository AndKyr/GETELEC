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



class BandEmitter : public ODESolver{
private:

    TunnelingFunctionBase* barrier;
    TransmissionSolver& transmissionCalculator;

public:
    



    BandEmitter(TunnelingFunctionBase* tunnelFunctionPtr, 
                            int systemDimension = 2, 
                            double relativeTolerance = 1.e-4,
                            double absoluteTolerance = 1.e-4,
                            const gsl_odeiv2_step_type* stepType = gsl_odeiv2_step_rk8pd
                        );


    int solveDifferentialSystem();

};  



#endif /* BANDEMITTER_H_ */
