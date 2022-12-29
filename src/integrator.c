#include <math.h>

//precompiler defenitions to make life easier on referencing the dataArray
#define energy dataArray[0]
#define workFunction dataArray[1]
#define kT dataArray[2]
#define minimumEnergyDepth dataArray[3]
#define maximumEnergyDepth dataArray[4]
#define dGmin dataArray[5]
#define dGmax dataArray[6]
#define effectiveMass dataArray[7]
#define energyBandLimit dataArray[8]
#define polynomial(i) (dataArray[dataArrayLength - i - 1])
#define polynomialLength dataArrayLength - 9


/**
 * Calculates the logarithm Fermi-Dirac distribution
*/
double logFermiDirac(double energyOverkT){

    if (energyOverkT > 20.)
        return exp(-energyOverkT);
    else if(energyOverkT < -20.)
        return -energyOverkT + exp(energyOverkT);
    else
        return log(1. + exp(-energyOverkT));
}


double gamowFunctionForEnergy(int dataArrayLength, double *dataArray, double energyDepth){
    
    
    double Gamow = polynomial(0);   
    // set energyDepth within the limits of validity of the polynomial
    if (energyDepth > maximumEnergyDepth)
        energyDepth = maximumEnergyDepth;
    else if (energyDepth < minimumEnergyDepth)
        energyDepth = minimumEnergyDepth;
    
    //evaluate the polynomial
    double power = energyDepth;
    for (int i = 1; i < polynomialLength; i++){
        Gamow += polynomial(i) * power;
        power *= energyDepth;
    }

    //extrapolate in case it is out of bounds
    if (energyDepth > maximumEnergyDepth)
        Gamow += dGmax * (energyDepth - maximumEnergyDepth);
    else if (energyDepth < minimumEnergyDepth)
        Gamow += dGmin * (energyDepth - minimumEnergyDepth);

    return Gamow;
}

/**
 * Calculates the gamow factor at energy by evaluating the polynomial
*/
double gamowFunction(int dataArrayLength, double *dataArray){
    
    double energyDepth = workFunction - energy;
    return gamowFunctionForEnergy(dataArrayLength, dataArray, energyDepth);

}


/**
 * @brief  Calculates transmission coefficient for a given Gamow factor
 * @param Gamow The Gamow factor
 * @return double transmission coefficient
 */
double transmissionCoefficientForGamow(double Gamow){
    if (Gamow > 40.)
        return exp(-Gamow);
    else
        return 1. / (1. + exp(Gamow));
}

/**
 * @brief Calculates the derivative of the inner integral g(E) in the current density integration. See Andreas' notes.
 * @param dataArrayLength the length of the input parameters
 * @param dataArray array of input parameters
 * @return the inner integrand derivative value
 */

double innerIntegralDerivative(int dataArrayLength, double *dataArray){
    double Gamow = gamowFunction(dataArrayLength, dataArray);
    double output = transmissionCoefficientForGamow(Gamow);
    if (effectiveMass != 1.){
        Gamow = gamowFunctionForEnergy(dataArrayLength, dataArray, workFunction - energyBandLimit - (1. - effectiveMass) * (energy - energyBandLimit));
        output -= (1. - effectiveMass) * transmissionCoefficientForGamow(Gamow);
    }
    return output;
}

/**
 * Evaluates the integrand of the current density (g'(E) * lFD(E))
*/
double currentDensityIntegrand(int dataArrayLength, double *dataArray){
    return logFermiDirac(energy / kT) * innerIntegralDerivative(dataArrayLength, dataArray);
}

/**
 * @brief Calculates the Nottingham heat integrand
 * 
 * @param dataArrayLength the length of the input parameters
 * @param dataArray array of input parameters
 * @return the integrand value
 */
double nottinghamHeatIntegrand(int dataArrayLength, double *dataArray){

    extern double dilog_();
    double exponential = -exp(-energy / kT);

    return (logFermiDirac(energy / kT) * energy - kT * dilog_(&exponential)) * 
        innerIntegralDerivative(dataArrayLength, dataArray);
}



