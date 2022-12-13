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
#define Ec dataArray[8]
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
 * Calculates the gamow factor at energy by evaluating the polynomial
*/
double gamowFunctionAtReducedEnergy(int dataArrayLength, double *dataArray){

    double aBar = 1. - effectiveMass;
    double energyDepth = workFunction - Ec - aBar * (energy - Ec) ;
    
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
 * Evaluates the integrand of the current density (D(E) * lFD(E))
*/
double currentDensityPerNormalEnergy(int dataArrayLength, double *dataArray){

    double Gamow = gamowFunction(dataArrayLength, dataArray);
        
    if (effectiveMass == 1.)
        return logFermiDirac(energy / kT) * transmissionCoefficientForGamow(Gamow);
    else{
        double GamowReduced = gamowFunctionAtReducedEnergy(dataArrayLength, dataArray);
        return logFermiDirac(energy / kT) * (transmissionCoefficientForGamow(Gamow) - (1. - effectiveMass) * transmissionCoefficientForGamow(GamowReduced));
    }
}

/**
 * @brief Calculates the Nottingham heat integrand
 * 
 * @param dataArrayLength the length of the input parameters
 * @param dataArray array of input parameters
 * @return double the integrand value
 */
double nottinghamHeatInegrand(int dataArrayLength, double *dataArray){

    double Gamow = gamowFunction(dataArrayLength, dataArray);
    double transmissionCoefficient;
    extern double dilog_();

    if (Gamow > 40.)
        transmissionCoefficient = exp(-Gamow);
    else
        transmissionCoefficient = 1. / (1. + exp(Gamow));

    double exponent = -exp(-energy / kT);

    return transmissionCoefficient * (logFermiDirac(energy / kT) * energy - kT * dilog_(&exponent)) ;
}



