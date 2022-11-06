#include <math.h>

//precompiler defenitions to make life easier on referencing the dataArray
#define energy dataArray[0]
#define workFunction dataArray[1]
#define kT dataArray[2]
#define minimumEnergyDepth dataArray[3]
#define maximumEnergyDepth dataArray[4]
#define dGmin dataArray[5]
#define dGmax dataArray[6]
#define polynomial(i) (dataArray[dataArrayLength - i - 1])
#define polynomialLength dataArrayLength - 7

/**
 * Calculates the logarithm Fermi-Dirac distribution
*/
double logFermiDirac(double energyOverkT){
    if (energyOverkT > 20.)
        return exp(-energyOverkT);
    else if(energyOverkT < -20.)
        return -energyOverkT;
    else
        return log(1. + exp(-energyOverkT));
}

/**
 * Calculates the gamow factor by evaluating the polynomial
*/
double gamowFunction(int dataArrayLength, double *dataArray){
    double Gamow = polynomial(0);
    double energyDepth = workFunction - energy;
    
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
 * Evaluates the integrand of the current density (D(E) * lFD(E))
*/
double currentDensityPerNormalEnergy(int dataArrayLength, double *dataArray){

    double Gamow = gamowFunction(dataArrayLength, dataArray);
    
    double transmissionCoefficient;
    if (Gamow > 40.)
        transmissionCoefficient = exp(-Gamow);
    else
        transmissionCoefficient = 1. / (1. + exp(Gamow));

    return logFermiDirac(energy / kT) * transmissionCoefficient;
}

double nottinghamHeatInegrand(int dataArrayLength, double *dataArray){

    double Gamow = gamowFunction(dataArrayLength, dataArray);
    double transmissionCoefficient;
    extern double dilog_() ;

    if (Gamow > 40.)
        transmissionCoefficient = exp(-Gamow);
    else
        transmissionCoefficient = 1. / (1. + exp(Gamow));

    double exponent = -exp(-energy / kT);

    return transmissionCoefficient * (logFermiDirac(energy / kT) * energy - kT * dilog_(&exponent)) ;
}



