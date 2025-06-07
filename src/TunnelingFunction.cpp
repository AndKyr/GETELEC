#include "TunnelingFunction.h"

namespace getelec{

void ModifiedSNBarrier::setRandomParameters(){
    if (!generator)
        throw runtime_error("Random number generator is empty");
    field = Utilities::getUniformRandomDouble(1.e-5, 20., *generator);
    double curvature = Utilities::getUniformRandomDouble(1.e-5, 10., *generator);
    radius = 1. / curvature;
    gamma = Utilities::getUniformRandomDouble(1., 20., *generator);
}


void ModifiedSNBarrier::setField(double f){
    if (f < 1.e-5){
        cout << "Warning: Field is too small. Setting it to 1.e-5 V/nm" << endl;
        f = 1.e-5;
    }
    field = f;
}

void ModifiedSNBarrier::setRadius(double R){
    if (R < 0.1){
        cout << "Warning: Radius is too small. Setting it to 0.1 nm" << endl;
        R = 0.1;
    }
    radius = R;
}

void ModifiedSNBarrier::setGamma(double g){
    if (g < 1.0001){
        cout << "Warning: Gamma is too small. Setting it to 1.0001" << endl;
        g = 1.0001;
    }
    gamma = g;
}

double ModifiedSNBarrier::findRightXLimit(double maxPotentialDepth){
    double b = field * radius * (gamma - 1.) - gamma * maxPotentialDepth;
    double c = - maxPotentialDepth * radius * (gamma - 1.);
    double result = 0.5 * (-b + sqrt(b*b - 4*field*c)) / field;
    return max(result, 1.);
}


double ModifiedSNBarrierWithDftXC::XCPotential(double z){

    if (z <= xcFunctionParams.extensionStartPoint){
        cerr << "Serious warning: z is out of range of validity of XC function. Returning 0." << endl;
        return 0.;
    }
    double imagePotential = 0.;
    if (z > 0.) imagePotential = ModifiedSNBarrier::XCPotential(z);
    if (imagePotential > 20.) imagePotential = 20.;

    double dftXCPotential = 0.;
    if (z < xcFunctionParams.polynomialRange[0])
        dftXCPotential = xcFunctionParams.extensionPrefactor / (z - xcFunctionParams.extensionStartPoint);
    else if (z < xcFunctionParams.polynomialRange[1])
        dftXCPotential = gsl_poly_eval(xcFunctionParams.dftXcPolynomial.data(), xcFunctionParams.dftXcPolynomial.size(), z);

    double transitionFunction = .5 * gsl_sf_erfc((z - xcFunctionParams.transitionPoint) / xcFunctionParams.transisionWidth);
    return transitionFunction * dftXCPotential + (1. - transitionFunction) * imagePotential;
}


double ModifiedSNBarrierWithDftXC::XCPotentialDerivative(double z){
    if (z <= xcFunctionParams.extensionStartPoint){
        cerr << "Serious warning: z is out of range of validity of XC function. Returning 0." << endl;
        return 0.;
    }

    double imagePotential = 0., imagePotentialDerivative = 0.;
    if (z > 0.) {
        imagePotential = ModifiedSNBarrier::XCPotential(z);
        imagePotentialDerivative = ModifiedSNBarrier::XCPotentialDerivative(z);
    }

    if (imagePotential > 20.){
        imagePotential = 20.;
        imagePotentialDerivative = 0.;
    }

    double dftXCPotentialAndDerivative[2];
    if (z < xcFunctionParams.polynomialRange[0]){
        dftXCPotentialAndDerivative[0] = xcFunctionParams.extensionPrefactor / (z - xcFunctionParams.extensionStartPoint);
        dftXCPotentialAndDerivative[1] = - xcFunctionParams.extensionPrefactor / ((z - xcFunctionParams.extensionStartPoint) * (z - xcFunctionParams.extensionStartPoint));
    } else if (z < xcFunctionParams.polynomialRange[1]){
        gsl_poly_eval_derivs(xcFunctionParams.dftXcPolynomial.data(), xcFunctionParams.dftXcPolynomial.size(), z, dftXCPotentialAndDerivative, 2);
    }

    double transitionFunction = .5 * gsl_sf_erfc((z - xcFunctionParams.transitionPoint) / xcFunctionParams.transisionWidth);
    double transitionFunctionDerivative = - exp(- (z - xcFunctionParams.transitionPoint) * (z - xcFunctionParams.transitionPoint) / (xcFunctionParams.transisionWidth * xcFunctionParams.transisionWidth)) / (sqrt(M_PI) * xcFunctionParams.transisionWidth);
    return transitionFunction * dftXCPotentialAndDerivative[1] + dftXCPotentialAndDerivative[0] * transitionFunctionDerivative + 
            (1. - transitionFunction) * imagePotentialDerivative - imagePotential * transitionFunctionDerivative;
}

double ModifiedSNBarrierWithDftXC::findLeftXLimit(double maxPotentialDepth){
    double potentialAtXmin = XCPotential(xcFunctionParams.polynomialRange[0]);
    if (potentialAtXmin > maxPotentialDepth) 
        return xcFunctionParams.polynomialRange[0];
    else
        return xcFunctionParams.extensionStartPoint + xcFunctionParams.extensionPrefactor / maxPotentialDepth;
}

}