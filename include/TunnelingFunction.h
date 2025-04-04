#ifndef TUNNELINGFUNCTION_H_
#define TUNNELINGFUNCTION_H_

#include <cmath>
#include "Utilities.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_erf.h>
#include "ConfigGetelec.h"
#include <cassert>

namespace getelec{
using namespace std;


/**
 * @class TunnelingFunction
 * @brief A base class for tunneling-related calculations.
 * 
 * This class provides a framework for computing potential functions and their derivatives,
 * along with methods to calculate limits and tunneling parameters. A specific implementation 
 * of a barrier can inherit from this base class.
 */
class TunnelingFunction {
private:
    double energy = 0.; /**< Energy level for the tunneling function (eV). */
    bool doFindBarrierTop = false; /**< Flag to find the maximum kappa value. */
    double minKappaSquared = numeric_limits<double>::infinity();; /**< Maximum kappa value. */

protected:
    /** @brief A pointer to a random number generator to use for testing purposes */
    mt19937* generator = NULL;

public:
    /**
     * @brief Default constructor for TunnelingFunction.
     */
    TunnelingFunction() {}

    /**
     * @brief Constructor that initializes the energy level.
     * @param E Energy level (eV).
     */
    TunnelingFunction(double E) : energy(E) {}

    /**
     * @brief Virtual destructor for TunnelingFunction.
     */
    virtual ~TunnelingFunction() {}

    /**
     * @brief Computes the potential at a given position.
     * @param x Position coordinate.
     * @return Potential value at the specified position.
     * 
     * @note This method should be overridden in derived classes for specific potential models.
     */
    virtual double potentialFunction(double x) {
        cout << "This is the base class potential function (virtual). It should not be used directly. Just returns 0" << endl;
        return 0. * x; 
    }

    /**
     * @brief Computes the derivative of the potential function.
     * @param x Position coordinate.
     * @return Derivative of the potential at the specified position.
     * 
     * @note Override this method in derived classes for custom behavior.
     */
    virtual double potentialFunctionDerivative(double x) {
        cout << "This is the base class potential derivative function (virtual). It should not be used directly. Just returns 0" << endl;
        return 0. * x;
    }

    /**
     * @brief Finds the left boundary of the tunneling region.
     * @param maxPotentialDepth The maximum potential depth.
     * @return The left boundary of the potential well.
     */
    virtual double findLeftXLimit(double maxPotentialDepth) { 
        cout << "This is a base class function (virtual). It should not be used directly. Just returns 0" << endl;
        return 0.; 
    }

    /**
     * @brief Finds the right boundary of the tunneling region.
     * @param maxPotentialDepth The maximum potential depth.
     * @return The right boundary of the potential well.
     */
    virtual double findRightXLimit(double maxPotentialDepth) { 
        cout << "This is the base class potential function (virtual). It should not be used directly. Just returns 0" << endl;

        return 0.;
    }

    /**
     * @brief Computes the squared kappa value for tunneling calculations.
     * @param x Position coordinate.
     * @return The squared kappa value at the specified position.
     * @note If the flag doFindBarrierTop is set, the maximum kappa value is also calculated.
     */
    double kappaSquared(double x) {
        double kappaSquared = CONSTANTS.kConstant * (energy - this->potentialFunction(x));
        if (doFindBarrierTop && kappaSquared < minKappaSquared)
            minKappaSquared = kappaSquared;
        return kappaSquared;
    }

    /**
     * @brief sets the flag to find the maximum kappa value and resets the current value.
     * @param flag Flag value.
     */
    void setBarrierTopFinder(bool flag) { 
        doFindBarrierTop = flag; 
        if (flag)
            minKappaSquared = numeric_limits<double>::infinity();
    }
    
    double getBarrierTop() const { 
        assert(isfinite(minKappaSquared) && "Attempt to get the barrier top without finding it first");
        return energy - minKappaSquared / CONSTANTS.kConstant;
    }

    /**
     * @brief Computes the deriative of the squared kappa value for tunneling.
     * @param x Position coordinate.
     * @return The derivative value.
     */
    double kappaSquaredDerivative(double x){
        return -CONSTANTS.kConstant * this->potentialFunctionDerivative(x);
    }

    /**
     * @brief Sets the energy level for the tunneling function.
     * @param E Energy level (eV).
     */
    void setEnergy(double E) { energy = E; }

    /**
     * @brief Retrieves the energy level of the tunneling function.
     * @return The energy level (eV).
     */
    double getEnergy() const { return energy; }

    /**
     * @brief Sets the random number generator for the barrier (used to create random barrier for testing).
     * @param generator_ A pointer to the random number generator.
     */
    void setGenerator(mt19937* generator_){
        generator = generator_;
    }
};

/**
 * @class ModifiedSNBarrier
 * @brief A derived class representing the modified Schottky-Nordheim tunneling barrier. (see 10.1016/j.commatsci.2016.11.010 for details)
 */
class ModifiedSNBarrier : public TunnelingFunction {
protected:
    double radius = 1.e3; /**< Local radius of curvature of the emitter at the point of interest (nm). */
    double field = 5.; /**< Local electric field at the emission point (V/nm) */
    double gamma = 10.; /**< Gamma parameter controlling the far-away shape of the barrier */



    /**
     * @brief Calculates the image potential at a point z.
     * @param z coordinate - distance from electrical surface (nm)
     */
    virtual double XCPotential(double z){
        return CONSTANTS.imageChargeConstant /  (z + .5 * z * z / radius);
    }

    /**
     * @brief Calculates the derivative of the image potential at a point z.
     * @param z coordinate - distance from electrical surface (nm)
     */
    virtual double XCPotentialDerivative(double z){
        double denominator = (z + .5 * z * z / radius);
        return -CONSTANTS.imageChargeConstant * (1. + z / radius) / (denominator * denominator);
    }

    /**
     * @brief Calculates the electrostatic potential (component of the barrier) at a point z.
     * @param z coordinate - distance from electrical surface (nm)
     */
    double electrostaticPotential(double z){
        if (z < 0.) return 0.;
        return field * (radius * (gamma - 1.) * z + z*z) / (gamma * z + radius * (gamma - 1.));
    }

    /**
     * @brief Calculates the derivative of the electrostatic potential at a point z.
     * @param z coordinate - distance from electrical surface (nm)
     */
    double electrostaticPotentialDerivative(double z){
        if (z < 0.) return 0.;
        return field * ((gamma-1)*(gamma-1) * radius*radius + 2*(gamma-1)*radius*z + gamma*z*z) / 
        (((gamma - 1)*radius + gamma*z) * ((gamma - 1)*radius + gamma*z));
    }

public:
    /**
     * @brief Default constructor for ModifiedSNBarrier.
     */
    ModifiedSNBarrier() {}

    /**
     * @brief Constructor that initializes barrier parameters.
     * @param r Radius of the emitter (nm).
     * @param f Electric field strength (V/nm).
     * @param g Gamma parameter controlling the barrier shape.
     */
    ModifiedSNBarrier(double r, double f, double g) : radius(r), field(f), gamma(g) {}


    /**
     * @brief Setter for all Modified SN barrier parameters.
     * @param r Radius of the emitter (nm).
     * @param f Electric field strength (V/nm).
     * @param g Gamma parameter controlling the barrier shape.
     */
    void setBarrierParameters(double f, double r, double g){
        setField(f);
        setRadius(r);
        setGamma(g);
    } 

    /**
     * @brief sets all the parameter of the barrier to (meaningful) random values
     */
    void setRandomParameters();

    /**
     * @brief Setter for the field parameter.
     * @param f Electric field strength (V/nm).
     */
    void setField(double f);

    /**
     * @brief Setter for the radius parameter.
     * @param R Local radius of curvature (nm).
     */
    void setRadius(double R);

    /**
     * @brief Setter for the gamma parameter.
     * @param g Gamma
     */
    void setGamma(double g);

    /**
     * @brief Getter for the field.
     * @return Field (GV/m).
     */
    double getField() const { return field; }

    /**
     * @brief Getter for the radius.
     * @return Radius (nm).
     */
    double getRadius() const { return radius; }

    /**
     * @brief Getter for the gamma parameter.
     * @return Gamma parameter.
     */
    double getGamma() const { return gamma; }

    /**
     * @brief Computes the potential at a given position.
     * @param x Position coordinate.
     * @return Potential value at the specified position.
     */
    double potentialFunction(double z){
        return -this->XCPotential(z) - electrostaticPotential(z);
    }

    /**
     * @brief Computes the derivative of the potential function.
     * @param x Position coordinate.
     * @return Derivative of the potential at the specified position.
     */
    double potentialFunctionDerivative(double z){
        return -this->XCPotentialDerivative(z) - electrostaticPotentialDerivative(z);
    }

    /**
     * @brief Finds the left limit of numerical integration for calculating the barrier. 
     * Estimates the left position where the left of the barrier has a  potential depth of @param maxPotentialDepth.
     * 
     * @return x position where the potential depth is @param maxPotentialDepth.
     */
    virtual double findLeftXLimit(double maxPotentialDepth){
        return 1. / (maxPotentialDepth / CONSTANTS.imageChargeConstant + 0.5 / radius);
    }

    /**
     * @brief Finds the right limit of numerical integration for calculating the barrier. 
     * Estimates the right position where the left of the barrier has a  potential depth of @param maxPotentialDepth.
     * 
     * @return x position where the potential depth is @param maxPotentialDepth.
     */
    double findRightXLimit(double maxPotentialDepth);

};

class ModifiedSNBarrierWithDftXC : public ModifiedSNBarrier {

public:
    ModifiedSNBarrierWithDftXC(Config::XCFunctionParams xcParams = Config().xcFunctionParams) : xcFunctionParams(xcParams) {}

    /**
     * @brief Calculates the image potential at a point z.
     * @param z coordinate - distance from electrical surface (nm)
     */
    double XCPotential(double z) override;

    /**
     * @brief Calculates the derivative of the image potential at a point z.
     * @param z coordinate - distance from electrical surface (nm)
     */
    double XCPotentialDerivative(double z) override;

    /**
     * @brief Finds the left limit of numerical integration for calculating the barrier. 
     * Estimates the left position where the left of the barrier has a  potential depth of @param maxPotentialDepth.
     * 
     * @return x position where the potential depth is @param maxPotentialDepth.
     */
    double findLeftXLimit(double maxPotentialDepth) override;


private:
    Config::XCFunctionParams xcFunctionParams;

};

}

#endif /* TUNNELINGFUNCTION_H_ */
