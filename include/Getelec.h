#ifndef GETELEC_H
#define GETELEC_H

#include <vector>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include "BandEmitter.h"
#include "TunnelingFunction.h"
#include "ConfigGetelec.h"

namespace getelec{

using namespace std;


enum class CalculationFlags : unsigned int {
    None            = 0,
    CurrentDensity  = 1 << 0, // 1
    NottinghamHeat  = 1 << 1, // 2
    TotalEnergyDistribution         = 1 << 2, // 4
    NormalEnergyDistribution        = 1 << 3, // 8
    ParallelEnergyDistribution      = 1 << 4 //16
};

inline bool operator&(CalculationFlags a, CalculationFlags b) {
    using T = std::underlying_type_t<CalculationFlags>;
    return (static_cast<T>(a) & static_cast<T>(b)) != 0;
}

inline CalculationFlags operator|(CalculationFlags a, CalculationFlags b) {
    using T = std::underlying_type_t<CalculationFlags>;
    return static_cast<CalculationFlags>(static_cast<T>(a) | static_cast<T>(b));
}

inline CalculationFlags& operator|=(CalculationFlags& a, CalculationFlags b) {
    a = a | b;
    return a;
}


class Getelec {
public:
    /**
     * @brief Construct a new Getelec object
     */
    Getelec(string configFileName = "GetelecConfig.txt", string barrierType = "modifiedSN") :
        config(Config(configFileName)),
        threadLocalBarrier([this, barrierType]() -> unique_ptr<ModifiedSNBarrier> {
            if (!barrierType.starts_with("dftXC")) {
                return make_unique<ModifiedSNBarrier>();
            } else {
                config.xcFunctionParams.name = barrierType;
                int nReadParams = config.readParamGroup(&config.xcFunctionParams);
                if (nReadParams < config.xcFunctionParams.keyMap.size()) {
                    cerr << "WARNING: No parameters found for the DFT XC function named " << barrierType 
                            << " ; Using default values, which are for W(110)." << std::endl;
                    return make_unique<ModifiedSNBarrierWithDftXC>();
                } else {
                    return make_unique<ModifiedSNBarrierWithDftXC>(config.xcFunctionParams);
                }
            }
        }), 
        threadLocalSolver([this] { return TransmissionSolver(threadLocalBarrier.local().get(), config.transmissionSolverParams, 10., 1); }),
        threadLocalEmitter([this] { return BandEmitter(threadLocalSolver.local(), config.bandEmitterParams); })
    {}


    ~Getelec(){}

    /**
     * @brief Set the field parameters of the calculation at a single value
     * @param field_ The electric field in V/nm.
     */
    void setField(double field_ = 5.) { 
        fieldsVector.resize(1);
        fieldsVector[0] = field_;
     }

    /**
     * @brief Set the field parameters of the calculation at multiple values
     * @param fieldsVector_ The electric field in V/nm, multiple values to iterate over.
     */
    void setField(std::vector<double>& fieldsVector_) { fieldsVector = vector<double>(fieldsVector_.begin(), fieldsVector_.end()); }

    /**
     * @brief Set the field parameters of the calculation at multiple values
     * @param fieldsVector_ The electric field in V/nm, multiple values to iterate over.
     * @note This method uses move semantics and is more efficient than the other setField method, but it invalidates the input vector.
     */
    void setFieldMove(std::vector<double>&& fieldsVector_) { fieldsVector = std::move(fieldsVector_); }

    /**
     * @brief Set the field parameters of the calculation at multiple values with a pointer and size input (C-style)
     * @param fieldsArray_ The electric field in V/nm, multiple values to iterate over.
     * @param size The size of the array
     */
    void setField(const double* fieldsArray_, size_t size) { fieldsVector = std::vector<double>(fieldsArray_, fieldsArray_ + size); }

    /**
     * @brief Set the radius of the emitter at a single value
     * @param radius_ The radius of the emitter in nm.
     */
    void setRadius(double radius_ = 1.e4) {
        radiiVector.resize(1);
        radiiVector[0] = radius_;
    }

    /**
     * @brief Set the radius of the emitter at multiple values
     * @param radiiVector_ The radius of the emitter in nm, multiple values to iterate over.
     */
    void setRadius(std::vector<double>& radiiVector_) { radiiVector = vector<double>(radiiVector_.begin(), radiiVector_.end()); }

    /**
     * @brief Set the radius of the emitter at multiple values
     * @param radiiVector_ The radius of the emitter in nm, multiple values to iterate over.
     * @note This method uses move semantics and is more efficient than the other setRadius method, but it invalidates the input vector.
     */
    void setRadiusMove(std::vector<double>&& radiiVector_) { radiiVector = std::move(radiiVector_); }

    /**
     * @brief Set the radius of the emitter at multiple values with a pointer and size input (C-style)
     * @param radiiArray_ The radius of the emitter in nm, multiple values to iterate over.
     * @param size The size of the array
     */
    void setRadius(const double* radiiArray_, size_t size) { radiiVector = std::vector<double>(radiiArray_, radiiArray_ + size); }

    /**
     * @brief Set the gamma parameter of the general barrier model at a single value
     * @param gamma_ The gamma parameter of the general barrier model.
     */
    void setGamma(double gamma_ = 10.){
        gammasVector.resize(1);
        gammasVector[0] = gamma_;
    }

    /**
     * @brief Set the gamma parameter of the general barrier model at multiple values
     * @param gammasVector_ The gamma parameter of the general barrier model, multiple values to iterate over.
     */
    void setGamma(std::vector<double>& gammasVector_) { gammasVector = vector<double>(gammasVector_.begin(), gammasVector_.end()); }

    /**
     * @brief Set the gamma parameter of the general barrier model at multiple values
     * @param gammasVector_ The gamma parameter of the general barrier model, multiple values to iterate over.
     * @note This method uses move semantics and is more efficient than the other setGamma method, but it invalidates the input vector.
     */
    void setGammaMove(std::vector<double>&& gammasVector_) { gammasVector = std::move(gammasVector_); }

    /**
     * @brief Set the gamma parameter of the general barrier model at multiple values with a pointer and size input (C-style)
     * @param gammasArray_ The gamma parameter of the general barrier model, multiple values to iterate over.
     * @param size The size of the array
     */
    void setGamma(const double* gammasArray_, size_t size) { gammasVector = std::vector<double>(gammasArray_, gammasArray_ + size); }

    /**
     * @brief Set the temperature at a single value
     * @param kT_ Temperature in eV.
     */
    void setkT(double kT_ = .025){
        kTVector.resize(1);
        kTVector[0] = kT_;
    }

    /**
     * @brief Set the temperature at multiple values
     * @param kTVector_ Temperature in eV, multiple values to iterate over.
     */
    void setkT(std::vector<double>& kTVector_) { kTVector = vector<double>(kTVector_.begin(), kTVector_.end()); }

    /**
     * @brief Set the temperature at multiple values
     * @param kTVector_ Temperature in eV, multiple values to iterate over.
     * @note This method uses move semantics and is more efficient than the other setkT method, but it invalidates the input vector.
     */
    void setkTMove(std::vector<double>&& kTVector_) { kTVector = std::move(kTVector_); }

    /**
     * @brief Set the temperature at multiple values with a pointer and size input (C-style)
     * @param kTArray_ Temperature in eV, multiple values to iterate over.
     * @param size The size of the array
     */
    void setkT(const double* kTArray_, size_t size) { kTVector = std::vector<double>(kTArray_, kTArray_ + size); }

    /**
     * @brief Set the work function at a single value
     * @param workFunction_ Work function in eV.
     */
    void setWorkFunction(double workFunction_ = 4.5){
        workFunctionVector.resize(1);
        workFunctionVector[0] = workFunction_;
    }
    
    /**
     * @brief Set the work function at multiple values
     * @param workFunctionVector_ Work function in eV, multiple values to iterate over.
     */
    void setWorkFunction(std::vector<double>& workFunctionVector_) { workFunctionVector = vector<double>(workFunctionVector_.begin(), workFunctionVector_.end()); }

    /**
     * @brief Set the work function at multiple values
     * @param workFunctionVector_ Work function in eV, multiple values to iterate over.
     * @note This method uses move semantics and is more efficient than the other setWorkFunction method, but it invalidates the input vector.
     */
    void setWorkFunctionMove(std::vector<double>&& workFunctionVector_) { workFunctionVector = std::move(workFunctionVector_); }

    /**
     * @brief Set the work function at multiple values with a pointer and size input (C-style)
     * @param workFunctionArray_ Work function in eV, multiple values to iterate over.
     * @param size The size of the array
     */
    void setWorkFunction(const double* workFunctionArray_, size_t size) { workFunctionVector = std::vector<double>(workFunctionArray_, workFunctionArray_ + size); }

    /**
     * @brief Set the band depth at a single value
     * @param bandDepth_ Depth of the electronic band in eV.
     */
    void setBandDepth(double bandDepth_ = 10.) {
        bandDepthVector.resize(1);
        bandDepthVector[0] = bandDepth_;
    }
    
    /**
     * @brief Set the band depth at multiple values
     * @param bandDepthVector_ Depth of the electronic band in eV, multiple values to iterate over.
     */
    void setBandDepth(std::vector<double>& bandDepthVector_) { bandDepthVector = vector<double>(bandDepthVector_.begin(), bandDepthVector_.end()); } 

    /**
     * @brief Set the band depth at multiple values
     * @param bandDepthVector_ Depth of the electronic band in eV, multiple values to iterate over.
     * @note This method uses move semantics and is more efficient than the other setBandDepth method, but it invalidates the input vector.
     */
    void setBandDepthMove(std::vector<double>&& bandDepthVector_) { bandDepthVector = std::move(bandDepthVector_); }
    /**
     * @brief Set the band depth at multiple values with a pointer and size input (C-style)
     * @param bandDepthArray_ Depth of the electronic band in eV, multiple values to iterate over.
     * @param size The size of the array
     */
    void setBandDepth(const double* bandDepthArray_, size_t size) { bandDepthVector = std::vector<double>(bandDepthArray_, bandDepthArray_ + size); }


    /**
     * @brief Set the effective mass at a single value
     * @param effectiveMass_ Effective mass of the electron.
     */
    void setEffectiveMass(double effectiveMass_ = 1.){
        effectiveMassVector.resize(1);
        effectiveMassVector[0] = effectiveMass_;
    }
    /**
     * @brief Set the effective mass at multiple values
     * @param effectiveMassVector_ Effective mass of the electron, multiple values to iterate over.
     */
    void setEffectiveMass(const vector<double>& effectiveMassVector_) { effectiveMassVector = vector<double>(effectiveMassVector_.begin(), effectiveMassVector_.end()); } 
    
    /**
     * @brief Set the effective mass at multiple values
     * @param effectiveMassVector_ Effective mass of the electron, multiple values to iterate over.
     * @note This method uses move semantics and is more efficient than the other setEffectiveMass method, but it invalidates the input vector.
     */
    void setEffectiveMassMove(std::vector<double>&& effectiveMassVector_) { effectiveMassVector = std::move(effectiveMassVector_); }

    /**
     * @brief Set the effective mass at multiple values with a pointer and size input (C-style)
     * @param effectiveMassArray_ Effective mass of the electron, multiple values to iterate over.
     * @param size The size of the array
     */
    void setEffectiveMass(const double* effectiveMassArray_, size_t size) { effectiveMassVector = std::vector<double>(effectiveMassArray_, effectiveMassArray_ + size); }

    /**
     * @brief Calculate for the i-th element of the array of inputs
     * @param i The index of the element to calculate
     * @param calculateSpectra If true, calculate the spectra
     */
    void runIteration(size_t i = 0, CalculationFlags flags = CalculationFlags::CurrentDensity);

    /**
     * @brief Run the calculation
     * @param calculateSpectra If true, calculate the spectra
     * @return The number of iterations
     */
    size_t run(CalculationFlags flags = CalculationFlags::CurrentDensity);


    /**
     * @brief Calculate the transmission coefficient for a specific energy
     * @param energy The energy level (eV)
     * @param paramsIndex The index of the parameter vector space to use (optional)
     * @return The transmission coefficient
     * @note This method is relevant for a single or a few calculations of the transmission coefficient. Don't use it for multiple calculations on the same barrier as it resets the barrier which might be slow. Use calculateTransmissionCoefficientForEnergies instead
     */
    double calculateTransmissionCoefficientForEnergy(double energy, size_t paramsIndex = 0);


    /**
     * @brief Calculate the transmission coefficient for multiple energies
     * @param energies The energy levels (eV)
     * @return The transmission coefficients
     * @note This method is relevant for multiple calculations of the transmission coefficient. It is faster than calculateTransmissionCoefficientForEnergy for multiple calculations on the same barrier. However, if you are iterating over many many energies, it might be better to use calculateTransmissionCoefficientForManyEnergies, which prepares the interpolator and then just interpolates.
     */
    vector<double> calculateTransmissionCoefficientForEnergies(const vector<double>& energies, size_t paramsIndex = 0);

    /**
     * @brief Calculate the transmission coefficient for multiple energies
     * @param energies The energy levels (eV)
     * @return The transmission coefficients
     * @note This method is relevant for multiple calculations of the transmission coefficient. It is faster than calculateTransmissionCoefficientForEnergy for multiple calculations on the same barrier. However, if you are iterating over many many energies, it might be better to use calculateTransmissionCoefficientForManyEnergies, which prepares the interpolator and then just interpolates.
     */
    vector<double> calculateTransmissionCoefficientForManyEnergies(const vector<double>& energies, size_t paramsIndex = 0);

    /**
     * @brief Get the current density at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The current density in A/nm^2
     */
    double getCurrentDensity(unsigned i = 0) const { return currentDensityVector[i]; }
    
    /**
     * @brief Get the Nottingham heat at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The Nottingham heat in W/nm^2
     */
    double getNottinghamHeat(unsigned i = 0) const { return nottinghamHeatVector[i]; }
    
    /**
     * @brief Get the spectra at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The spectra in A/nm^2/eV
     */
    const tuple<vector<double>, vector<double>, vector<double>>& getSpectra(unsigned i = 0) const { return spectra[i]; }
    
    /**
     * @brief Get all the spectra
     * @param i The index of the element to get
     * @return The spectra in A/nm^2/eV
     */
    vector<tuple<vector<double>, vector<double>, vector<double>>> getAllSpectra() const { 
        return spectra; 
    }
    
    /**
     * @brief Get all the resulting current densities as a vector
     * @return The vector
     * @note This is a copy of the internal vector, so it is safe to modify it.
     */
    vector<double> getCurrentDensities() const { return vector<double>(currentDensityVector.begin(), currentDensityVector.end()); }
    
    /**
     * @brief Get all the resulting current densities as a vector (used for C and python interface)
     * @param size The size of the vector to be output
     * @return Pointer to the first element of the vector
     */
    const double* getCurrentDensities(size_t* size) const { 
        *size = currentDensityVector.size();
        return currentDensityVector.data(); 
    }

    /**
     * @brief Get all the resulting Nottingham heats as a vector
     * @return The vector
     * @note This is a copy of the internal vector, so it is safe to modify it.
     */
    vector<double> getNottinghamHeats() const { return vector<double>(nottinghamHeatVector.begin(), nottinghamHeatVector.end()); }

    /**
     * @brief Get all the resulting Nottingham heats as a vector (used for C and python interface)
     * @param size The size of the vector to be output
     * @return Pointer to the first element of the vector
     */
    const double* getNottinghamHeats(size_t* size) const { 
        *size = nottinghamHeatVector.size();
        return nottinghamHeatVector.data(); 
    }

    /**
     * @brief Get the spectra abssicae (energies in eV) of the i-th iteration in parameter space
     * @param i The index of the iteration
     * @param length The length of the spectra array to be output
     * @return Pointer to the first element of the spectra array
     */
    const double* getSpectraEnergies(size_t i, size_t* length) const {
        auto& spectra_i = get<0>(spectra[i]);
        *length = spectra_i.size();
        return spectra_i.data(); 
    }

    /**
     * @brief Get the spectra values (ordinates) (values in A/nm^2/eV) of the i-th iteration in parameter space
     * @param i The index of the iteration
     * @param length The length of the spectra array to be output
     * @return Pointer to the first element of the spectra array
     */
    const double* getSpectraValues(size_t i, size_t* length) const { 
        *length = get<1>(spectra[i]).size();
        return get<1>(spectra[i]).data(); 
    }

    /**
     * @brief Get the spectra derivatives (values in A/nm^2/eV^2) of the i-th iteration in parameter space
     * @param i The index of the iteration
     * @param length The length of the spectra array to be output
     * @return Pointer to the first element of the spectra array
     */
    const double* getSpectraDerivatives(size_t i, size_t* length) const {
        auto& spectra_i = get<2>(spectra[i]);
        *length = spectra_i.size();
        return spectra_i.data(); 
    }

    /**
     * @brief Get the potential values (eV) of the barrier of the i-th iteration in parameter space
     * @param i The index of the iteration
     * @param length The length of the spectra array to be output
     * @return Pointer to the first element of the spectra array
     */
    vector<double> getBarrierValues(const vector<double>& x, size_t paramsIndex = 0);

    void getBarrierValues(const double* x, double* potential, size_t size, size_t paramsIndex = 0);

    const double* getBarrierValues(const double* x, size_t size, size_t paramsIndex = 0) {
        return getBarrierValues(vector<double>(x, x + size), paramsIndex).data();
    }

    pair<double, double> getBarrierIntegrationLimits(size_t paramIndex = 0);

private:

    Config config; ///< Configuration object for the calculation
    /**
     * @brief The parameters for the calculation used in a certain iteration
     */
    struct ParamsForIteration{
        double workFunction = 4.5; ///< Work function in eV.
        double kT = .025; ///< Temperature in eV.
        double effectiveMass = 1.; ///< Effective mass of the electron.
        double bandDepth = 10.; ///< The depth of the electronic band in eV.
        double field = 5.; ///< The electric field in V/nm.
        double radius = 1.e5; ///< The radius of the emitter in nm.
        double gamma = 10.; ///< The gamma parameter of the general barrier model.
    };

    vector<double> fieldsVector = {5.}; ///< The electric field in V/nm, multiple values to iterate over.
    vector<double> radiiVector = {1.e5}; ///< The radius of the emitter in nm, multiple values to iterate over.
    vector<double> gammasVector = {10.}; ///< Gamma parameter of the general barrier model, multiple values to iterate over.
    vector<double> kTVector = {.025}; ///< Temperature in eV, multiple values to iterate over.
    vector<double> workFunctionVector = {4.5}; //< Work function in eV, multiple values to iterate over.
    vector<double> bandDepthVector = {10.}; //< Depth of the electronic band in eV, multiple values to iterate over.
    vector<double> effectiveMassVector = {1.}; ///< Effective mass of the electron, multiple values to iterate over.

    // double currentDensity = 0.; //< The current density (output) in A/nm^2.
    vector<double> currentDensityVector; ///< The current density (output) in A/nm^2, multiple values to iterate over.

    // double nottinghamHeat = 0.; //< The Nottingham heat (output) in W/nm^2.
    vector<double> nottinghamHeatVector; ///< The Nottingham heat (output) in W/nm^2, multiple values to iterate over.

    vector<tuple<vector<double>, vector<double>, vector<double>>> spectra; /**< The spectra (output) in A/nm^2/eV, multiple values to iterate over. */

    tbb::enumerable_thread_specific<unique_ptr<ModifiedSNBarrier>> threadLocalBarrier; ///< Thread-local instances of ModifiedSNBarrier
    tbb::enumerable_thread_specific<TransmissionSolver> threadLocalSolver; ///< Thread-local instances of TransmissionSolver
    tbb::enumerable_thread_specific<BandEmitter> threadLocalEmitter; ///< Thread-local instances of BandEmitter
    tbb::enumerable_thread_specific<ParamsForIteration> threadLocalParams; ///< Thread-local instances of ParamsForIteration


    /**
     * @brief Get the maximum number of iterations
     * @return The maximum number of iterations
     */
    size_t getMaxIterations();

    /**
     * @brief Set the parameters for a specific iteration
     * @param i The index of the iteration
     */
    void setParamsForIteration(size_t i = 0);
};

}

#endif // GETELEC_H