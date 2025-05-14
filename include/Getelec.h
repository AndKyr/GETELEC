#ifndef GETELEC_H
#define GETELEC_H

#include <vector>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include "BandEmitter.h"
#include "TunnelingFunction.h"
#include "ConfigGetelec.h"
#include <numeric>


namespace getelec{

using namespace std;


enum class CalculationFlags : unsigned int {
    None            = 0,
    CurrentDensity  = 1 << 0, // 1
    NottinghamHeat  = 1 << 1, // 2
    TotalEnergyDistribution         = 1 << 2, // 4
    NormalEnergyDistribution        = 1 << 3, // 8
    ParallelEnergyDistribution      = 1 << 4, // 16
    TotalEnergyDistributionDerivatives = 1 << 5, //32
    All = CurrentDensity | NottinghamHeat | TotalEnergyDistribution | NormalEnergyDistribution | ParallelEnergyDistribution | TotalEnergyDistributionDerivatives
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
     * @param configFileName The name of the configuration input file
     */
    Getelec(string configFileName = "GetelecConfig.txt", string barrierType = "modifiedSN", mt19937* generator_ = NULL, int seed = -1) :
        config(Config(configFileName)),
        threadLocalBarrier([this, barrierType]() -> unique_ptr<ModifiedSNBarrier> {
            if (!barrierType.starts_with("dftXC")) {
                return make_unique<ModifiedSNBarrier>();
            } else {
                config.xcFunctionParams.name = barrierType;
                size_t nReadParams = config.readParamGroup(&config.xcFunctionParams);
                if (nReadParams < config.xcFunctionParams.keyMap.size()) {
                    cout << "WARNING: No parameters found for the DFT XC function named " << barrierType 
                            << " ; Using default values, which are for W(110)." << std::endl;
                    return make_unique<ModifiedSNBarrierWithDftXC>();
                } else {
                    return make_unique<ModifiedSNBarrierWithDftXC>(config.xcFunctionParams);
                }
            }
        }), 
        threadLocalSolver([this] { return TransmissionSolver(threadLocalBarrier.local().get(), config.transmissionSolverParams, 10., 1); }),
        threadLocalEmitter([this] { return BandEmitter(threadLocalSolver.local(), config.bandEmitterParams); })
    {
        if (generator_){
            generator = generator_;
            generatorIsInternallyConstructed = false;
        } else {
            if (seed < 0)
                generator = new mt19937(chrono::system_clock::now().time_since_epoch().count());
            else
                generator = new mt19937(seed);
            
            generatorIsInternallyConstructed = true;
        }
    }

    ~Getelec(){
        if (generator && generatorIsInternallyConstructed) {
            delete generator;
            generator = nullptr;
        }
    }

    void setGenerator(mt19937* generator_){
        if (generator && generatorIsInternallyConstructed) delete generator;
        generator = generator_;
        generatorIsInternallyConstructed = false;
    }

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
     * @brief Sets (within range) random parameters to all inputut variables. This is usable for testing purposes
     * @param numberOfParameters The number of different parameter sets to be set (length of vectors)
     */
    void setRandomParameters(unsigned numberOfParameters = 1);

    /**
     * @brief Calculate for the i-th element of the array of inputs
     * @param i The index of the element to calculate
     * @param flags Determines what quantities are to be calculated
     */
    void runIteration(size_t i = 0, CalculationFlags flags = CalculationFlags::CurrentDensity);

    /**
     * @brief Run the calculation
     * @param flags Determines what quantities are to be calculated
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
    // double calculateTransmissionProbability(double energy, double waveVector = -1., size_t paramsIndex = 0);

    /**
     * @brief Calculate the transmission probabilities for multiple energies and waveVectors
     * @param waveVectors The list of the waveVectors to calculate for (if empty, it defaults to 12. / nm)
     * @return The transmission probabilities
     * @note This method requires the calculateTransmissionForEnergies() to have been called in prior
     */
    vector<double> getTransmissionProbabilities(const vector<double>& waveVectors = {}) const;

    /**
     * @brief Calculate the transmission coefficients for multiple energies and waveVectors
     * @param waveVectors The list of the waveVectors to calculate for (if empty, it defaults to 12. / nm)
     * @return The transmission coefficients
     * @note This method requires the calculateTransmissionForEnergies() to have been called in prior
     */
    vector<gsl_complex> getTransmissionCoefficients(const vector<double>& waveVectors = {}) const;


    /**
     * @brief Calculate the transmission coefficient for multiple energies
     * @param energies The energy levels (eV)
     * @return The transmission coefficients
     * @note This method is relevant for multiple calculations of the transmission coefficient. It is faster than calculateTransmissionCoefficientForEnergy for multiple calculations on the same barrier. However, if you are iterating over many many energies, it might be better to use calculateTransmissionCoefficientForManyEnergies, which prepares the interpolator and then just interpolates.
     */
    vector<double> interpolateTransmissionProbabilities(const vector<double>& energies, const vector<double>& waveVectors = {}, size_t paramsIndex = 0);

    /**
     * Raw C-style getters for all input data arrays
     */
    const double* getFields(size_t* size) const { *size = fieldsVector.size(); return fieldsVector.data(); }
    const double* getRadii(size_t* size) const { *size = radiiVector.size(); return radiiVector.data(); }
    const double* getGammas(size_t* size) const { *size = gammasVector.size(); return gammasVector.data(); }
    const double* getkT(size_t* size) const { *size = kTVector.size(); return kTVector.data(); }
    const double* getWorkFunction(size_t* size) const { *size = workFunctionVector.size(); return workFunctionVector.data(); }
    const double* getBandDepth(size_t* size) const { *size = bandDepthVector.size(); return bandDepthVector.data(); }
    const double* getEffectiveMass(size_t* size) const { *size = effectiveMassVector.size(); return effectiveMassVector.data(); }

    /**
     * @brief Get the current density at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The current density in A/nm^2
     */
    double getCurrentDensity(unsigned i = 0) const { 
        assert(calculationStatusFlags & CalculationFlags::CurrentDensity && "Current density not calculated yet");
        assert(i < currentDensityVector.size() && "Index out of bounds");
        return currentDensityVector[i]; 
    }
    
    /**
     * @brief Get the Nottingham heat at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The Nottingham heat in W/nm^2
     */
    double getNottinghamHeat(unsigned i = 0) const { 
        assert(calculationStatusFlags & CalculationFlags::NottinghamHeat && "Nottingham heat not calculated yet");
        assert(i < nottinghamHeatVector.size() && "Index out of bounds");
        return nottinghamHeatVector[i]; 
    }
    
    /**
     * @brief Get the total energy distribution at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The spectra in A/nm^2/eV
     */
    const pair<vector<double>, vector<double>> getTotalEnergyDistribution(unsigned i = 0) const { 
        assert(calculationStatusFlags & CalculationFlags::TotalEnergyDistribution && "Total energy distribution not calculated yet");
        assert(i < totalEnergyDistributions.size() && "Index out of bounds");
        return totalEnergyDistributions[i]; 
    }

    /**
     * @brief Get the total energy distribution at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The spectra in A/nm^2/eV
     */
    const pair<vector<double>, vector<double>> getNormalEnergyDistribution(unsigned i = 0) const { 
        assert(calculationStatusFlags & CalculationFlags::NormalEnergyDistribution && "Normal energy distribution not calculated yet");
        assert(i < normalEnergyDistributions.size() && "Index out of bounds");
        return normalEnergyDistributions[i]; 
    }

    /**
     * @brief Get the total energy distribution at the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The spectra in A/nm^2/eV
     */
    const pair<vector<double>, vector<double>> getParallelEnergyDistribution(unsigned i = 0) const { 
        assert(calculationStatusFlags & CalculationFlags::ParallelEnergyDistribution && "Total energy distribution not calculated yet");
        assert(i < parallelEnergyDistributions.size() && "Index out of bounds");
        return parallelEnergyDistributions[i]; 
    }

    /**
     * @brief Get the total energy distribution derivatives of the i-th element of the array of inputs
     * @param i The index of the element to get
     * @return The spectra derivative in A/nm^2/eV/eV
     */
    const vector<double> getTotalEnergyDistributionDerivatives(unsigned i = 0) const { 
        assert(calculationStatusFlags & CalculationFlags::TotalEnergyDistributionDerivatives && "Total energy distribution derivatives not calculated yet");
        assert(i < totalEnergyDistributionsDerivatives.size() && "Index out of bounds");
        return totalEnergyDistributionsDerivatives[i]; 
    }
    
    
    /**
     * @brief Get all the resulting current densities as a vector
     * @return The vector
     * @note This is a copy of the internal vector, so it is safe to modify it.
     */
    vector<double> getCurrentDensities() const { 
        assert(calculationStatusFlags & CalculationFlags::CurrentDensity && "Current density not calculated yet");
        return currentDensityVector; 
    }
    
    /**
     * @brief Get all the resulting current densities as a vector (used for C and python interface)
     * @param size The size of the vector to be output
     * @return Pointer to the first element of the vector
     */
    const double* getCurrentDensities(size_t* size) const { 
        assert(calculationStatusFlags & CalculationFlags::CurrentDensity && "Current density not calculated yet");
        *size = currentDensityVector.size();
        return currentDensityVector.data(); 
    }

    /**
     * @brief Get all the resulting Nottingham heats as a vector
     * @return The vector
     * @note This is a copy of the internal vector, so it is safe to modify it.
     */
    vector<double> getNottinghamHeats() const { 
        assert(calculationStatusFlags & CalculationFlags::NottinghamHeat && "Nottingham heat not calculated yet");
        return nottinghamHeatVector;
    }

    /**
     * @brief Get all the resulting Nottingham heats as a vector (used for C and python interface)
     * @param size The size of the vector to be output
     * @return Pointer to the first element of the vector
     */
    const double* getNottinghamHeats(size_t* size) const { 
        assert(calculationStatusFlags & CalculationFlags::NottinghamHeat && "Nottingham heat not calculated yet");
        *size = nottinghamHeatVector.size();
        return nottinghamHeatVector.data(); 
    }

    /**
     * @brief Get the spectra abssicae (energies in eV) of the i-th iteration in parameter space
     * @param i The index of the iteration
     * @param length The length of the spectra array to be output
     * @param spectraType The type of spectra to be returned. 'T' for TED, 'N' for NED, 'P' for PED
     * @return Pointer to the first element of the spectra array
     */
    const double* getSpectraEnergies(size_t i, size_t* length, char spectraType = 'T') const; 

    /**
     * @brief Get the spectra values (A / nm^2 / eV) of the i-th iteration in parameter space
     * @param i The index of the iteration
     * @param length The length of the spectra array to be output
     * @param spectraType The type of spectra to be returned. 'T' for TED, 'N' for NED, 'P' for PED, 'D' for TED derivatives (if D in A/nm^2/eV/eV)
     * @return Pointer to the first element of the spectra array
     */
    const double* getSpectraValues(size_t i, size_t* length, char spectraType = 'T') const;

    /** @TODO: Consider exposing the barrier object rather than writing endless wrapper functions here */
    vector<double> getBarrierValues(const vector<double>& x, size_t paramsIndex = 0);

    void getBarrierValues(const double* x, double* potential, size_t size, size_t paramsIndex = 0);

    const double* getBarrierValues(const double* x, size_t size, size_t paramsIndex = 0) {
        return getBarrierValues(vector<double>(x, x + size), paramsIndex).data();
    }

    pair<double, double> getBarrierIntegrationLimits(size_t paramIndex = 0);

    const CalculationFlags& getCalculationStatusFlags() const { return calculationStatusFlags; }

    /**
     * @brief writes all the calculates spectra for the calculation paramIndex into corresponding files
     * @param paramIndex The index of the data to be written
     */
    void writeSpectraToFiles(size_t paramIndex) const;

    /**
     * @brief set the writing flag (if true it forces all plotting data to be written to files)
     * @param flag The flag to set
     * @note Set this flag to true for debugging only
     */
    void setFileWriteFlag(bool flag) { doWritePlotFiles = flag; }

    /**
     * @brief calculates the transmission solution for a list of energies
     * @param energies The list of (normal) energies to calculate transmission for (eV, counting from fermi level)
     * @param paramsIndex The index of the parameter list to be used
     * @param forceCalculate Flag to force full calculation for every energy rather than inteprolation (default is fause)
     * @note If the number of requested energie sis > 32 , an interpolator is prepared and interpolated
     */
    void calculateTransmissionForEnergies(const vector<double>& energies, size_t paramsIndex = 0, bool forceCalculate = false);

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

    vector<pair<vector<double>, vector<double>>> totalEnergyDistributions; ///< The total energy distributions (output) in A/nm^2/eV.
    vector<pair<vector<double>, vector<double>>> normalEnergyDistributions; ///< The normal energy distributions (output) in A/nm^2/eV.
    vector<pair<vector<double>, vector<double>>> parallelEnergyDistributions; ///< The parallel energy distributions (output) in A/nm^2/eV.
    vector<vector<double>> totalEnergyDistributionsDerivatives; ///< The total energy distributions derivatives (A / nm^2 / eV / eV)
    vector<vector<double>> transmissionSolutions; ///< Stores the solutions to be used for transmission calculations of a single paramsIndex

    tbb::enumerable_thread_specific<unique_ptr<ModifiedSNBarrier>> threadLocalBarrier; ///< Thread-local instances of ModifiedSNBarrier
    tbb::enumerable_thread_specific<TransmissionSolver> threadLocalSolver; ///< Thread-local instances of TransmissionSolver
    tbb::enumerable_thread_specific<BandEmitter> threadLocalEmitter; ///< Thread-local instances of BandEmitter
    tbb::enumerable_thread_specific<ParamsForIteration> threadLocalParams; ///< Thread-local instances of ParamsForIteration
    CalculationFlags calculationStatusFlags; ///< Flag that shows which output data is available
    mt19937* generator = NULL; ///< Random number generator for setting random parameters for testing.
    bool generatorIsInternallyConstructed; ///< Keeps track of whether getelec is responsible of freeing the RNG
    bool doWritePlotFiles = false; ///< Determines whether all output is written.



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


    /**
     * @brief Runs the calculation iteration in the case when effectiveMass ~= 1. (within a given tolerance determined in config)
     * @param i The iteration to run for
     * @param flags The flags determining which quantities are to be calculated
     */
    void runIterationEffectiveMassUnity(size_t i = 0, CalculationFlags flags = CalculationFlags::CurrentDensity);

    /**
     * @brief Runs the calculation iteration in the case when effectiveMass != 1. (out of a given tolerance determined in config)
     * @param i The iteration to run for
     * @param flags The flags determining which quantities are to be calculated
     */
    void runIterationEffectiveMassNonUnity(size_t i = 0, CalculationFlags flags = CalculationFlags::CurrentDensity);
};

}

#endif // GETELEC_H