#include "TransmissionSplines.h"
#include <gsl/gsl_min.h>

namespace getelec{

void TransmissionSpline::sampleUniform(double Emin, double Emax, int nPoints){
    if (nPoints < 2)
        throw std::invalid_argument("The number of points for the spline must be at least 2.");

    
    sampleEnergies = Utilities::linspace(Emin, Emax, nPoints);
    solutionSamples = vector<vector<double>>(3, vector<double>(nPoints, 0.));
    solutionDerivativeSamples = solutionSamples;


    for (size_t i = 0; i < sampleEnergies.size(); i++){
        solver.setEnergyAndInitialValues(sampleEnergies[i]);
        solver.solveNoSave();
        const auto& solutionVector = solver.getSolution();
        for (int j = 0; j < 3; j++){
            solutionSamples[j][i] = solutionVector[j];
            solutionDerivativeSamples[j][i] = solutionVector[j + 3] * CONSTANTS.kConstant;
        }
    }

    initializeMultiple(sampleEnergies, solutionSamples, solutionDerivativeSamples);    
}

void TransmissionSpline::writeSplineSolution(string filename, int nPoints){
    ofstream file(filename);

    vector<double> energyPoints = Utilities::linspace(getMinimumSampleEnergy(), getMaximumSampleEnergy(), nPoints);
    
    for (size_t i = 0; i < energyPoints.size(); i++){
        solver.setEnergyAndInitialValues(energyPoints[i]);
        solver.solveNoSave();
        const auto& solutionVector = solver.getSolution();
        auto interpolatedValues = evaluateMultiple(energyPoints[i]);
        file << energyPoints[i] << " " << solutionVector[0] << " " << solutionVector[1] << " " << solutionVector[2] << " ";
        file << interpolatedValues[0] << " " << interpolatedValues[1] << " " << interpolatedValues[2] << endl;
    }
    file.close();
}


void TransmissionSpline::writeSplineNodes(string filename){
    ofstream file(filename);
    for (size_t i = 0; i < sampleEnergies.size(); i++){
        file << sampleEnergies[i] << " ";
        for (int j = 0; j < 3; j++){
            file << solutionSamples[j][i] << " " << solutionDerivativeSamples[j][i] << " ";
        }
        file << endl; 
    }
    file.close();
}

void TransmissionSpline::sampleEnergyPoint(double energy, size_t index, bool bisect){

    // ensure that the energy is within the limits of the solver
    if (energy < minimumValidEnergyForSolver){
        double newBarrierDepth = -energy + 5.;
        while(energy < minimumValidEnergyForSolver){
            solver.setXlimits(newBarrierDepth); // set the limits for the transmission solver.
            newBarrierDepth += 5.;
            minimumValidEnergyForSolver = solver.minimumValidEnergy();
            if (newBarrierDepth > 100.){
                throw std::runtime_error("The solver's barrier depth was reduced below 100 eV without satisfying WKB. Something is wrong with the barrier shape");
            }
        }
    }
   
    solver.setEnergyAndInitialValues(energy);
    solver.solveNoSave();
    auto solutionVector = solver.getSolution();

    if (index < 1 || index > sampleEnergies.size()) index = sampleEnergies.size();

    sampleEnergies.emplace(sampleEnergies.begin() + index, energy);
    bisectList.emplace(bisectList.begin() + index, bisect);
    for (int i = 0; i < 3; i++){
        solutionSamples[i].emplace(solutionSamples[i].begin() + index, solutionVector[i]);
        solutionDerivativeSamples[i].emplace(solutionDerivativeSamples[i].begin() + index, solutionVector[i + 3] * CONSTANTS.kConstant);
    }
}



void TransmissionSpline::smartSampling(){

    assert((solver.getEnergyDerivativeLevel() > 0) && "The solver should be set with energyDerivatigveLevel>0 to get Hermite spline interpolation");

    // make sure that the solver has energy derivatives
    solver.setXlimits(workFunction + 5.); // set the limits for the transmission solver. 
    minimumValidEnergyForSolver = solver.minimumValidEnergy(); // get the minimum valid energy for the solver


    solver.getBarrier()->setBarrierTopFinder(true); // set barrier top finder
    sampleEnergyPoint(-workFunction); //sample the fermi level
    double topOfBarrier = solver.getBarrier()->getBarrierTop(); //get the barrier top
    solver.getBarrier()->setBarrierTopFinder(false); // reset the barrier top finder to off. No need any more


    //the extent by which the transmission coefficient decays below Fermi level
    double dLogD_dE_atFermi = -2. * solutionDerivativeSamples[2].back();
    double logD_atFermi = -solutionSamples[2].back();
    double lowDecayLength = 1. / dLogD_dE_atFermi;
    
    //TODO: take into account the absolute tolerance in determining the numberOfDecayLengths
    double numberOfDecayLengths = -log(relativeTolerance);



    if (abs(lowDecayLength) > 1.e-2){ //it decays slow enough to be worth it
        double pointToSample = -workFunction - min(numberOfDecayLengths * lowDecayLength, -(log(absoluteTolerance) - logD_atFermi)/dLogD_dE_atFermi);
        sampleEnergyPoint(pointToSample);
    }

    double highDecayRate = 1./ kT - dLogD_dE_atFermi;
    double fermiToBarrier = topOfBarrier + workFunction;

    
    //If the spectra don't decay fast enough to have decayed before the barrier top (or don't decay at all), sample the top of barrier and a point beyond
    if (highDecayRate < numberOfDecayLengths / fermiToBarrier){ 
        sampleEnergyPoint(topOfBarrier);
        sampleEnergyPoint(topOfBarrier + numberOfDecayLengths* kT);
    }
    //If the spectra decay slow enough to be worth it, sample a above Ef. Force sample it we have only one sample point
    else if (highDecayRate < 100. ||  sampleEnergies.size() < 2){
        double pointToSample = -workFunction + numberOfDecayLengths / highDecayRate;
        sampleEnergyPoint(pointToSample);
    }

    sortAndInitialize();
}

int TransmissionSpline::findMaximumCurrentEstimate(int maxIterations, double rTol){
    auto lambdaMinimizer = [](double energy, void* params) -> double {
        TransmissionSpline* thisObject = (TransmissionSpline*) params;
        return - thisObject->emissionCurrentEstimate(energy);
    }; // a lambda that gives the - of the NED
    double (* functionPointer)(double, void*) = lambdaMinimizer;//convert the lambda into raw function pointer
    gsl_function F = {functionPointer, this}; //define gsl function to be minimized
    
    // define and set the minimizer objects (GSL)
    const gsl_min_fminimizer_type* type = gsl_min_fminimizer_brent;
    gsl_min_fminimizer* minimizer = gsl_min_fminimizer_alloc(type);
    gsl_min_fminimizer_set(minimizer, &F, 0.5 * (getMinimumSampleEnergy() + getMaximumSampleEnergy()), getMinimumSampleEnergy(), getMaximumSampleEnergy());

    //do minimization iterations
    for (int i = 0; i < maxIterations; i++){
        int status = gsl_min_fminimizer_iterate(minimizer); //perform minimization iteration
        
        if (status == GSL_EBADFUNC || status == GSL_FAILURE) // something is messed up. Throw exception
            throw std::runtime_error("GSL minimizer failed to iterate. Something is wrong with the ");

        //check if we have converged
        double minValue = gsl_min_fminimizer_f_minimum(minimizer);
        double maxValueInInterval = max(gsl_min_fminimizer_f_lower(minimizer), gsl_min_fminimizer_f_upper(minimizer));
        if (abs(minValue / maxValueInInterval - 1.) < rTol){ //if we have reached tolerance
            mamximumCurrentEstimate = - minValue;
            maximumCurrentPosition = gsl_min_fminimizer_minimum(minimizer);
            gsl_min_fminimizer_free(minimizer);
            return i;
        }
    }

    //case when tolerance never met within the maxIterations
    gsl_min_fminimizer_free(minimizer);
    throw std::runtime_error("GSL minimizer failed to converge with relative tolerance " + std::to_string(rTol) + " and maxIterations " + std::to_string(maxIterations) + ".");
    return maxIterations;

}

int TransmissionSpline::refineSampling(){
    double testWaveVector = sqrt(7.) * CONSTANTS.sqrt2mOverHbar;
    int noInsertedSamples = 0;
    for (int i = 0; i < bisectList.size(); i++){
        if (bisectList[i]){
            double newEnergy = .5*(sampleEnergies[i] + sampleEnergies[i+1]);
            noInsertedSamples++;

            double calculatedCurrentEstimate = Utilities::logFermiDiracFunction(newEnergy + workFunction, kT) * solver.getTransmissionProbabilityforWaveVector(testWaveVector);
            double interpolatedCurrentEstimate = emissionCurrentEstimate(newEnergy, testWaveVector);
            double error = abs(calculatedCurrentEstimate - interpolatedCurrentEstimate);
            double tolerance = absoluteTolerance + relativeTolerance * mamximumCurrentEstimate;
            if (error > tolerance){
                sampleEnergyPoint(newEnergy, i+1, true);
            } else{
                sampleEnergyPoint(newEnergy, i+1, false);
                bisectList[i] = false;
            }
        }
    }
    return noInsertedSamples;
}

void TransmissionSpline::sortAndInitialize(){

    // Lambda to compare the indices based on the sample energies
    // This lambda captures 'this' to access the member variables of the class
    auto comparisonLambda = [this](size_t i1, size_t i2) { return this->sampleEnergies[i1] < this->sampleEnergies[i2]; };

    vector<size_t> indices(sampleEnergies.size());
    iota(indices.begin(), indices.end(), 0); // Fill indices with 0, 1, ..., n-1
    // Sort the indices based on the sample energies
    sort(indices.begin(), indices.end(), comparisonLambda);
    
    auto savedSolutions = solutionSamples;
    auto savedDerivatives = solutionDerivativeSamples;
    auto savedEnergies = sampleEnergies;
    // Reorder the solution and derivative samples based on the sorted indices

    for (size_t j = 0; j < sampleEnergies.size(); j++){
        sampleEnergies[j] = savedEnergies[indices[j]];
        for (int i = 0; i < 3; i++){
            solutionSamples[i][j] = savedSolutions[i][indices[j]];
            solutionDerivativeSamples[i][j] = savedDerivatives[i][indices[j]];
        }
    }

    initializeMultiple(sampleEnergies, solutionSamples, solutionDerivativeSamples);
}


}