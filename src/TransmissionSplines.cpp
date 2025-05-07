#include "TransmissionSplines.h"
#include <iomanip>

namespace getelec{

vector<double> TransmissionSpline::evaluateSolution(double energy) const{
    double minEnergy = getMinimumSampleEnergy();
    double maxEnergy = getMaximumSampleEnergy();
    if (energy < minEnergy){
        vector<double> solutionOut = evaluateMultiple(minEnergy);
        vector<double> solutionDerivatives = evaluateDerivativeMultiple(minEnergy);
        for (int i = 0; i < solutionOut.size(); i++)
            solutionOut[i] += solutionDerivatives[i] * (energy - minEnergy);         
        return solutionOut;          
    } else if (energy > maxEnergy){
        vector<double> solutionOut = evaluateMultiple(maxEnergy);
        vector<double> solutionDerivatives = evaluateDerivativeMultiple(maxEnergy);
        for (int i = 0; i < solutionOut.size(); i++)
            solutionOut[i] += solutionDerivatives[i] * (energy - maxEnergy);
        return solutionOut;
    } else
        return evaluateMultiple(energy); 
}

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

int TransmissionSpline::writeSplineSolution(string filename, int nPoints){
    ofstream file(filename);

    vector<double> energyPoints = Utilities::linspace(getMinimumSampleEnergy(), getMaximumSampleEnergy(), nPoints);
    
    const int columnWidth = 16;
    vector<string> header = {"#energy", "sol[0]", "sol[1]", "sol[2]", "sol[3]", "sol[4]", "sol[5]", "interp[0]", "interp[1]", "interp[2]", "D[calc]", "D[interp]", "est. NED[calc]", "est.NED[interp]"};
    for(auto& s : header) file << setw(columnWidth) << s;
    file << endl;

    for (size_t i = 0; i < energyPoints.size(); i++){
        solver.setEnergyAndInitialValues(energyPoints[i]);
        solver.solveNoSave();
        const auto solutionVector = solver.getSolution();
        auto interpolatedValues = evaluateMultiple(energyPoints[i]);

        file << setw(columnWidth) << energyPoints[i];
        for(const auto& s : solutionVector) file << setw(columnWidth) << s;
        for(const auto& s : interpolatedValues) file << setw(columnWidth) << s;

        file << setw(columnWidth) << solver.getTransmissionProbabilityforWaveVector(12.) << setw(columnWidth) << getTransmissionProbability(energyPoints[i], 12.);
        file << setw(columnWidth) << solver.getEmissionEstimate(12., workFunction, kT) << setw(columnWidth) << normalEnergyDistributionEstimate(energyPoints[i]) << endl;
    }
    file.close();
    return 1;
}


int TransmissionSpline::writeSplineNodes(string filename){
    ofstream file(filename);
    file << "# energy, solutions[0], solutions[1], solutions[2], derivatives[0], derivatives[1], derivatives[2]" << endl;
    for (size_t i = 0; i < sampleEnergies.size(); i++){
        file << sampleEnergies[i] << " ";
        for (int j = 0; j < 3; j++){
            file << solutionSamples[j][i] << " " << solutionDerivativeSamples[j][i] << " ";
        }
        file << normalEnergyDistributionEstimate(sampleEnergies[i]) << endl; 
    }
    file.close();
    return 1;
}

void TransmissionSpline::sampleEnergyPoint(double energy, size_t index, bool bisect){

    if (index < 1 || index > sampleEnergies.size()) index = sampleEnergies.size(); //if the index is out of bounds, set it to the end
    sampleEnergies.emplace(sampleEnergies.begin() + index, energy);
    bisectList.emplace(bisectList.begin() + index, bisect);

    for (int i = 0; i < 3; i++){
        solutionSamples[i].emplace(solutionSamples[i].begin() + index, 0.);
        solutionDerivativeSamples[i].emplace(solutionDerivativeSamples[i].begin() + index, 0.);
    }

    if (solver.ensureBarrierDeepEnough(energy)){ // make sure the barrier is deep enough for the energy
        for(int i = 0; i < sampleEnergies.size(); i++){ // if the barrier depth has been changed, recalculate all the points
            calculateAndSetSampleValues(i);
        }
    } else // if nothing has changed, just calculate the new point
        calculateAndSetSampleValues(index); // calculate the values for the new point

}

void TransmissionSpline::calculateAndSetSampleValues(size_t index){
    double energy = sampleEnergies[index];

    const vector<double>& solutionVector = solver.calculateSolution(energy);
    for (int i = 0; i < 3; i++){
        solutionSamples[i][index] = solutionVector[i];
        solutionDerivativeSamples[i][index] = solutionVector[i + 3] * CONSTANTS.kConstant;
    }

    //check if solution is continuous and write the solution if not. Call it with plotSolutionData=true if you want data written and warning given.
    assert((checkSolutionSanity(solutionVector, energy, 0.5, false), true)); 
}

void TransmissionSpline::smartInitialSampling(double bandDepth, double effectiveMass){

    reset();

    assert((solver.getEnergyDerivativeLevel() > 0) && "The solver should be set with energyDerivatigveLevel>0 to get Hermite spline interpolation");

    //choose initial sampling point To be the Fermi, except if m*<0, (valence band), choose top of band
    double initialSampleEnergy = effectiveMass > 0. ? -workFunction : -workFunction - bandDepth;

    // make sure that the solver has energy derivatives
    solver.setXlimits(-initialSampleEnergy + 5.); // set the limits for the transmission solver. 

    if (solver.getXInitial() > 10.){
        writeSplineSolution();
        writeSplineNodes();
    }

    solver.getBarrier()->setBarrierTopFinder(true); // set barrier top finder
    sampleEnergyPoint(initialSampleEnergy); //sample the fermi level
    double topOfBarrier = solver.getBarrier()->getBarrierTop(); //get the barrier top
    solver.getBarrier()->setBarrierTopFinder(false); // reset the barrier top finder to off. No need any more

    double dLogD_dE_atInitialSample = -2. * solutionDerivativeSamples[2].back(); //Estimation for the transmission coefficient decay rate down from initial soample
    double logD_atInitialSample =  log(solver.getTransmissionProbabilityforWaveVector(12.)) ; //Estimation of log(D) at the initial sample 
    double numberOfDecayLengthsForRtol = -log(relativeTolerance);

    assert(logD_atInitialSample < 0. && "Transmission probability > 1. Something is off with the solver");
    solver.writeBarrierPlottingData();

    bool noTunnelingRegime = false; // it is true if  D(E < barrier top) ~= 0, i.e. no tunneling

    if (logD_atInitialSample > log(absoluteTolerance)){ //if D large enough to be worth the extra sample
        double lowDecayLength = 1. / dLogD_dE_atInitialSample;
        double lengthToDecay = min(numberOfDecayLengthsForRtol * lowDecayLength,  (logD_atInitialSample - log(absoluteTolerance)) * lowDecayLength);
        sampleEnergyPoint( initialSampleEnergy - lengthToDecay);
    } else {
        noTunnelingRegime = true;
    }

    maximumCurrentPosition = initialSampleEnergy;


    if (effectiveMass > 0.){ //if effective mass < 0., there is no need to sample above the top of the band
        double highDecayRate = 1./ kT - dLogD_dE_atInitialSample; //the estimated upwards decay rate at initialSample (if negative it means it goes upwards)
        double fermiToBarrier = topOfBarrier + workFunction;
        
        //If the spectra don't decay fast enough to have decayed before the barrier top (or don't decay at all), or Fermi is above barrier, sample the top of barrier and a point beyond
        if (highDecayRate < numberOfDecayLengthsForRtol / fermiToBarrier || fermiToBarrier < 1.e-3 ){ 

            if (fermiToBarrier > 0.2) sampleEnergyPoint(topOfBarrier); //sample the barrier top

            double logFDatBarrierTop = Utilities::logFermiDiracFunction(topOfBarrier + workFunction, kT); //estimate the current density at the barrier top
            if (logFDatBarrierTop < absoluteTolerance || noTunnelingRegime){ //if it is not worth going far or nasty noTunnelingRegime, sample nearby
                sampleEnergyPoint(topOfBarrier + 0.2);
            }

            if (logFDatBarrierTop >= absoluteTolerance){ //if sampling far beyond top barrier is worthit
                double lengthToDecay = min(numberOfDecayLengthsForRtol * kT,  (logFDatBarrierTop - log(absoluteTolerance)) * kT);
                sampleEnergyPoint(max(topOfBarrier, -workFunction) + lengthToDecay);
            }
            maximumCurrentPosition = topOfBarrier;
        }
        //If the spectra decay slow enough to be worth it, sample a above Ef. Force sample it we have only one sample point
        else if (highDecayRate < 100. ||  sampleEnergies.size() < 2){
            double pointToSample = -workFunction + numberOfDecayLengthsForRtol / highDecayRate;
            sampleEnergyPoint(pointToSample);
        }
    }

    sortAndInitialize();
}

int TransmissionSpline::findMaximumCurrentEstimate(int maxIterations, double rTol){
    auto lambdaMinimizer = [](double energy, void* params) -> double {
        TransmissionSpline* thisObject = (TransmissionSpline*) params;
        return -thisObject->normalEnergyDistributionEstimate(energy);
    }; // a lambda that gives the - of the NED
    double (* functionPointer)(double, void*) = lambdaMinimizer;//convert the lambda into raw function pointer
    gsl_function F = {functionPointer, this}; //define gsl function to be minimized
    

    gsl_set_error_handler_off(); // turn off GSL error handler
    int status = gsl_min_fminimizer_set(minimizer, &F, maximumCurrentPosition, getMinimumSampleEnergy(), getMaximumSampleEnergy());

    if (status == GSL_EINVAL){
        if (checkAllBelowTolerance())
            return 0; // all values are below the tolerance. No need to minimize
        else{
            return -1;
        }
    }

    //do minimization iterations
    for (int i = 0; i < maxIterations; i++){
        int status = gsl_min_fminimizer_iterate(minimizer); //perform minimization iteration
        
        if (status == GSL_EBADFUNC || status == GSL_FAILURE){ // something is messed up. Throw exception
            assert(false && "Step for finding max estimated current energy failed to iterate");
            throw std::runtime_error("Step for finding max estimated current energy failed to iterate");
        }
        //check if we have converged
        double minValue = gsl_min_fminimizer_f_minimum(minimizer);
        double maxValueInInterval = max(gsl_min_fminimizer_f_lower(minimizer), gsl_min_fminimizer_f_upper(minimizer));
        if (abs(minValue / maxValueInInterval - 1.) < rTol){ //if we have reached tolerance
            mamximumCurrentEstimate = - minValue;
            maximumCurrentPosition = gsl_min_fminimizer_minimum(minimizer);
            return i + 1;
        }
    }

    assert((writeSplineSolution(), writeSplineNodes(), false));
    return maxIterations;

}

int TransmissionSpline::refineSampling(){
    double testWaveVector = sqrt(7.) * CONSTANTS.sqrt2mOverHbar;
    int noInsertedSamples = 0;
    for (int i = 0; i < sampleEnergies.size() - 1; i++){
        if (bisectList[i]){
            double newEnergy = .5*(sampleEnergies[i] + sampleEnergies[i+1]);
            noInsertedSamples++;
            sampleEnergyPoint(newEnergy, ++i, true);

            double calculatedCurrentEstimate = solver.getEmissionEstimate(testWaveVector, workFunction, kT);
            double interpolatedCurrentEstimate = normalEnergyDistributionEstimate(newEnergy, testWaveVector);
            double error = abs(calculatedCurrentEstimate - interpolatedCurrentEstimate);
            double tolerance = absoluteTolerance + relativeTolerance * mamximumCurrentEstimate;
            if (error < tolerance){
                bisectList[i] = false;
                bisectList[i-1] = false;
            }
        }
    }
    return noInsertedSamples;
}

bool TransmissionSpline::checkSolutionSanity(const vector<double> &solutionVector, double energy, double absoluteTolerance, bool writeSolverPlottingData){
    // Check if the solution vector is finite
    if (!all_of(solutionVector.begin(), solutionVector.end(), [](double x) { return isfinite(x); })) {
        return false;
    }

    if (!isInitialized)
        return true; // if not initialized, we can't check the solution

    vector<double> expectedSolution = evaluateMultiple(energy);

    for (int i = 0; i < 3; i++){
        double diff = abs(solutionVector[i] - expectedSolution[i]);

        /** The following commented block was used for debugging the sampling. The bug is fixed, but it might be useful at some point again, so I leave it. */
        if (diff > absoluteTolerance && writeSolverPlottingData){
            cout << "GETELEC WARNING: The solution vector is not continuous. The difference is " << diff << " at energy " << energy << endl;
            cout << "Writing the problematic solution to file 'problematicWaveFunctionSolution.dat' and barrier to 'problematicBarrier.dat' " << endl;
            solver.solve(true);
            solver.writeSolution("problematicWaveFunctionSolution.dat");
            solver.writeBarrierPlottingData("problematicBarrier.dat", 0);

            solver.setEnergyAndInitialValues(energy + 0.1);
            solver.solve(true);
            solver.writeSolution("problematicWaveFunctionSolution_neighborEnergy.dat");

            writeSplineSolution("problematicSplineSolution.dat", 64);
            writeSplineNodes("problematicSplineNodes.dat");
            return false;
        }
    }

    return true;
}

void TransmissionSpline::refineSamplingToTolerance(int maxRefineSteps){
    int maxFoundStatus = findMaximumCurrentEstimate();

    if (maxFoundStatus == 0) 
        return; // all values are below the tolerance. No need to refine
        
    // unsigned attemptCounter = 0;
    // while(maxFoundStatus < 0){
    //     sampleEnergyPoint(getMaximumSampleEnergy() + kT);
    //     sortAndInitialize();
    //     maxFoundStatus = findMaximumCurrentEstimate();
    //     if (attemptCounter++ > 5){
    //         writeSplineSolution("problematicSplineSolution.dat", 64);
    //         writeSplineNodes("problematicSplineNodes.dat");
    //         assert(false && "Current density estimation sequence failed to locate maximum even after adding 5 points");
    //         throw std::runtime_error("Current density estimation sequence failed to locate maximum even after adding 5 points");
    //     }
    // }

    if(maxFoundStatus < 0){
        sampleEnergyPoint(getMaximumSampleEnergy() + kT);
        sortAndInitialize();
        assert((writeSplineSolution(), writeSplineNodes(), solver.writeBarrierPlottingData(), false) && "Failed to locate estimated current density maximum"); //if in debug mode, then write the splines to check them.
    }

    for (int i = 0; i < maxRefineSteps; i++){
        int insertedSamples = refineSampling();
        if (!insertedSamples)
            return;
        else
            sortAndInitialize();
    }
    cout << "GETELEC WARNING: Max allowed refine steps " << maxRefineSteps << " reached in spline interpolation without satisfying tolerance rtol: " 
        << relativeTolerance << " atol: " << absoluteTolerance << endl;
    
    assert((writeSplineSolution(), writeSplineNodes(), true)); //if in debug mode, then write the splines to check them.
}

void TransmissionSpline::sortAndInitialize(){

    // Lambda to compare the indices based on the sample energies
    // This lambda captures 'this' to access the member variables of the class
    auto comparisonLambda = [this](size_t i1, size_t i2) { return this->sampleEnergies[i1] < this->sampleEnergies[i2]; };

    vector<size_t> indices(sampleEnergies.size());
    iota(indices.begin(), indices.end(), 0); // Fill indices with 0, 1, ..., n-1
    // Sort the indices based on the sample energies
    sort(indices.begin(), indices.end(), comparisonLambda);

    bool isSorted = true;
    for (size_t i = 0; i < indices.size(); i++){
        if (indices[i] != i){
            isSorted = false;
            break;
        }
    }

    if (!isSorted){
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
    }

    initializeMultiple(sampleEnergies, solutionSamples, solutionDerivativeSamples);
    isInitialized = true;
}


}