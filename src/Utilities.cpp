#include "Utilities.h"

vector<double> Utilities::linspace(double start, double end, int n) {
    vector<double> result(n);
    if (n == 1) {
        result[0] = start;
        return result;
    }

    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

double Utilities::fermiDiracFunction(double energy, double kT){
    if (energy > CONSTANTS.exponentLimit * kT)
        return exp(-energy / kT);
    else if(energy < - CONSTANTS.exponentLimit * kT)
        return 1. - exp(energy / kT);
    else
        return 1. / (1. + exp(energy / kT));
}

double Utilities::logFermiDiracFunction(double energy, double kT){
    if (energy > CONSTANTS.exponentLimit * kT)
        return exp(-energy / kT);
    else if(energy < - CONSTANTS.exponentLimit * kT)
        return -energy / kT + exp(energy / kT);
    else
        return log(1. + exp(-energy / kT));
}

int FunctionInterpolator::updateSpline(){
    gsl_spline_free(spline);
    gsl_interp_accel_free(accelerator);

    accelerator = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, samplingList.size());

    vector<double> x(samplingList.size());
    vector<double> y(samplingList.size());

    int i = 0;
    for (auto& val : samplingList){
        x[i] = val.x;
        y[i++] = val.y;
    }

    return gsl_spline_init(spline, x.data(), y.data(), x.size());
}

void FunctionInterpolator::initialize(double xInit, double xFinal, int numberOfInitialElements){
    samplingList.clear();
    if (spline) gsl_spline_free(spline);
    if (accelerator) gsl_interp_accel_free(accelerator);

    vector<double> xPoints = Utilities::linspace(xInit, xFinal, numberOfInitialElements);
    vector<double> yPoints(xPoints);
    for (int i = 0; i < xPoints.size(); i++){
        yPoints[i] = this->calculateYforX(xPoints[i]);
        samplingList.push_back(SplineElement(xPoints[i], yPoints[i], true));
    }
    samplingList.front().bisect = false;
    spline = gsl_spline_alloc(gsl_interp_cspline, xPoints.size());
    accelerator = gsl_interp_accel_alloc();
    gsl_spline_init(spline, xPoints.data(), yPoints.data(), xPoints.size());
}   

int FunctionInterpolator::refineSampling(){
    if (!spline) throw std::runtime_error("the class FunctionInterpolator is not properly initialized. Make sure you call initialize()");
    //TODO: check if initialized
    int numberOfAddedNodes = 0;
    for(auto it = samplingList.begin(); it != samplingList.end(); it++){
        if (it->bisect){
            numberOfAddedNodes++;
            double xNew = .5*(it->x + prev(it)->x);
            double yNew = this->calculateYforX(xNew);
            double error = this->calculateError(xNew, yNew);
            double tolerance = this->calculateTolerance(xNew, yNew);
            if (error < tolerance){
                samplingList.emplace(it, SplineElement(xNew, yNew, false));
                it->bisect = false;
            } else
                samplingList.emplace(it, SplineElement(xNew, yNew, true));
        }
    }
    return numberOfAddedNodes;
}

int FunctionInterpolator::refineToTolerance(int maxRefiningSteps){
    int numberOfSteps;
    for (numberOfSteps = 0; numberOfSteps < maxRefiningSteps; numberOfSteps++)
        if(refineSampling())
            updateSpline();
        else
            break;
    return numberOfSteps;
}

void FunctionInterpolator::writeSplineNodes(string filename){
    ofstream outFile(filename, ios::out);        
    for (SplineElement& element : samplingList)
        outFile << element.x << " " << element.y << endl;
}