
#include <fstream>
#include <chrono>
#include <string>
#include <sstream>
#include "Getelec.h"

using namespace std;

static string error;
static getelec::Getelec* globalGetelecObj = nullptr;
static ofstream logFile;
static bool verbose = false;

static char* allocatedErrorStrForOutput = nullptr;

void logTimeStamp(){
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    logFile << "Time = " << std::ctime(&now_time);
}

int terminateGetelec(){
    if (!logFile.is_open()) {
        error = "log file is not open";
        return 0;
    }

    logTimeStamp();
    logFile << "Terminating: closing log file and deleting getelec object and error string" << endl << endl;

    logFile.close();
    if (allocatedErrorStrForOutput){
        free(allocatedErrorStrForOutput);
        allocatedErrorStrForOutput = nullptr;
    }
    
    if (globalGetelecObj)
        delete globalGetelecObj;
    globalGetelecObj = nullptr;
    return 1;
}

extern "C" {


// Initialization function for COMSOL wrapper
int init(const char *str) {
    string inputString = str;

    string inputStringLower = inputString;
    transform(inputStringLower.begin(), inputStringLower.end(), inputStringLower.begin(), ::tolower);
    if (inputStringLower.find("verbose") != std::string::npos)
        verbose = true;
    else
        verbose = false;

    stringstream ss(inputString);

    string logFilePath = "ComsolExternalLibrary.log";
    string barrierType = "modifiedSN";

    
    getline(ss, logFilePath, ',');
    getline(ss, barrierType, ',');
    if (logFilePath.empty())
        logFilePath = "ComsolExternalLibrary.log";
    
    if(barrierType.empty())
        barrierType = "modifiedSN";

    logFile.open(logFilePath, ios::app);
    if (!logFile.is_open()) {
        error = "Error opening log file " + logFilePath + " . Opening default log file: ./ComsolExternalLibrary.log ";
        logFile.open("ComsolExternalLibrary.log", ios::app);
        if (!logFile.is_open()) {
            error += "Error opening standard log file ./ComsolExternalLibrary.log \n";
            return 0;
        }
        logFile << error << endl;
    }


    logTimeStamp();
    if (!globalGetelecObj){
        globalGetelecObj = new getelec::Getelec();
        logFile << "Initialized GETELEC object at" << globalGetelecObj << " with barrier type " << barrierType << endl;
    } else{
        logFile << "Getelec Object already initialized at RAM address" << globalGetelecObj << endl;
    }

    if (verbose)
        logFile << "Verbosity mode: Verbose" << endl;
    else
        logFile << "Verbosity mode: silent" << endl;

    error.clear();

    logFile << endl;
    return 1;
}


// Error message function for COMSOL wrapper
const char* getLastError() {

    if (logFile.is_open()) {
        logTimeStamp();
        logFile << "Called getLastError. ERROR: " << error << endl << endl;
    } else {
        error += " ERROR: log file is not open";
    }

    if (allocatedErrorStrForOutput)
        free(allocatedErrorStrForOutput);

    allocatedErrorStrForOutput = (char*) malloc(error.size() + 1);
    std::strcpy(allocatedErrorStrForOutput, error.c_str());

    return allocatedErrorStrForOutput;
}
    
int eval(const char *func,
                    int nArgs,
                    const double **inReal,
                    const double **inImag,
                    int blockSize,
                    double *outReal,
                    double *outImag) 
{
    if (!logFile.is_open()) {
        error = "error opening log file";
        return 0;
    }
    string functionStr = func;


    if (verbose){
        logTimeStamp();
        logFile << "Called eval with function string: " << functionStr << ", nArgs= " << nArgs << "and blockSize= " << blockSize << endl;
        logFile << "Input data pointers are: inReal: " << inReal << " inImag: " << inImag << " outReal: " << outReal << " outImag: " << outImag << endl;

    }

    if (functionStr == "terminate"){
        return terminateGetelec();
    }

    if (functionStr != "getelec") {
        error = "Unknown function name string: " + functionStr  + "\n";
        logFile << "error message = " << error << endl;
        return 0;
    }

    if (nArgs > 7 || nArgs < 0) {
        error = "Invalid number of arguments. Expected is between 0 and 7";
        logFile << "Error with message: " << error << endl;
        return 0;
    }

    if (verbose){
        logFile << "Printing input data..." << endl;
        for (int i = 0; i < nArgs; i++){
            logFile << "inReal[" << i << "] @ " << inReal[i];
            if (inImag) logFile << "\t inImag[" << i << "] @ " << inImag[i] <<" \t";
        }
        logFile << endl;
        for (int j = 0; j < blockSize; j++){
            for (int i = 0; i < nArgs; i++){
                logFile << inReal[i][j];
                if (inImag && inImag[i]) logFile << " \t " << inImag[i][j] << " \t";
            }
            logFile << endl;
        }
        logFile << endl;
    }

    if (nArgs >= 1){
        if (verbose) logFile << "setting field and other parameters with blockSize: " << blockSize << endl;
        globalGetelecObj->setField(inReal[0], blockSize);
    }
    if (nArgs >= 2)
        globalGetelecObj->setRadius(inReal[1], blockSize);
    if (nArgs >= 3)
        globalGetelecObj->setGamma(inReal[2], blockSize);
    if (nArgs >= 4)
        globalGetelecObj->setkT(inReal[3], blockSize);
    if (nArgs >= 5)
        globalGetelecObj->setWorkFunction(inReal[4], blockSize);
    if (nArgs >= 6)
        globalGetelecObj->setEffectiveMass(inReal[5], blockSize);
    if (nArgs >= 7)
        globalGetelecObj->setBandDepth(inReal[6], blockSize);
    
    if (verbose)
        logFile << "running GETELEC and extracting current density and Nottingham heat" << endl;

    globalGetelecObj->run(getelec::CalculationFlags::CurrentDensity | getelec::CalculationFlags::NottinghamHeat);
    
    size_t outSize;
    const double* currentDensity = globalGetelecObj->getCurrentDensities(&outSize);
    const double* nottinghamHeat = globalGetelecObj->getNottinghamHeats(&outSize);

    if (outSize != blockSize) {
        error = "Output size does not match input size";
        logFile << "error message = " << error << endl;
        return 0;
    }

    if (verbose) logFile << "copying current density data into real output" << endl;

    for (size_t i = 0; i < blockSize; i++)
        outReal[i] = currentDensity[i];

    if (verbose) logFile << "copying Nottingham heat data into imaginary output" << endl;
    for (size_t i = 0; i < blockSize; i++)
        outImag[i] = nottinghamHeat[i];

    if (verbose){
        logFile << "outReal " << "outImag" << endl;
        for (size_t i = 0; i < blockSize; i++)
            logFile << outReal[i] << " " << outImag[i] << endl;
    }
    
    logFile << endl;
    return 1;
}

} // extern "C"