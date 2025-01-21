#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>
#include <vector>


using namespace std;

/** @brief Class to initialize and read configuration parameters from configuration file */
class Config {
public:

    /** 
     * @brief initializes configuration parameters 
     * @param fileName path to the configuration file
    */
    Config(string fileName = "GetelecConfig.txt"){
        transmissionSolverParams = TransmissionSolverParams();
        bandEmitterParams = BandEmitterParams();
        read_all(fileName);
    }

    /**
     * @brief Read the configuration parameters from input script 
     * @param file_name path to the configuration file
    */
    void read_all(const string& file_name);

    /** 
     * @brief Read configuration parameter of type string
     * @param paramName the name of the parameter 
     * @param arg reference to the parameter to be read
    */
    int read_command(string paramName, string& arg);

    /** 
     * @brief Read configuration parameter of type boolean
     * @param paramName the name of the parameter 
     * @param arg reference to the parameter to be read  
    */
    int read_command(string paramName, bool& arg);

    /** 
     * @brief Read configuration parameter of type unsigned int
     * @param paramName the name of the parameter 
     * @param arg reference to the parameter to be read  
    */
    int read_command(string paramName, unsigned int& arg);

    /** 
     * @brief Read configuration parameter of type int
     * @param paramName the name of the parameter 
     * @param arg reference to the parameter to be read  
    */
    int read_command(string paramName, int& arg);

    /** 
     * @brief Read configuration parameter of type double
     * @param paramName the name of the parameter 
     * @param arg reference to the parameter to be read  
    */
    int read_command(string paramName, double& arg);


    /** 
     * @brief Read two configuration parameter of type double
     * @param paramName the name of the parameter 
     * @param arg reference to the first parameter to be read 
     * @param arg2 reference to the second parameter to be read 
    */
    int read_command(string paramName, double& arg, double& arg2);


    /** 
     * @brief Read a configuration parameter that is a series of strings
     * @param paramName the name of the parameter 
     * @param args reference to the vector of parameters to be read  
    */
    int read_command(string paramName, vector<string>& args);

    /** @brief Look up the configuration parameter with several double arguments */
    int read_command(string paramName, vector<double>& args);

    /** @brief Print the stored commands and parameters */
    void print_data();

    /** @brief Struct that contains all the configuration parameters related to the TransmissionSolver class */
    struct TransmissionSolverParams {
        double relativeTolerance = 1.e-5;
        double absoluteTolerance = 1.e-5;
        string stepType = "rk8pd";
        int maxSteps = 4096;
        int minSteps = 64;
        int stepExpectedForInitialStep = 64;
    } transmissionSolverParams;

    /** @brief Parameters releated to the BandEmitter Class */
    struct BandEmitterParams {
        double relativeTolerance = 1.e-4;
        double absoluteTolerance = 1.e-12;
        int maxSteps = 4096;
        int minSteps = 16;
        int stepExpectedForInitialStep = 256;
        int maxAllowedRefiningSteps = 10;
    } bandEmitterParams;

private:
    vector<vector<string>> data; ///< commands and their arguments found from the input script
    string file_name;            ///< path to configuration file

    const string comment_symbols = "!#%";
    const string data_symbols = "+-/*_.0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ()";

    /** @brief  Check for the obsolete commands from the buffered commands */
    void check_obsolete(const string& file_name);

    /** @brief Check for the obsolete commands that are similar to valid ones */
    void check_changed(const string& command, const string& substitute);

    /** @brief Read the commands and their arguments from the file and store them into the buffer */
    void parse_file(const string& file_name);

    /** @brief Remove the noise from the beginning of the string */
    void trim(string& str);
};

#endif /* CONFIG_H_ */
