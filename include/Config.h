#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>
#include <vector>
#include <map>
#include <variant>
#include <sstream>


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
        transmissionSolverParams.initializeKeyMap();

        bandEmitterParams = BandEmitterParams();
        bandEmitterParams.initializeKeyMap();
        read_all(fileName);
    }

    /**
     * @brief Read the configuration parameters from input script 
     * @param file_name path to the configuration file
    */
    void read_all(const string& file_name);

    /** @brief Print all the configuration parameters on a file with name
     * @param file_name path to the file where the configuration parameters are printed. If not provided, the parameters are printed on the console.
     */
    void print_all_params(const string& file_name = "printedGetelecParams.txt");

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

    // Function to handle a std::variant
    int read_command(string paramName, variant<string*, vector<string>*, bool*, int*, unsigned*, double*, vector<double>*>& var) {
        return std::visit([this, paramName](auto* arg)-> int {
            return read_command(paramName, *arg); // Calls the correct overload
        }, var);
    }

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

    struct ParamGroup {
        map<string, variant<string*, vector<string>*, bool*, int*, unsigned*, double*, vector<double>*>> keyMap;
        string name;

        string printParams(){
            ostringstream oss;
            for (auto [key, value] : keyMap) {
                visit([this, &key, &oss](auto* arg) {
                    oss << name << "." <<  key << " = ";
                    using T = remove_cvref_t<decltype(*arg)>; // Remove const, volatile, and reference qualifiers
                    if constexpr (is_same_v<T, vector<double>>)
                        for (const auto& val : *arg) oss << val << " ";
                    else if constexpr (is_same_v<T, vector<string>>)
                        for (const auto& val : *arg) oss << val << " ";
                    else
                        oss << *arg;
                    oss << endl;
                }, value);    
            }
            oss << endl;
            return oss.str();
        }
    };

    /** @brief Struct that contains all the configuration parameters related to the TransmissionSolver class */
    struct TransmissionSolverParams  : public ParamGroup {
        double relativeTolerance = 1.e-5;
        double absoluteTolerance = 1.e-5;
        string stepType = "rk8pd";
        int maxSteps = 4096;
        int minSteps = 64;
        int stepExpectedForInitialStep = 64;

        /** @brief Map of keywords to param references releated to the TransmissionSolver Class and their keywords. */
        void initializeKeyMap(){ 
            name = "transmissionSolver";
            if (keyMap.size() == 0)
                keyMap = {
                    {"relativeTolerance", &relativeTolerance},
                    {"absoluteTolerance", &absoluteTolerance},
                    {"stepType", &stepType},
                    {"maxSteps", &maxSteps},
                    {"minSteps", &minSteps},
                    {"stepExpectedForInitialStep", &stepExpectedForInitialStep}
                };
        }
    } transmissionSolverParams;

    /** @brief Parameters releated to the BandEmitter Class */
    struct BandEmitterParams : public ParamGroup {
        double relativeTolerance = 1.e-4;
        double absoluteTolerance = 1.e-12;
        int maxSteps = 4096;
        int minSteps = 16;
        int stepExpectedForInitialStep = 256;
        int maxAllowedRefiningSteps = 10;
        
        void initializeKeyMap(){
            name = "bandEmitter";
            if (keyMap.size() == 0)
                keyMap = {
                    {"relativeTolerance", &relativeTolerance},
                    {"absoluteTolerance", &absoluteTolerance},
                    {"maxSteps", &maxSteps},
                    {"minSteps", &minSteps},
                    {"stepExpectedForInitialStep", &stepExpectedForInitialStep},
                    {"maxAllowedRefiningSteps", &maxAllowedRefiningSteps}
                };
        }
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
