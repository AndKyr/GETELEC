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

        xcFunctionParams = XCFunctionParams();
        xcFunctionParams.initializeKeyMap();

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

        /** @brief Parameters releated to the BandEmitter Class */
    struct XCFunctionParams : public ParamGroup {
        vector<double> dftXcPolynomial = {5.090468630345640e+00, -5.444183453319768e+01, 3.094936706449328e+02, -2.269318513321200e+02, 3.278425496065674e+03, -1.847834599204311e+05, 1.470254693085155e+06,  1.851826095135927e+06, -8.716343484652430e+07,  4.478932537960991e+08, 3.808767564452839e+07, -9.668478210868515e+09, 4.575620066264593e+10, -1.034493233247187e+11, 1.200702820177568e+11, -5.751842394314753e+10}; /**< Coefficients of the polynomial describing the DFT XC function  */
        vector<double> polynomialRange = {-0.136561638911387,  0.39}; /**< Range of validity of the polynomial  */
        double extensionPrefactor = 2.73378181454067;
        double extensionStartPoint = -0.24846083021481852;
        double transitionPoint = 0.19;
        double transisionWidth = 0.1;
        
        void initializeKeyMap(){
            name = "XCFunction";
            if (keyMap.size() == 0)
                keyMap = {
                    {"dftXcPolynomial", &dftXcPolynomial},
                    {"polynomialRange", &polynomialRange},
                    {"extensionPrefactor", &extensionPrefactor},
                    {"extensionStartPoint", &extensionStartPoint},
                    {"transitionPoint", &transitionPoint},
                    {"transisionWidth", &transisionWidth}
                };
        }
    } xcFunctionParams;

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
