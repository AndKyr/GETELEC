/*
 * Config.cpp
 *
 *  Created on: 11.11.2016
 *      Author: veske
 */

#include <fstream>
#include <algorithm>

#include "Config.h"
#include <cassert>
#include <iostream>
#include <sstream>


using namespace std;

// Remove the noise from the beginning of the string
void Config::trim(string& str) {
    str.erase(0, str.find_first_of(comment_symbols + data_symbols));
}

void Config::read_all(const string& fname) {
    if (fname == "") return;
    file_name = fname;

    // Store the commands and their arguments
    parse_file(fname);

    for (auto& [key, value] : transmissionSolverParams.keyMap) {
        read_command("transmissionSolver." + key, value);
    }

    // Modify the parameters that are specified in input script
    for (auto& [key, value] : bandEmitterParams.keyMap) {
        read_command("bandEmitter." + key, value);
    }

}

void Config::print_all_params(){
    for (auto [key, value] : transmissionSolverParams.keyMap) {
        std::visit([key](auto* arg) {
            cout << "transmissionSolver." <<  key << " = ";
            using T = remove_cvref_t<decltype(*arg)>; // Remove const, volatile, and reference qualifiers
            if constexpr (is_same_v<T, vector<double>>)
                for (const auto& val : *arg) cout << val << " ";
            else if constexpr (is_same_v<T, vector<string>>)
                for (const auto& val : *arg) cout << val << " ";
            cout << endl;
        }, value);    
    }

    for (auto [key, value] : bandEmitterParams.keyMap) {
        std::visit([key](auto* arg) {
            cout << "transmissionSolver." <<  key << " = ";
            using T = remove_cvref_t<decltype(*arg)>; // Remove const, volatile, and reference qualifiers
            if constexpr (is_same_v<T, vector<double>>)
                for (const auto& val : *arg) cout << val << " ";
            else if constexpr (is_same_v<T, vector<string>>)
                for (const auto& val : *arg) cout << val << " ";
            cout << endl;
        }, value);    
    }
}

void Config::parse_file(const string& file_name) {
    ifstream file(file_name);

    if(!file.is_open()){
        std::cerr << "Config file " << file_name << " not found." << endl;
        return;
    }    
    string line;
    data.clear();

    // loop through the lines in a file
    while (getline(file, line)) {
        line += " "; // needed to find the end of line

        bool line_started = true;
        // store the command and its parameters from non-empty and non-pure-comment lines
        while(line.size() > 0) {
            trim(line);
            int i = line.find_first_not_of(data_symbols);
            if (i <= 0) break;

            if (line_started && line.substr(0, i) == "femocs_end") return;
            if (line_started) data.push_back({});
            if (line_started) {
                // force all the characters in a command to lower case
                string command = line.substr(0, i);
                std::transform(command.begin(), command.end(), command.begin(), ::tolower);
                data.back().push_back(command);
            } else {
                data.back().push_back( line.substr(0, i) );
            }

            line = line.substr(i);
            line_started = false;
        }
    }
}

void Config::check_obsolete(const string& command) {
    for (const vector<string>& cmd : data)
        if (cmd[0] == command) {
            cout << "Command '"  << command  << "' is obsolete! You can safely remove it!" << endl;
            return;
        }
}

void Config::check_changed(const string& command, const string& substitute) {
    for (const vector<string>& cmd : data)
        if (cmd[0] == command) {
            cout << "Command '" << command << "' has changed!"
                    " It is similar yet different to the command '" << substitute << "'!" << endl;
            return;
        }
}

int Config::read_command(string param, string& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            arg = str[1]; return 0;
        }
    return 1;
}

int Config::read_command(string param, bool& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is1(str[1]);
            istringstream is2(str[1]);
            bool result;
            // try to parse the bool argument in text format
            if (is1 >> std::boolalpha >> result) { arg = result; return 0; }
            // try to parse the bool argument in numeric format
            else if (is2 >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, unsigned int& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); int result;
            if (is >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, int& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); int result;
            if (is >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, double& arg) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); double result;
            if (is >> result) { arg = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, double& arg, double& arg2) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param) {
            istringstream is(str[1]); double result;
            if (is >> result) { arg = result; arg2 = result; return 0; }
            return 1;
        }
    return 1;
}

int Config::read_command(string param, vector<string>& args) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    int n_read_args = 0;
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param)
            for (unsigned i = 0; i < args.size() && i < (str.size()-1); ++i) {
                args[i] = str[i+1];
		    n_read_args++;
            }
    return n_read_args;
}

int Config::read_command(string param, vector<double>& args) {
    // force the parameter to lower case
    std::transform(param.begin(), param.end(), param.begin(), ::tolower);
    // loop through all the commands that were found from input script
    int n_read_args = 0;
    for (const vector<string>& str : data)
        if (str.size() >= 2 && str[0] == param)
            for (unsigned i = 0; i < args.size() && i < (str.size()-1); ++i) {
                istringstream is(str[i+1]);
                double result;
                if (is >> result) { args[i] = result; n_read_args++; }
            }
    return n_read_args;
}

void Config::print_data() {
    const int cmd_len = 20;

    for (const vector<string>& line : data) {
        for (const string& ln : line) {
            int str_len = ln.length();
            int whitespace_len = max(1, cmd_len - str_len);
            cout << ln << string(whitespace_len, ' ');
        }

        cout << endl;
    }
}