function current_density = current_metal(field, radius, gamma, workFunction, temperature)
    
    %Checking and initialising Python environment
    pe = pyenv;
    if pe.Status ~= 'Loaded'
        disp("Initialising Python enviroment");
        pyenv("Version","/usr/bin/python3")
        pyenv("ExecutionMode","OutOfProcess");
    end
    

    %Importing GETELEC
    filePath = fileparts(mfilename('fullpath'));
    splitted = strsplit(filePath, filesep);
    getelecPath = append(join(splitted(1:end-2), filesep), filesep, 'src');
    getelecPath = getelecPath{1};
    pyrun(["import sys", sprintf("if '%s' not in sys.path: sys.path.append('%s')", getelecPath, getelecPath)]);
    pyrun(["matlabPaths = [x for x in sys.path if 'MATLAB' in x]",  "for x in matlabPaths: sys.path.remove(x)"]);
    
    getelecModule = py.importlib.import_module('getelecModel');

    gtm = getelecModule.GETELECModel();
    %py.importlib.reload(py.getelecModel);

    gtm.setParameters(field = py.numpy.array(field), radius = py.numpy.array(radius), ...
        gamma = py.numpy.array(gamma), workFunction = py.numpy.array(workFunction), ...
        temperature = py.numpy.array(temperature), emitterType = "metal");

    gtm.run(calculateCurrent = true, calculateNottinghamHeat = true, calculateSpectrum = false);

    
    current_density = gtm.getCurrentDensity();
end
