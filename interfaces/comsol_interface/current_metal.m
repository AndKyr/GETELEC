function current_density = current_metal(field,radius,gamma,workf,t)
    
    %Checking and initialising Python environment
    pe = pyenv;
    if pe.Status == 'Loaded'
        % pass
    else 
        disp("Initialising Python enviroment");
        %pyenv("Version","/home/salva/Documents/getelec_priv/python/venv/bin/python3.8")
        pyenv("Version","/usr/bin/python3.10")
        pyenv("ExecutionMode","OutOfProcess");
    end
    
    %Pointing MATLAB to the right Python path
    if count(py.sys.path, '/home/salva/Documents/getelec_priv/python') == 0
        disp('Loading path')
        insert(py.sys.path, int64(0),'/home/salva/Documents/getelec_priv/python');
    else
        %pass
    end

    %Importing GETELEC
    getelec = py.importlib.import_module('getelec_tabulator');

    %Calculating current density from GETELEC
    disp('Calculating J metal emitter')
    
    % field is transformed from a double to np, so GETELEC can handle it
    % GETELEC outcomes are transformed to double for COMSOL
    % in A/m^2
    current_density = double(getelec.current_metal_emitter(py.numpy.array(field),py.numpy.array(radius),py.numpy.array(gamma),py.numpy.array(workf), py.numpy.array(t))*1E18);

end
