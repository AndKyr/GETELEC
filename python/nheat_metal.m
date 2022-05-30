function n_heat = nheat_metal(field)

    %Checking and initialising Python environment
    pe = pyenv;
    if pe.Status == 'Loaded'
        % pass
    else 
        disp("Initialising Python enviroment");
        pyenv("Version","/home/salva/Documents/getelec_priv/python/venv/bin/python3.8")
        pyenv("ExecutionMode","OutOfProcess");
    end
    
    %Pointing MATLAB to the right Python path
    if count(py.sys.path, '/home/salva/Documents/getelec_priv/python') == 0
        disp('Loading path')
        insert(py.sys.path, int64(0),'/home/salva/Documents/getelec_priv/python');
    else
        %pass
    end

    % WARNING: parfor returns a warming about missing libraries/functions +
    %and error from using dummy as indexing is not supported by the that
    %variable
    %Checking and initiliasing parallel processing
    %if isempty(gcp('nocreate')) == 1
    %    disp('Initialising parpool');
    %    parpool
    %else
        %pass
    %end
    
    %Importing GETELEC
    getelec = py.importlib.import_module('getelec_tabulator');

    %Initialising array to store our data
    n_heat = zeros(size(field));

    %Calculating current density from GETELEC
    disp('Calculating J metal emitter')
    for i = 1:length(field)
        n_heat(i) = getelec.heat_metal_emitter(field(i), 10, 10, 4.5, 300);
        %disp(size(Field(1,i)))
        %disp(current_density(i))
    end
    
end