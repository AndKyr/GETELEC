function current_density = current_metal(field)
    
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
    current_density = zeros(size(field));
    
    %Changing field magnitude (COMSOL - V/m) for GETELEC (V/A)
    %field = 2*field/1E8;

    %Calculating current density from GETELEC
    %disp('Calculating J metal emitter')
    for i = 1:length(field)
        current_density(1) = getelec.current_metal_emitter(field(i), 20, 10, 4.5, 300);
        %disp(size(Field(1,i)))
        %disp(current_density(i))
    end
    
    %Changing current density magnitude (GETELEC - A/nm^2) for COMSOL
    %(A/m^2)

    %current_density = current_density/1E18;
    %disp("TIME")
    %disp(toc)
end