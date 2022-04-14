function [Current_Density,Nottingham_Heat] = Getelec_Matlab_Wrapper(material,Field, Radius, Gamma, Ec, Ef, Eg, T, me, mp)
    % Wrappers Geletec (in Python) with Matlab enviroment so it is callable
    % from COMSOL
    % material should be a number
    % input has to be an array
    % get rid of constant parameter
        %might be able to be e includ as text file 
        %but can be expensive. so inlcude constant parameters as constants
    %check for hdf5 format out
    % get rid material Eg gamma - make function just for semiconductor
    % 1) make the code first for metals (Field, Rad and T)
    % 2) make sure the whole thing runs on one machine (install Comsol)
    % 3) Tuesday 26th Veronika shows how is it done

    % self consistency should be done automatically, but we need to know how to do it

    %Checking and initialising Python environment
    pe = pyenv;
    if pe.Status == 'Loaded'
        % pass
    else 
        disp("Initialising Python enviroment");
        pyenv("ExecutionMode","OutOfProcess");
    end
    
    %Checking and initiliasing parallel processing
    if isempty(gcp('nocreate')) == 1
        disp('Initialising parpool');
        parpool
    else
        %pass
    end
    
    %Pointing MATLAB to the right files
    if count(py.sys.path, '/home/salva/getelec_remote/getelec_priv/python') == 0
        disp('Loading path')
        insert(py.sys.path, int64(0),'/home/salva/getelec_remote/getelec_priv/python');
    else
        %pass
    end
    
    %Importing GETELEC
    getelec = py.importlib.import_module('getelec_tabulator');
    
    %Reload environment if GETELEC modifed while pyenv loaded
    %py.importlib.reload(getelec)
    
    %Material selection
    if material == 'm'
        disp('Calculating J and Pn from metal emitter')
        getelec_data = cell(getelec.metal_emitter(Field, Radius, Gamma, Ef, T));
        Current_Density = getelec_data(1);
        Nottingham_Heat = getelec_data(2);
    elseif material == 's'
        disp('Calculating J and Pn from semiconductor emitter')
        getelec_data = cell(getelec.semiconductor_emitter(Field, Radius, Gamma, Ec, Ef, Eg, T, me, mp));
        Current_Density = getelec_data(1);
        Nottingham_Heat = getelec_data(2);
    else
        disp('Wrong definition. Insert "m" for metal emitters or "s" for semiconductor emitters');
    end
    
    %disp(Current_Density, Nottingham_Heat);

end