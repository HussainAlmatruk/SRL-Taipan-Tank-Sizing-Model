function P = getInputs(vehicle_name, C)
%{
---------------------------------------------------------------------------
getInputs
  Router function that evaluates the vehicle-specific input file.
  The vehicle name is passed in from the main script.

Input:
  vehicle_name - [string] Name of the vehicle configuration to use (e.g., 'Taipan')
  C            - [struct] Constants struct from getConstants()

Output:
  P - [struct] Struct containing all vehicle and design input parameters
---------------------------------------------------------------------------
%}

    % Evaluate the selected vehicle function by its string name
    func = str2func(vehicle_name);
    P = func(C);
    
    % Store the vehicle name in the struct so output display plots and
    % logs can dynamically format the vehicle name.
    P.vehicle_name = vehicle_name;

end
