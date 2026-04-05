function writeResults( ...
    m_ox_kg, m_fuel_kg, m_pressurant_gas_kg,                                ...
    m_empty_ox_tank_kg, m_empty_fuel_tank_kg, m_empty_pressurant_tank_kg,   ...
    m_misc_and_plumbing_kg, m_final_kg, m_total_kg,                         ...
    d_ox_tank_m, t_cyl_ox_m, l_ox_tank_total_m,                            ...
    d_fuel_tank_m, t_cyl_fuel_m, l_fuel_tank_total_m,                      ...
    d_pressurant_tank_outer_m, r_pressurant_internal_m, r_pressurant_outer_m, ...
    t_pressurant_tank_m, l_pressurant_tank_total_m,                         ...
    l_total_vehicle_m, twr_ratio, delta_v_ms, f_drag_max_n, maxq_pa, h_max_m, ...
    acceleration_ms2, velocity_ms, altitude_m, t_s,                                          ...
    l_cyl_ox_m, l_cyl_fuel_m,                                                          ...
    v_ox_m3, v_fuel_m3, v_total_ox_tank_m3, v_total_fuel_tank_m3, v_pressurant_tank_m3, ...
    m_pressurant_gas_isothermal_kg, v_storage_isothermal_L,                             ...
    m_pressurant_gas_thermal_kg, v_storage_thermal_L,                                   ...
    m_pressurant_gas_v1_thermal_kg, v_storage_v1_thermal_L,                             ...
    T_ox_history_K, T_fuel_history_K, t_history_s,                                      ...
    T_pt_history_K, T_wall_history_K,                                                   ...
    p_pt_history_pa, m_pt_history_kg,                                                   ...
    m_ox_ullage_history_kg, m_fuel_ullage_history_kg,                                   ...
    P, C)

%{
---------------------------------------------------------------------------
writeResults
  Writes outputs from displayResults to an excel spreadsheet

Inputs: (all computed values from the calculation functions)
  See variable names - each maps directly to a computed result.
  P  - [struct] Design inputs from getInputs()   (uses constraint limits)
  C  - [struct] Physical constants from getConstants()  (uses unit conversions)

Outputs:
  None (writes to excel sheet)
---------------------------------------------------------------------------
%}

% --- 0.0 Setup Excel Template ---
if ~exist('Outputs', 'dir')
    mkdir('Outputs');
end
P.writefilename = fullfile('Outputs', sprintf('%s_Outputs.xlsx', P.vehicle_name));
copyfile('Outputs.xlsx', P.writefilename, 'f');

% --- Vehicle Configuration Title ---
writecell({sprintf('Vehicle Configuration: %s', P.vehicle_name)}, P.writefilename, 'Sheet', 1, 'Range', 'D2');
writecell({sprintf('Vehicle Configuration: %s', P.vehicle_name)}, P.writefilename, 'Sheet', 2, 'Range', 'C2');

% --- Mass Summary ---
writematrix([m_ox_kg;m_fuel_kg;m_pressurant_gas_kg;m_ox_kg + m_fuel_kg + m_pressurant_gas_kg;...
m_empty_ox_tank_kg;m_empty_fuel_tank_kg;m_empty_pressurant_tank_kg;m_misc_and_plumbing_kg;...
m_final_kg;m_ox_kg+m_empty_ox_tank_kg;m_fuel_kg+m_empty_fuel_tank_kg;...
m_empty_pressurant_tank_kg;m_total_kg],P.writefilename,'Sheet',1,'Range','D5:D17');

% --- Oxidiser Tank Summary ---
writematrix([d_ox_tank_m;d_ox_tank_m-2*t_cyl_ox_m;t_cyl_ox_m;l_ox_tank_total_m],P.writefilename,'Sheet',1,'Range','H5:H8');

% --- Fuel Tank Summary ---
writematrix([d_fuel_tank_m;d_fuel_tank_m-2*t_cyl_ox_m;t_cyl_fuel_m;l_fuel_tank_total_m],P.writefilename,'Sheet',1,'Range','H12:H15');

% --- Pressurant Tank Summary ---
writematrix([d_pressurant_tank_outer_m;2*r_pressurant_internal_m;t_pressurant_tank_m;l_pressurant_tank_total_m],P.writefilename,'Sheet',1,'Range','H19:H22');

% --- Flight Performance ---
writematrix([twr_ratio;delta_v_ms;f_drag_max_n;maxq_pa;h_max_m/.3048],P.writefilename,'Sheet',2,'Range','C5:C9');

runlength = find(acceleration_ms2(2:end)==0,1,'first')-1;
writematrix(runlength,P.writefilename,'Sheet',2,'Range','V3:V3');

runstring = strcat('W3:W',num2str(runlength+3));
writematrix(t_s',P.writefilename,'Sheet',2,'Range',runstring);
runstring = strcat('X3:X',num2str(runlength+3));
writematrix(acceleration_ms2',P.writefilename,'Sheet',2,'Range',runstring);
runstring = strcat('Z3:Z',num2str(runlength+3));
writematrix(velocity_ms',P.writefilename,'Sheet',2,'Range',runstring);
runstring = strcat('AA3:AA',num2str(runlength+3));
writematrix(altitude_m',P.writefilename,'Sheet',2,'Range',runstring);
% --- 3.0 Post-Processing (COM Server Image Injection & Sheet Organization) ---
try
    % Full paths required for COM server
    full_xlsx_path = fullfile(pwd, P.writefilename);
    plot_thermal   = fullfile(pwd, 'Outputs', sprintf('%s_Pressurant_Thermal_Analysis.png', P.vehicle_name));
    plot_stackup   = fullfile(pwd, 'Outputs', sprintf('%s_Vehicle_Stackup.png', P.vehicle_name));
    
    Excel = actxserver('Excel.Application');
    Excel.DisplayAlerts = false;
    Workbook = Excel.Workbooks.Open(full_xlsx_path);
    
    % Get existing sheets
    Sheet1 = Workbook.Sheets.Item(1);
    Sheet2 = Workbook.Sheets.Item(2);
    
    % Rename sheets clearly
    Sheet1.Name = 'Mass and Structure Results';
    Sheet2.Name = 'Flight Performance';
    
    % AutoFit text columns to prevent cutoff
    Sheet1.Range('D:D').EntireColumn.AutoFit();
    Sheet2.Range('C:C').EntireColumn.AutoFit();
    Sheet2.Range('B:B').EntireColumn.AutoFit();
    
    % Create a new dedicated sheet for Plots at the end of the workbook
    LastSheet = Workbook.Sheets.Item(Workbook.Sheets.Count);
    PlotSheet = Workbook.Sheets.Add([], LastSheet);
    PlotSheet.Name = 'Plots';
    
    % Insert shapes into the isolated Plot sheet (stacked vertically)
    PlotSheet.Shapes.AddPicture(plot_thermal, 0, 1, 10, 10, 650, 400);
    PlotSheet.Shapes.AddPicture(plot_stackup, 0, 1, 10, 430, 650, 450);
    
    % Reorder sheets to: Flight Performance, Mass and Structure, Plots
    Sheet2.Move(Sheet1);
    
    % Select Flight Performance so it opens first natively
    Sheet2.Activate();
    
    Workbook.Save();
    Workbook.Close();
    Excel.Quit();
    delete(Excel);
    disp('Successfully grouped and formatted the Excel Output.');
catch ME
    disp('Note: Could not finalize Excel via COM (Excel may not be installed or file is locked).');
    if exist('Excel', 'var') && isprop(Excel, 'Quit')
        Excel.Quit();
        delete(Excel);
    end
end

end