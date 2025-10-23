%{
---------------------------------------------------------------------------
CU Sounding Rocket Lab - Liquid Propulsion Team
Tank Sizing and Vehicle Mass Model

Description:
This script performs a trade study for the design of the propellant and
pressurant tanks for the Taipan liquid rocket engine.

Authors: Hussain Almatruk, Jonathan Forte
Last Updated: 10/21/2025
---------------------------------------------------------------------------
%}

%% 0.0 - SETUP
% Clears workspace, command window, and closes all figures.
clear; clc; close all;

%% 0.1 - CONSTANTS
g_earth_ms2         = 9.80665;       % [m/s^2] Standard gravitational acceleration
r_universal_jmolk   = 8.3145;      % [J/(mol*K)] Universal gas constant

%% 1.0 - INPUTS

% --- 1.1 Engine & Performance ---
f_thrust_n      = 750*4.44822;  % [N] Engine thrust
i_sp_s          = 287.458;      % [s] Specific impulse
t_burn_s        = 10;           % [s] Total burn time
o_f_ratio       = 2.5;          % [unitless] Oxidizer-to-Fuel mass ratio

% --- 1.2 Propellants & Pressurant ---
ox_density_kgm3   = 1141;                   % [kg/m^3] Density of Liquid Oxygen
fuel_density_kgm3 = 820;                    % [kg/m^3] Density of Jet-A
pressurant_molar_mass_kgmol = 0.0280134;    % [kg/mol] Molar mass of Nitrogen (N2)
residual_fraction = 0.1;                    % %%% TEMPORARY VALUE %%% [~] Percentage of fuel left over in tank after burnout

p_op_ox_tank_pa = 850*6894.76;              % %%% TEMPORARY VALUE %%% [Pa] Operating pressure of LOX tank TODO: Define this value
p_op_fuel_tank_pa = 850*6894.76;            % %%% TEMPORARY VALUE %%% [Pa] Operating pressure of Fuel tank TODO: Define this value
p_storage_pressurant_pa = 3000*6894.76;     % %%% TEMPORARY VALUE %%% [Pa] The pressure the Nitrogen is stored in dedicated tank (MEOP). TODO: Define this value
ox_temp_k         = 90;                     % %%% TEMPORARY VALUE %%% [K] Temperature of LOX in tank TODO: Define this value
fuel_temp_k       = 294;                    % %%% TEMPORARY VALUE %%% [K] Temperature of Jet-A in tank (Ambient probably)
pressurant_temp_k = 294;                    % %%% TEMPORARY VALUE %%% [K] Temperature of pressurant gas in its tank

% --- 1.3 Vehicle Geometry & Materials ---
% Assumption: Both LOX and Fuel tanks are cylinders of the same diameter
D = 6*2.54/100;                             %  %%% TEMPORARY VALUE %%% [m] Outer diameter of vehicle
d_ox_tank_m     = 5*2.54/100;               %  %%% TEMPORARY VALUE %%% [m] Outer diameter of the LOX tank
d_fuel_tank_m   = 5*2.54/100;               %  %%% TEMPORARY VALUE %%% [m] Outer diameter of the Fuel tank
d_pressurant_tank_allowed_m = 5*2.54/100;   %  %%% TEMPORARY VALUE %%% [m] Maximum allowed outer diameter of the Pressurant tank 
C_D = .5;                                   %  %%% TEMPORARY VALUE %%% [~] Vehicle drag coefficient

% Material Properties for LOX Tank 
material_density_ox_kgm3         = 2840;      % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
material_allowable_stress_ox_pa  = 2.90*10^8; % %%% TEMPORARY VALUE %%%  [Pa] TODO: Define this value

% Material Properties for Fuel Tank 
material_density_fuel_kgm3       = 2840; % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
material_allowable_stress_fuel_pa= 2.90*10^8; % %%% TEMPORARY VALUE %%%  [Pa] TODO: Define this value

% Material Properties for Pressurant Tank 
material_density_liner_kgm3             = 2840;     % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
t_liner_m                                 = 0.003;    % %%% TEMPORARY VALUE %%%  [m] Thickness of COPV liner 
material_density_pressurant_kgm3        = 1800;     % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
material_allowable_stress_pressurant_pa = 3.5*10^9; % %%% TEMPORARY VALUE %%%  [Pa] TODO: Define this value

% --- 1.4 Design Margins & Factors ---
safety_factor   = 1.5;        % [unitless] Safety factor for pressure vessels
% Joint efficiency can differ based on tank material and welding process
joint_efficiency_ox_tank        = 0.8; % %%% TEMPORARY VALUE %%%  [unitless] TODO: Define this value 
joint_efficiency_fuel_tank      = 0.8; % %%% TEMPORARY VALUE %%%  [unitless] TODO: Define this value 
joint_efficiency_pressurant_tank= 0.8; % %%% TEMPORARY VALUE %%%  [unitless] TODO: Define this value 

corrosion_allowance_m = 0.001;       % %%% TEMPORARY VALUE %%%  [m] Extra thickness for material degradation
ullage_fraction_ox       = 0.1;      % %%% TEMPORARY VALUE %%%  [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)
ullage_fraction_fuel     = 0.1;      % %%% TEMPORARY VALUE %%%  [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)

% --- 1.5 Estimated Masses (Non-Calculated) ---
m_misc_kg     = 25; % [kg]  %%% TEMPORARY VALUE %%% TODO: Estimate mass of payload, structure, fins, avionics, recovery
m_plumbing_kg = 10; % [kg]  %%% TEMPORARY VALUE %%% TODO: Estimate mass of valves and plumbing

% --- 1.6 Design Constraints ---

l_airframe_max_m = 3.048; %[m] Maximum allowable vehicle length TODO: Define this value 

%% 2.0 - CALCULATIONS
% This section should not be modified unless equations are being updated or the initial code is done and is now being itirated on to optimize values

% --- 2.1 - Propellant & Flow Rate Analysis ---

% Calculate Prop Mass Flow Rates (Eq. 1,2,3)
m_dot_total_kgs = f_thrust_n/(i_sp_s*g_earth_ms2);          % [kg/s] Total propellant mass flow rate
m_dot_ox_kgs = m_dot_total_kgs*(o_f_ratio/(1+o_f_ratio));   % [kg/s] Oxidizer mass flow rate 
m_dot_fuel_kgs = m_dot_total_kgs/(1+o_f_ratio);             % [kg/s] Fuel mass flow rate 

% Calculate Prop Masses (Eq. 4)
m_ox_kg = m_dot_ox_kgs*t_burn_s/(1-residual_fraction);        % [kg] Total mass of oxidizer
m_fuel_kg = m_dot_fuel_kgs *t_burn_s/(1-residual_fraction);   % [kg] Total mass of fuel

% --- 2.2 - Tank Sizing & Mass ---

% Calculate Volume of Props (Eq. 5)
v_ox_m3 = m_ox_kg/ox_density_kgm3;          % [m^3] Volume of liquid oxidizer
v_fuel_m3 = m_fuel_kg/fuel_density_kgm3;    % [m^3] Volume of liquid fuel

% Calculate Internal Volume of Prop Tanks (Volume of prop in addition to ullage) (Eq. 6)
v_total_ox_tank_m3 = v_ox_m3/(1-ullage_fraction_ox);       % [m^3] Total internal volume of LOX tank
v_total_fuel_tank_m3 = v_fuel_m3/(1-ullage_fraction_fuel);   % [m^3] Total internal volume of Fuel tank

% Calculate tank outer radii (Basic Geometry)
r_ox_tank_m = d_ox_tank_m/2;     % [m] Outer radius of cylindrical ox tank
r_fuel_tank_m = d_fuel_tank_m/2; % [m] Outer radius of cylindrical fuel tank
r_pressurant_tank_allowed_m = d_pressurant_tank_allowed_m/2; % [m] Maximum allowed outer radius of pressurant tank

% Calculate tank design pressures (Eq. 10)
p_design_ox_pa = p_op_ox_tank_pa * safety_factor;     % [Pa] Design pressure for LOX tank
p_design_fuel_pa = p_op_fuel_tank_pa * safety_factor; % [Pa] Design pressure for Fuel tank

% Calculate end cap wall thicknesses (Eq. 12, Spherical Shell Stress)
t_caps_ox_m = (p_design_ox_pa * d_ox_tank_m) / (4 * material_allowable_stress_ox_pa * joint_efficiency_ox_tank) + corrosion_allowance_m; % [m] LOX tank end cap thickness
t_caps_fuel_m = (p_design_fuel_pa * d_fuel_tank_m) / (4 * material_allowable_stress_fuel_pa * joint_efficiency_fuel_tank) + corrosion_allowance_m; % [m] Fuel tank end cap thickness

% Calculate cylinder wall thicknesses (Eq. 11, ASME Hoop Stress)
t_cyl_ox_m = (p_design_ox_pa * d_ox_tank_m) / (2 * material_allowable_stress_ox_pa * joint_efficiency_ox_tank) + corrosion_allowance_m; % [m] LOX tank cylinder thickness
t_cyl_fuel_m = (p_design_fuel_pa * d_fuel_tank_m) / (2 * material_allowable_stress_fuel_pa * joint_efficiency_fuel_tank) + corrosion_allowance_m; % [m] Fuel tank cylinder thickness

% Calculate end cap volumes (Eq. 7)
v_caps_ox_m3 = (4/3)*pi*(r_ox_tank_m - t_caps_ox_m)^3;     % [m^3] Combined volume of the two hemispherical end caps for the LOX tank
v_caps_fuel_m3 = (4/3)*pi*(r_fuel_tank_m - t_caps_fuel_m)^3; % [m^3] Combined volume of the two hemispherical end caps for the Fuel tank

% Calculate cylindrical section volumes (Eq. 8)
v_cyl_ox_m3 = v_total_ox_tank_m3 - v_caps_ox_m3;     % [m^3] Volume of the LOX tank's cylindrical section
v_cyl_fuel_m3 = v_total_fuel_tank_m3 - v_caps_fuel_m3; % [m^3] Volume of the Fuel tank's cylindrical section

% Calculate cylindrical section lengths (Eq. 9)
l_cyl_ox_m = v_cyl_ox_m3 / (pi * (r_ox_tank_m-t_cyl_ox_m)^2);     % [m] Required length of the LOX tank's cylindrical section
l_cyl_fuel_m = v_cyl_fuel_m3 / (pi * (r_ox_tank_m-t_cyl_ox_m)^2); % [m] Required length of the Fuel tank's cylindrical section

% Calculate empty tank masses (Eq. 13)
m_empty_ox_tank_kg = material_density_ox_kgm3 * (pi * (r_ox_tank_m^2 - (r_ox_tank_m - t_cyl_ox_m)^2)*l_cyl_ox_m + 4*pi/3 * (r_ox_tank_m^3 - (r_ox_tank_m - t_caps_ox_m)^3)); % [kg] Mass of the empty LOX tank
m_empty_fuel_tank_kg = material_density_fuel_kgm3 * (pi * (r_fuel_tank_m^2 - (r_fuel_tank_m - t_cyl_fuel_m)^2)*l_cyl_fuel_m + 4*pi/3 * (r_fuel_tank_m^3 - (r_fuel_tank_m - t_caps_ox_m)^3)); % [kg] Mass of the empty Fuel tank

% --- 2.3 - Pressurant System Sizing --- 

% Calculate moles of pressurant required to displace propellants (Eq. 14, Ideal Gas Law)
% Assumes pressurant gas takes the temperature in the propellant tanks
n_ox_mol = (p_op_ox_tank_pa * v_total_ox_tank_m3) / (r_universal_jmolk * ox_temp_k);       % [mol] Moles of gas to displace oxidizer
n_fuel_mol = (p_op_fuel_tank_pa * v_total_fuel_tank_m3) / (r_universal_jmolk * fuel_temp_k); % [mol] Moles of gas to displace fuel

% Calculate total moles of pressurant (Eq. 15)
n_total_mol = n_ox_mol + n_fuel_mol; % [mol] Total moles of pressurant gas needed

% Calculate total mass of pressurant gas (Eq. 16)
m_pressurant_gas_kg = n_total_mol * pressurant_molar_mass_kgmol; % [kg] Total Mass of pressurant gas needed

% Calculate internal volume of pressurant storage tank (Eq. 17)
v_pressurant_tank_internal_m3 = (n_total_mol * r_universal_jmolk * pressurant_temp_k) / p_storage_pressurant_pa; % [m^3] Internal volume of the high-pressure storage tank

% Calculate internal radius of spherical pressurant tank (Eq. 18)
r_pressurant_tank_internal_m = ((3 * v_pressurant_tank_internal_m3) / (4 * pi))^(1/3); % [m] Internal radius of the pressurant tank

% Calculate wall thickness of spherical pressurant tank (Eq. 19)
p_design_pressurant_pa = p_storage_pressurant_pa * safety_factor; % [Pa] Design pressure for pressurant tank
t_pressurant_tank_m = (p_design_pressurant_pa * r_pressurant_tank_internal_m) / (2 * material_allowable_stress_pressurant_pa * joint_efficiency_pressurant_tank) + corrosion_allowance_m; % [m] Wall thickness of pressurant tank

% Calculate mass of empty pressurant tank (Eq. 20)
v_shell_pressurant_m3 = (4/3) * pi * ((r_pressurant_tank_internal_m + t_liner_m + t_pressurant_tank_m)^3 - (r_pressurant_tank_internal_m + t_liner_m)^3); % [m^3] Volume of the tank material
v_shell_liner_m3 = (4/3) * pi * ((r_pressurant_tank_internal_m + t_liner_m)^3 - (r_pressurant_tank_internal_m )^3); % [m^3] Volume of the tank liner material
m_empty_pressurant_tank_kg = material_density_pressurant_kgm3 * v_shell_pressurant_m3 + material_density_liner_kgm3 * v_shell_liner_m3; % [kg] Mass of the empty pressurant tank

% Calculate pressurant tank outer diameter
d_pressurant_tank_outer_m = 2 * (r_pressurant_tank_internal_m + t_liner_m + t_pressurant_tank_m);     % [m] Outer diameter of pressurant tank

if d_pressurant_tank_outer_m > d_pressurant_tank_allowed_m
    warning('Outer diameter of pressurant COPV exceeds maximum allowed diameter. COPV converted to cylinder.');
end

% Change pressurant tank to cylinder if outer diameter too large
if d_pressurant_tank_outer_m > d_pressurant_tank_allowed_m
    t_pressurant_tank_m = (p_design_pressurant_pa * (r_pressurant_tank_internal_m + t_liner_m)) / (material_allowable_stress_pressurant_pa * joint_efficiency_pressurant_tank) + corrosion_allowance_m; % [m] Wall thickness of pressurant tank
    l_cyl_pressurant_m = (v_pressurant_tank_internal_m3 - 4*pi/3*(r_pressurant_tank_allowed_m - t_liner_m - t_pressurant_tank_m)^3)/(pi * (r_pressurant_tank_allowed_m - t_liner_m - t_pressurant_tank_m)^2); % [m] Length of cylindrical portion of pressurant tank
    m_empty_pressurant_tank_kg = 4*pi/3*(material_density_pressurant_kgm3 * (r_pressurant_tank_allowed_m^3 - (r_pressurant_tank_allowed_m - t_pressurant_tank_m)^3) + material_density_liner_kgm3 * ((r_pressurant_tank_allowed_m - t_pressurant_tank_m)^3 - (r_pressurant_tank_allowed_m - t_pressurant_tank_m - t_liner_m)^3)) + pi*(material_density_pressurant_kgm3*(r_pressurant_tank_allowed_m^2 - (r_pressurant_tank_allowed_m - t_pressurant_tank_m)^2) + material_density_liner_kgm3*((r_pressurant_tank_allowed_m - t_pressurant_tank_m)^2 - (r_pressurant_tank_allowed_m - t_pressurant_tank_m - t_liner_m)^2)) * l_cyl_pressurant_m; % [kg] Mass of the empty pressurant tank 
end

% --- 2.4 - Vehicle Mass Buildup ---
m_total_kg = m_empty_ox_tank_kg + m_empty_fuel_tank_kg + m_ox_kg + m_fuel_kg + m_empty_pressurant_tank_kg + m_pressurant_gas_kg + m_plumbing_kg + m_misc_kg; % [kg] Total vehicle liftoff mass
m_final_kg = m_total_kg - (1 - residual_fraction) * (m_ox_kg + m_fuel_kg); % [kg] Final mass after propellant is consumed

% --- 2.5 - Performance Analysis ---

twr_ratio = (f_thrust_n / m_total_kg) / g_earth_ms2;                                % [unitless] Thrust-to-Weight ratio unitless
delta_v_ms = i_sp_s * g_earth_ms2 * log(m_total_kg / m_final_kg);                   % [m/s] Ideal change in velocity

% 1D Flight Sim
t = 0:0.001:120;                % [s] Time interval and step to analyze
altitude_m = zeros(size(t));      % [m] Pre-allocate altitude array
velocity_ms = zeros(size(t));      % [m/s] Pre-allocate velocity array
acceleration_ms2 = zeros(size(t));  % [m/s^2] Pre-allocate acceleration array
F_drag_N = zeros(size(t));        % [N] Pre-allocate drag force array
thrust_N = zeros(size(t)); thrust_N(t<=t_burn_s) = f_thrust_n;  % [N] Make thrust array (should be zero after burnout)
m_total_sim = m_total_kg - m_dot_total_kgs .* t .* (thrust_N>0); m_total_sim(t>t_burn_s) = m_final_kg; % [kg] Array of vehicle total mass over time
dt = t(2) - t(1);   % [s] time step
A = pi * (D/2)^2;   % [m^2] Vehicle cross-sectional area

% Step through time to calculate altitue, velocity, and acceleration through time
for n = 2:max(size(t))
    velocity_ms(n) = velocity_ms(n-1) + acceleration_ms2(n-1) * dt;                       % [m/s] Velocity array
    altitude_m(n) = altitude_m(n-1) + velocity_ms(n-1) * dt;                           % [m] Altitude array
    [~,~,~,rho_kgm3,~,~] = atmosisa(altitude_m(n));                                    % [kg/m^3] Air density at current altitude
    F_drag_N(n) = - 0.5 * rho_kgm3 * C_D * A * velocity_ms(n) * abs(velocity_ms(n));         % [N] Vehicle drag
    acceleration_ms2(n) = (thrust_N(n) + F_drag_N(n)) / m_total_sim(n) - g_earth_ms2;   % [m/s^2] Vehicle acceleration
% end the simulation if the rocket reaches the ground
if altitude_m(n) < 0
    break
end
end

h_max_m = max(altitude_m);                                          % [m] Estimated apogee
F_drag_max_N = max(F_drag_N);                                       % [N] Maximum drag force
maxq = F_drag_max_N/(C_D * A);                                      % [Pa] Maximum dynamic pressure

endofflight = find(altitude_m(10:end)==0,1,'first'); % [~] Time step when flight ends (used for plotting)

% --- 2.6 - Validation and Error Checks ---
% This section implements the checks defined in the technical plan.

warning_messages = {}; % Initialize a cell array to store warning messages (There might be multiple)

% TWR Check
if twr_ratio < 5 % add this value to the top for the next version
    warning_messages{end+1} = sprintf('TWR is %.2f, which is below the minimum requirement of 5.', twr_ratio);
end

% The Calculations in this section should be moved up for the next version

% Geometric Fit Check
l_ox_tank_total_m = l_cyl_ox_m + d_ox_tank_m; % [m] Length of Lox cylinder + 2*radius of caps
l_fuel_tank_total_m = l_cyl_fuel_m + d_fuel_tank_m; % [m] Length of Lox cylinder + 2*radius of caps
l_pressurant_tank_total_m = d_pressurant_tank_outer_m;
if d_pressurant_tank_outer_m > d_pressurant_tank_allowed_m
    l_pressurant_tank_total_m = l_cyl_pressurant_m + d_pressurant_tank_allowed_m;
end
% NOTE: This is a simplified length sum. A real stackup would be more complex.
l_total_vehicle_m = l_ox_tank_total_m + l_fuel_tank_total_m + l_pressurant_tank_total_m; % Add other lengths later (engine, avionics, etc.)

if l_total_vehicle_m > l_airframe_max_m
    warning_messages{end+1} = sprintf('Estimated vehicle length (%.2f m) exceeds max airframe length (%.2f m).', l_total_vehicle_m, l_airframe_max_m);
end

% Cylindrical Section Length Check
if l_cyl_ox_m < 0
    warning_messages{end+1} = 'Invalid LOX tank geometry: Cylindrical section length is negative. End caps are too large for the required volume.';
end
if l_cyl_fuel_m < 0
    warning_messages{end+1} = 'Invalid Fuel tank geometry: Cylindrical section length is negative. End caps are too large for the required volume.';
end


% Input Sanity Checks (Example)
if ullage_fraction_fuel <= 0 || ullage_fraction_fuel >= 1 || ullage_fraction_ox <= 0 || ullage_fraction_ox >= 1
    warning_messages{end+1} = 'Ullage fraction must be between 0 and 1.';
end
if safety_factor <= 1
    warning_messages{end+1} = 'Safety factor should be greater than 1.';
end


%% 3.0 - OUTPUTS
% This section displays the final calculated values in a clean format.

% This section displays the final calculated values in a clean format.
fprintf('====================================================\n');
fprintf('      Taipan Vehicle Mass & Tank Design Results     \n');
fprintf('====================================================\n\n');

% --- Mass Breakdown ---
fprintf('--- Mass Breakdown ---\n');
fprintf('Oxidizer Mass (LOX):          %8.2f kg\n', m_ox_kg);
fprintf('Fuel Mass (Jet-A):            %8.2f kg\n', m_fuel_kg);
fprintf('Pressurant Gas Mass (N2):     %8.2f kg\n', m_pressurant_gas_kg);
fprintf('------------------------------------------\n');
fprintf('Total Propellant Mass:        %8.2f kg\n', m_ox_kg + m_fuel_kg + m_pressurant_gas_kg);
fprintf('------------------------------------------\n');
fprintf('Empty LOX Tank Mass:          %8.2f kg\n', m_empty_ox_tank_kg);
fprintf('Empty Fuel Tank Mass:         %8.2f kg\n', m_empty_fuel_tank_kg);
fprintf('Empty Pressurant Tank Mass:   %8.2f kg\n', m_empty_pressurant_tank_kg);
fprintf('Plumbing & Misc Mass:         %8.2f kg\n', m_plumbing_kg + m_misc_kg);
fprintf('------------------------------------------\n');
fprintf('Dry Mass (Final):             %8.2f kg\n', m_final_kg);
fprintf('Wet Mass (Liftoff):           %8.2f kg\n\n', m_total_kg);

% --- Tank Geometry ---
fprintf('--- Tank Geometry ---\n');
fprintf('LOX Tank Cyl. Length:         %8.3f m\n', l_cyl_ox_m);
fprintf('Fuel Tank Cyl. Length:        %8.3f m\n', l_cyl_fuel_m);
fprintf('Pressurant Tank Int. Radius:  %8.3f m\n', r_pressurant_tank_internal_m);
fprintf('Pressurant Tank Wall Thick:   %8.4f m (%5.2f mm)\n\n', t_pressurant_tank_m, t_pressurant_tank_m*1000);

% --- Performance Analysis ---
fprintf('--- Performance ---\n');
fprintf('Liftoff Thrust-to-Weight Ratio (TWR): %8.2f\n', twr_ratio);
fprintf('Ideal Delta-V:                %8.2f m/s\n', delta_v_ms);
fprintf('Maximum Drag Force:           %8.2f N\n', F_drag_max_N);
fprintf('Max Q:                        %8.2f Pa\n', maxq);
fprintf('Estimated Apogee:             %8.2f m    (%8.2f ft)\n\n', h_max_m, h_max_m/0.3048);

% --- Warnings & Validation ---
if ~isempty(warning_messages)
    fprintf('!!! --- WARNINGS --- !!!\n');
    for i = 1:length(warning_messages)
        fprintf('  - %s\n', warning_messages{i});
    end
    fprintf('!!! ------------------ !!!\n');
else
    fprintf('--- All design checks passed. ---\n');
end

%% 4.0 - Test Section

% you can use this section temporarly to test that github works for you and the changes you make are actually working



