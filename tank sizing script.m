%{
---------------------------------------------------------------------------
CU Sounding Rocket Lab - Liquid Propulsion Team
Tank Sizing and Vehicle Mass Model

Description:
This script performs a trade study for the design of the propellant and
pressurant tanks for the Taipan liquid rocket engine.

Authors: Hussain Almatruk
Last Updated: 10/19/2025
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
f_thrust_n      = 2891.34405;  % [N] Engine thrust
i_sp_s          = 298.9491987; % [s] Specific impulse
t_burn_s        = 10;          % [s] Total burn time
o_f_ratio       = 2.5;         % [unitless] Oxidizer-to-Fuel mass ratio

% --- 1.2 Propellants & Pressurant ---
ox_density_kgm3   = 1141;       % [kg/m^3] Density of Liquid Oxygen
fuel_density_kgm3 = 820;        % [kg/m^3] Density of Jet-A
pressurant_temp_k = 294;        % [K] Temperature of pressurant gas in its tank
pressurant_molar_mass_kgmol = 0.0280134; % [kg/mol] Molar mass of Nitrogen (N2)
p_storage_pressurant_pa =;      % [Pa] The pressure the Nitrogen is stored in dedicated tank (MEOP). TODO: Define this value
p_op_ox_tank_pa =; % [Pa] Operating pressure of LOX tank TODO: Define this value
p_op_fuel_tank_pa =; % [Pa] Operating pressure of Fuel tank TODO: Define this value

% --- 1.3 Vehicle Geometry & Materials ---
% Assumption: Both LOX and Fuel tanks are cylinders of the same diameter
d_ox_tank_m     = 0.127;      % [m] Diameter of the LOX tank
d_fuel_tank_m   = 0.127;      % [m] Diameter of the Fuel tank

% Material Properties for LOX Tank 
material_density_ox_kgm3         =; % [kg/m^3] TODO: Define this value
material_allowable_stress_ox_pa  =; % [Pa] TODO: Define this value

% Material Properties for Fuel Tank 
material_density_fuel_kgm3       =; % [kg/m^3] TODO: Define this value
material_allowable_stress_fuel_pa=; % [Pa] TODO: Define this value

% Material Properties for Pressurant Tank 
material_density_pressurant_kgm3       =; % [kg/m^3] TODO: Define this value
material_allowable_stress_pressurant_pa=; % [Pa] TODO: Define this value

% --- 1.4 Design Margins & Factors ---
safety_factor   = 1.5;        % [unitless] Safety factor for pressure vessels
% Joint efficiency can differ based on tank material and welding process
joint_efficiency_ox_tank        =; % [unitless] TODO: Define this value 
joint_efficiency_fuel_tank      =; % [unitless] TODO: Define this value 
joint_efficiency_pressurant_tank=; % [unitless] TODO: Define this value 

corrosion_allowance_m = 0.001;  % [m] Extra thickness for material degradation
ullage_fraction       = 0.1;    % [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)

% --- 1.5 Estimated Masses (Non-Calculated) ---
m_misc_kg     =; % [kg] TODO: Estimate mass of payload, structure, fins, avionics, recovery
m_plumbing_kg =; % [kg] TODO: Estimate mass of valves and plumbing

% --- 1.6 Design Constraints ---

l_airframe_max_m =; %[m] Maximum allowable vehicle length TODO: Define this value 

%% 2.0 - CALCULATIONS
% This section should not be modified unless equations are being updated or the initial code is done and is now being itirated on to optimize values

% --- 2.1 - Propellant & Flow Rate Analysis ---

% Calculate Prop Mass Flow Rates (Eq. 1,2,3)
m_dot_total_kgs = f_thrust_n/(i_sp_s*g_earth_ms2);          % [kg/s] Total propellant mass flow rate
m_dot_ox_kgs = m_dot_total_kgs*(o_f_ratio/(1+o_f_ratio));   % [kg/s] Oxidizer mass flow rate 
m_dot_fuel_kgs = m_dot_total_kgs/(1+o_f_ratio);             % [kg/s] Fuel mass flow rate 

% Calculate Prop Masses (Eq. 4)
m_ox_kg = m_dot_ox_kgs*t_burn_s;        % [kg] Total mass of oxidizer
m_fuel_kg = m_dot_fuel_kgs *t_burn_s;   % [kg] Total mass of fuel

% --- 2.2 - Tank Sizing & Mass ---

% Calculate Volume of Props (Eq. 5)
v_ox_m3 = m_ox_kg/ox_density_kgm3;          % [m^3] Volume of liquid oxidizer
v_fuel_m3 = m_fuel_kg/fuel_density_kgm3;    % [m^3] Volume of liquid fuel

% Calculate Internal Volume of Prop Tanks (Volume of prop in addition to ullage) (Eq. 6)
v_total_ox_tank_m3 = v_ox_m3/(1-ullage_fraction);       % [m^3] Total internal volume of LOX tank
v_total_fuel_tank_m3 = v_fuel_m3/(1-ullage_fraction);   % [m^3] Total internal volume of Fuel tank

% Calculate tank radii (Basic Geometry)
r_ox_tank_m = d_ox_tank_m/2;     % [m] Radius of cylindrical ox tank
r_fuel_tank_m = d_fuel_tank_m/2; % [m] Radius of cylindrical fuel tank

% Calculate end cap volumes (Eq. 7)
v_caps_ox_m3 = (4/3)*pi*r_ox_tank_m^3;     % [m^3] Combined volume of the two hemispherical end caps for the LOX tank
v_caps_fuel_m3 = (4/3)*pi*r_fuel_tank_m^3; % [m^3] Combined volume of the two hemispherical end caps for the Fuel tank

% Calculate cylindrical section volumes (Eq. 8)
v_cyl_ox_m3 = v_total_ox_tank_m3 - v_caps_ox_m3;     % [m^3] Volume of the LOX tank's cylindrical section
v_cyl_fuel_m3 = v_total_fuel_tank_m3 - v_caps_fuel_m3; % [m^3] Volume of the Fuel tank's cylindrical section

% Calculate cylindrical section lengths (Eq. 9)
l_cyl_ox_m = v_cyl_ox_m3 / (pi * r_ox_tank_m^2);     % [m] Required length of the LOX tank's cylindrical section
l_cyl_fuel_m = v_cyl_fuel_m3 / (pi * r_fuel_tank_m^2); % [m] Required length of the Fuel tank's cylindrical section

% Calculate tank design pressures (Eq. 10)
p_design_ox_pa = p_op_ox_tank_pa * safety_factor;     % [Pa] Design pressure for LOX tank
p_design_fuel_pa = p_op_fuel_tank_pa * safety_factor; % [Pa] Design pressure for Fuel tank

% Calculate cylinder wall thicknesses (Eq. 11, ASME Hoop Stress)
t_cyl_ox_m = (p_design_ox_pa * d_ox_tank_m) / (2 * material_allowable_stress_ox_pa * joint_efficiency_ox_tank) + corrosion_allowance_m; % [m] LOX tank cylinder thickness
t_cyl_fuel_m = (p_design_fuel_pa * d_fuel_tank_m) / (2 * material_allowable_stress_fuel_pa * joint_efficiency_fuel_tank) + corrosion_allowance_m; % [m] Fuel tank cylinder thickness

% Calculate end cap wall thicknesses (Eq. 12, Spherical Shell Stress)
t_caps_ox_m = (p_design_ox_pa * d_ox_tank_m) / (4 * material_allowable_stress_ox_pa * joint_efficiency_ox_tank) + corrosion_allowance_m; % [m] LOX tank end cap thickness
t_caps_fuel_m = (p_design_fuel_pa * d_fuel_tank_m) / (4 * material_allowable_stress_fuel_pa * joint_efficiency_fuel_tank) + corrosion_allowance_m; % [m] Fuel tank end cap thickness

% Calculate empty tank masses (Eq. 13)
m_empty_ox_tank_kg = material_density_ox_kgm3 * (pi * d_ox_tank_m * l_cyl_ox_m * t_cyl_ox_m + pi * d_ox_tank_m^2 * t_caps_ox_m); % [kg] Mass of the empty LOX tank
m_empty_fuel_tank_kg = material_density_fuel_kgm3 * (pi * d_fuel_tank_m * l_cyl_fuel_m * t_cyl_fuel_m + pi * d_fuel_tank_m^2 * t_caps_fuel_m); % [kg] Mass of the empty Fuel tank

% --- 2.3 - Pressurant System Sizing --- 

% Calculate moles of pressurant required to displace propellants (Eq. 14, Ideal Gas Law)
% Assumes pressurant gas cools to its storage temperature in the propellant tanks
n_ox_mol = (p_op_ox_tank_pa * v_ox_m3) / (r_universal_jmolk * pressurant_temp_k);       % [mol] Moles of gas to displace oxidizer
n_fuel_mol = (p_op_fuel_tank_pa * v_fuel_m3) / (r_universal_jmolk * pressurant_temp_k); % [mol] Moles of gas to displace fuel

% Calculate total moles of pressurant (Eq. 15)
n_total_mol = n_ox_mol + n_fuel_mol; % [mol] Total moles of pressurant gas needed

% Calculate total mass of pressurant gas (Eq. 16)
m_pressurant_gas_kg = n_total_mol * (pressurant_molar_mass_kgmol); % [kg] Total Mass of pressurant gas needed

% Calculate internal volume of pressurant storage tank (Eq. 17)
v_pressurant_tank_internal_m3 = (n_total_mol * r_universal_jmolk * pressurant_temp_k) / p_storage_pressurant_pa; % [m^3] Internal volume of the high-pressure storage tank

% Calculate internal radius of spherical pressurant tank (Eq. 18)
r_pressurant_tank_internal_m = ((3 * v_pressurant_tank_internal_m3) / (4 * pi))^(1/3); % [m] Internal radius of the pressurant tank

% Calculate wall thickness of spherical pressurant tank (Eq. 19)
p_design_pressurant_pa = p_storage_pressurant_pa * safety_factor; % [Pa] Design pressure for pressurant tank
t_pressurant_tank_m = (p_design_pressurant_pa * r_pressurant_tank_internal_m) / (2 * material_allowable_stress_pressurant_pa * joint_efficiency_pressurant_tank) + corrosion_allowance_m; % [m] Wall thickness of pressurant tank

% Calculate mass of empty pressurant tank (Eq. 20)
v_shell_pressurant_m3 = (4/3) * pi * ((r_pressurant_tank_internal_m + t_pressurant_tank_m)^3 - r_pressurant_tank_internal_m^3); % [m^3] Volume of the tank material
m_empty_pressurant_tank_kg = material_density_pressurant_kgm3 * v_shell_pressurant_m3; % [kg] Mass of the empty pressurant tank

% --- 2.4 - Vehicle Mass Buildup ---
% Placeholder for summing up all masses





% --- 2.5 - Performance Analysis ---
% Placeholder for TWR, Delta-V, and Apogee calculations





%% 3.0 - OUTPUTS
% This section displays the final calculated values in a clean format.


%% 4.0 - Test Section

% you can use this section temporarly to test that github works for you and the changes you make are actually working



