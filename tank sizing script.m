%{
---------------------------------------------------------------------------
CU Sounding Rocket Lab - Liquid Propulsion Team
Tank Sizing and Vehicle Mass Model

Description:
This script performs a trade study for the design of the propellant and
pressurant tanks for the Taipan liquid rocket engine.

Authors: Hussain Almatruk
Last Updated: 10/17/2025
---------------------------------------------------------------------------
%}

%% 0.0 - SETUP
% Clears workspace, command window, and closes all figures.
clear; clc; close all;

%% 0.1 - CONSTANTS
g_earth_ms2         = 9.80665;       % [m/s^2] Standard gravitational acceleration
r_universal_jmolk   = 8.3145;      % [J/(mol*K)] Universal gas constant

%% 1.0 - INPUTS

% --- Engine & Performance ---
f_thrust_n      = 2891.34405; % [N] Engine thrust
i_sp_s          = 298.9491987;% [s] Specific impulse
t_burn_s        = 10;         % [s] Total burn time
o_f_ratio       = 2.5;        % [unitless] Oxidizer-to-Fuel mass ratio

% --- Propellants & Pressurant ---
ox_density_kgm3   = 1141;       % [kg/m^3] Density of Liquid Oxygen
fuel_density_kgm3 = 820;        % [kg/m^3] Density of Jet-A
pressurant_temp_k = 294;        % [K] Temperature of pressurant gas in its tank
pressurant_molar_mass_gmol = 28.0134; % [g/mol] Molar mass of Nitrogen (N2)

% --- Vehicle Geometry & Materials ---
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

% --- Design Margins & Factors ---
safety_factor   = 1.5;        % [unitless] Safety factor for pressure vessels
% Joint efficiency can differ based on tank material and welding process
joint_efficiency_ox_tank        =; % [unitless] TODO: Define this value 
joint_efficiency_fuel_tank      =; % [unitless] TODO: Define this value 
joint_efficiency_pressurant_tank=; % [unitless] TODO: Define this value 

corrosion_allowance_m = 0.001;  % [m] Extra thickness for material degradation
ullage_fraction       = 0.1;    % [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)

% --- Estimated Masses (Non-Calculated) ---
m_misc_kg     =; % [kg] TODO: Estimate mass of payload, structure, fins, avionics, recovery
m_plumbing_kg =; % [kg] TODO: Estimate mass of valves and plumbing

%% 2.0 - CALCULATIONS
% This section should not be modified unless equations are being updated.
% It should only reference variables from the INPUTS and CONSTANTS sections.

% --- 2.1 - Propellant & Flow Rate Analysis ---
% Placeholder for propellant mass calculations


% --- 2.2 - Tank Sizing & Mass ---
% Placeholder for LOX, Fuel, and Pressurant tank dimension, thickness, and mass calculations



% --- 2.3 - Vehicle Mass Buildup ---
% Placeholder for summing up all masses


% --- 2.4 - Performance Analysis ---
% Placeholder for TWR, Delta-V, and Apogee calculations


%% 3.0 - OUTPUTS
% This section displays the final calculated values in a clean format.

%% 4.0 Test Section

% you can use this section temporarly to test that github works for you and the changes you make are actually working



