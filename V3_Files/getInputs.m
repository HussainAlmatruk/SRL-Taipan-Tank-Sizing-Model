function P = getInputs(C)
%{
---------------------------------------------------------------------------
getInputs
  Returns a struct containing all vehicle design inputs and parameters.
  This is the PRIMARY FILE TO EDIT when running trade studies or updating
  design values. Physical constants and unit conversions are not defined
  here — they live in getConstants.m and are passed in as the struct C.

  To use a value from this struct in a calculation:
      P.f_thrust_n, P.i_sp_s, P.d_tank_outer_m, etc.

Input:
  C - [struct] Constants struct from getConstants(), needed for unit
               conversions within input definitions (e.g. lbf to N)

Output:
  P - [struct] Struct containing all vehicle and design input parameters
---------------------------------------------------------------------------
%}

% --- 1.1 Engine & Performance ---
P.f_thrust_n    = 650 * C.LBF_TO_N;    % [N]        Engine thrust (Updated)
P.i_sp_s        = 291.29;              % [s]        Specific impulse (Derived from new Thrust and Mass Flow)
P.t_burn_s      = 9;                   % [s]        Total burn time
P.o_f_ratio     = 2.3;                 % [unitless] Oxidizer-to-Fuel mass ratio

% --- 1.2 Propellants & Pressurant ---
P.ox_density_kgm3             = 1141;          % [kg/m^3]  Density of Liquid Oxygen
P.fuel_density_kgm3           = 820;           % [kg/m^3]  Density of Jet-A
P.pressurant_molar_mass_kgmol = 0.0280134;     % [kg/mol]  Molar mass of Nitrogen (N2)
P.residual_fraction           = 0.1;           % %%% TEMPORARY VALUE %%% [~] Percentage of fuel left over in tank after burnout

P.p_op_ox_tank_pa         = 1000 * C.PSI_TO_PA;    % [Pa] Operating pressure of LOX tank (Updated from 8000 kPa)
P.p_op_fuel_tank_pa       = 1000 * C.PSI_TO_PA;    % [Pa] Operating pressure of Fuel tank (Updated from 5050 kPa)
P.p_storage_pressurant_pa = 3000 * C.PSI_TO_PA;    % [Pa] The pressure the Nitrogen is stored in dedicated tank (MEOP). TODO: Define this value
P.ox_temp_k               = 90;                    % %%% TEMPORARY VALUE %%% [K] Temperature of LOX in tank TODO: Define this value
P.fuel_temp_k             = 294;                   % %%% TEMPORARY VALUE %%% [K] Temperature of Jet-A in tank (Ambient probably)
P.pressurant_temp_k       = 294;                   % %%% TEMPORARY VALUE %%% [K] Temperature of pressurant gas in its tank

% --- 1.3 Tank Construction Properties ---
P.prescribe_t = 0;      % [unitless] 1 - Prescribe tank thicknesses, 0 - calculate optimal tank thickness

% Prescribed wall thicknesses — only applied inside calcMassAndSizing if prescribe_t = 1
P.t_cyl_ox_prescribed_m     = 0.120*2.54/100;  % [m] LOX tank cylinder thickness
P.t_cyl_fuel_prescribed_m   = 0.120*2.54/100;  % [m] Fuel tank cylinder thickness
P.t_pressurant_prescribed_m = 0.337*2.54/100;  % [m] Pressurant tank wall thickness

% --- 1.3b COTS Pressurant Tank ---
% Toggle between a purchased off-the-shelf pressure vessel (true) and the
% calculated composite overwrapped pressure vessel / COPV math (false).
% When use_cots_tank = true, the COPV material properties in Section 1.4
% are ignored and the placeholder specs below are used instead.
% TODO: Replace placeholder values with actual COTS tank spec sheet values.
P.use_cots_tank = true;                     % [unitless] true - use COTS tank specs, false - calculate COPV

P.cots_tank_mass_kg           = 5;          % %%% PLACEHOLDER %%%  [kg]  Empty mass of the COTS tank TODO: Update from spec sheet
P.cots_tank_volume_m3         = 0.010;      % %%% PLACEHOLDER %%%  [m^3] Internal volume of the COTS tank TODO: Update from spec sheet
P.cots_tank_outer_diameter_m  = 0.1143;     % %%% PLACEHOLDER %%%  [m]   Outer diameter of the COTS tank (4.5 in placeholder) TODO: Update from spec sheet
P.cots_tank_length_m          = 0.400;      % %%% PLACEHOLDER %%%  [m]   Total length of the COTS tank TODO: Update from spec sheet

% --- 1.4 Vehicle Geometry & Materials ---
% Assumption: Both LOX and Fuel tanks are cylinders of the same diameter

% Tank outer diameter is the primary input
P.d_tank_outer_m    = 4.5 * C.IN_TO_M; % [m] Tank outer diameter
P.wall_thickness_m  = 0*2.54/100;      % [m] Vehicle airframe wall thickness
P.inner_clearance_m = 0*2.54/100;      % [m] Clearance between tank and vehicle inner wall

% Material Properties for LOX Tank
P.material_density_ox_kgm3        = 2700;   % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
P.material_allowable_stress_ox_pa = 2.76e8; % %%% TEMPORARY VALUE %%%  [Pa]     TODO: Define this value

% Material Properties for Fuel Tank
P.material_density_fuel_kgm3          = 2700;   % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
P.material_allowable_stress_fuel_pa   = 2.76e8; % %%% TEMPORARY VALUE %%%  [Pa]     TODO: Define this value

% Material Properties for Pressurant Tank
P.material_density_liner_kgm3             = 2700;   % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
P.t_liner_m                               = 0;      % %%% TEMPORARY VALUE %%%  [m]      Thickness of COPV liner
P.material_density_pressurant_kgm3        = 2700;   % %%% TEMPORARY VALUE %%%  [kg/m^3] TODO: Define this value
P.material_allowable_stress_pressurant_pa = 2.76e8; % %%% TEMPORARY VALUE %%%  [Pa]     TODO: Define this value

% --- 1.5 Design Margins & Factors ---
P.safety_factor = 1.5;  % [unitless] Safety factor for pressure vessels

% Joint efficiency can differ based on tank material and welding process
P.joint_efficiency_ox_tank         = 0.8; % %%% TEMPORARY VALUE %%%  [unitless] TODO: Define this value
P.joint_efficiency_fuel_tank       = 0.8; % %%% TEMPORARY VALUE %%%  [unitless] TODO: Define this value
P.joint_efficiency_pressurant_tank = 0.8; % %%% TEMPORARY VALUE %%%  [unitless] TODO: Define this value

P.corrosion_allowance_m = 0.0002;  % %%% TEMPORARY VALUE %%%  [m]        Extra thickness for material degradation
P.ullage_fraction_ox    = 0.2;     % %%% TEMPORARY VALUE %%%  [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)
P.ullage_fraction_fuel  = 0.1;     % %%% TEMPORARY VALUE %%%  [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)

% --- 1.6 Mass Assumptions ---
% NOTE: The combined mass of miscellaneous items (payload, structure, fins,
% avionics, recovery) and plumbing (valves, lines) is assumed to equal the
% total dry tank mass (LOX tank + Fuel tank + Pressurant tank, all empty).
% m_misc_and_plumbing_kg is computed in sizeTanks after tank masses are
% known, and is available as an output of that function.

% --- 1.7 Design Constraints ---
P.l_airframe_max_m  = 3.048;   % [m]     Maximum allowable vehicle length
P.twr_minimum_ratio = 5;       % [ratio] Thrust to Weight ratio

% --- 1.8 Vehicle Aerodynamic Properties ---
rasaero_data = readtable("4in_RASAero_CD_data.CSV");                                                                                    % %%% TEMPORARY VALUES %%% [~] Load in vehicle aerodynamic data from RASAero
P.rasaero_data_0_AoA = rasaero_data(rasaero_data.Alpha == 0,:); P.rasaero_data_0_AoA = P.rasaero_data_0_AoA(1:end-1,:);  % [~] Seperate out data for 0 degree angle of attack (also trim to be same size as 4 deg array)
P.rasaero_data_2_AoA = rasaero_data(rasaero_data.Alpha == 2,:); P.rasaero_data_2_AoA = P.rasaero_data_2_AoA(1:end-1,:);  % [~] Seperate out data for 2 degree angle of attack (also trim to be same size as 4 deg array)
P.rasaero_data_4_AoA = rasaero_data(rasaero_data.Alpha == 4,:);                                                            % [~] Seperate out data for 4 degree angle of attack

end
