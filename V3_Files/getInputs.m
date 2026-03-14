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
P.residual_fraction           = 0.05;           % %%% TEMPORARY VALUE %%% [~] Percentage of fuel left over in tank after burnout

P.p_op_ox_tank_pa         = 1000 * C.PSI_TO_PA;    % [Pa] Operating pressure of LOX tank
P.p_op_fuel_tank_pa       = 1000 * C.PSI_TO_PA;    % [Pa] Operating pressure of Fuel tank
P.p_storage_pressurant_pa = 3000 * C.PSI_TO_PA;    % [Pa] The pressure the Nitrogen is stored in dedicated tank (MEOP)
P.ox_temp_k               = 90;                    % %%% TEMPORARY VALUE %%% [K] Temperature of LOX in tank TODO: Define this value
P.fuel_temp_k             = 294;                   % %%% TEMPORARY VALUE %%% [K] Temperature of Jet-A in tank (Ambient probably)

% Pressurant fill temperature — the temperature at which the tank is
% filled and stored pre-launch.  Used to compute how much gas mass
% fits in the storage tank and for the required storage volume check.

P.pressurant_fill_temp_k  = 294;                   % [K] Temperature of pressurant gas at time of filling (always ambient)

% Legacy alias — kept so that any old code referencing
% P.pressurant_temp_k still works.  Points to the fill temperature.
P.pressurant_temp_k       = P.pressurant_fill_temp_k; % [K] Alias for backwards compatibility

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
% are ignored and the specs below are used instead.
%
% Current tank: Philips Respironics UltraFill ME36, 3000 PSI Oxygen Cylinder
% This is an E-size aluminum 6061 cylinder rated for 3000 PSI (DOT-3AL).
% Specs from spec sheet and standard E-cylinder dimensions:
%   Height:         29.1 in  (0.7391 m) — includes valve
%   Outer Diameter: 4.38 in  (0.1113 m) — standard E-cylinder OD
%   Weight:         9.3 lbs  (4.22 kg)  — includes valve
%   Internal Volume: ~4.5 L  (0.0045 m^3) — lowest estimate
%
% Wall thickness and inner diameter are ESTIMATED from DOT-3AL minimum
% wall stress calculations for 3000 PSI service on Al 6061-T6.
% Standard E-size at 2015 PSI has 7.1 mm wall thickness.  At 3000 PSI
% the DOT formula gives ~8.0 mm.  Confirm with physical measurement.
% See notes below for derivation.
P.use_cots_tank = true;                      % [unitless] true - use COTS tank specs (ME36), false - calculate COPV

P.cots_tank_mass_kg           = 9.3 * C.LBM_TO_KG;    % [kg]  Empty mass of the ME36 tank (9.3 lbs from spec sheet, includes valve)
P.cots_tank_volume_m3         = 0.0045;                % [m^3] Internal volume of the ME36 (4.5 L lowest estimate) TODO: Confirm with physical measurement
P.cots_tank_outer_diameter_m  = 4.38 * C.IN_TO_M;     % [m]   Outer diameter (4.38 in, standard E-cylinder OD)
P.cots_tank_length_m          = 29.1 * C.IN_TO_M;     % [m]   Total height of the ME36 from spec sheet (29.1 in, includes valve)

% --- ESTIMATED dimensions (from DOT-3AL wall stress calculation) ---
% For a DOT-3AL cylinder at 3000 PSI:
%   Test pressure = (5/3) * 3000 = 5000 PSI
%   Al 6061-T6 UTS ~ 45,000 PSI, allowable wall stress = (2/3)*UTS = 30,000 PSI
%   Solving the DOT thick-wall formula: S = P*(1.3*D^2 + 0.4*d^2)/(D^2 - d^2)
%   gives inner diameter d ~ 3.75 in, wall thickness ~ 0.315 in ~ 8.0 mm.
% These are estimates, we need to confirm with physical measurement or manufacturer data.
P.cots_tank_wall_thickness_m  = 0.0080;                                                       % %%% ESTIMATED %%% [m] Wall thickness (~8.0 mm for 3000 PSI Al 6061)
P.cots_tank_inner_diameter_m  = P.cots_tank_outer_diameter_m - 2 * P.cots_tank_wall_thickness_m; % [m] Inner diameter (computed from OD and wall thickness)

% --- Internal surface area estimate ---
% Assumes hemispherical end caps (common for small high-pressure cylinders).
% Computed from inner diameter and internal volume.
%   Two hemispheres = one sphere: V_hemi = (4/3)*pi*ri^3
%   Cylinder length from remaining volume: L_cyl = (V_total - V_hemi) / (pi*ri^2)
%   Total area = 4*pi*ri^2 (sphere) + 2*pi*ri*L_cyl (cylinder wall)
ri_cots_m       = P.cots_tank_inner_diameter_m / 2;                                   % [m]   Inner radius
v_hemi_cots_m3  = (4/3) * pi * ri_cots_m^3;                                           % [m^3] Volume of two hemispherical end caps
l_cyl_cots_m    = (P.cots_tank_volume_m3 - v_hemi_cots_m3) / (pi * ri_cots_m^2);    % [m]   Length of cylindrical section (internal)
P.cots_tank_internal_area_m2 = 4*pi*ri_cots_m^2 + 2*pi*ri_cots_m*l_cyl_cots_m;      % %%% ESTIMATED %%% [m^2] Total internal surface area

% --- Wall thermal mass ---
% For the blowdown model, we need the mass and specific heat of the tank
% wall to track how much thermal energy the wall can supply to the cooling
% gas.  The total spec-sheet mass includes the valve (~0.3 kg).  We
% subtract the valve mass to get the wall mass alone.
P.cots_tank_valve_mass_kg = 0.3;                                                       % %%% ESTIMATED %%% [kg] Approximate mass of the cylinder valve
P.cots_tank_wall_mass_kg  = P.cots_tank_mass_kg - P.cots_tank_valve_mass_kg;          % [kg]  Tank wall mass (total - valve)

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
P.ullage_fraction_ox    = 0.05;     % %%% TEMPORARY VALUE %%%  [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)
P.ullage_fraction_fuel  = 0.05;     % %%% TEMPORARY VALUE %%%  [unitless] Percent of empty volume in tanks (e.g., 0.1 for 10%)

% --- 1.6 Mass Assumptions ---
% NOTE: The combined mass of miscellaneous items (payload, structure, fins,
% avionics, recovery) and plumbing (valves, lines) is assumed to equal the
% total dry tank mass (LOX tank + Fuel tank + Pressurant tank, all empty).
% m_misc_and_plumbing_kg is computed in sizeTanks after tank masses are
% known, and is available as an output of that function.

% --- 1.7 Design Constraints ---
P.l_airframe_max_m  = 3.048;   % [m]     Maximum allowable vehicle length
P.twr_minimum_ratio = 5;       % [ratio] Thrust to Weight ratio

% --- 1.8 Pressurant Thermal Analysis Parameters ---
% These inputs control the transient ullage gas thermal model

% Toggle: true = use transient thermal model (more accurate, lower mass)
%         false = use isothermal model only (conservative worst-case)
P.use_thermal_model = true; % [unitless] true - use transient thermal model, false - isothermal only

% Convective heat transfer coefficient [W/(m^2*K)]
% This is the most important and most uncertain parameter in the analysis.
% It controls how fast the ullage gas loses heat to the propellant surface.
% Recommended values:
%   Conservative (quiescent):       h = 50   W/(m^2*K) — use for structural sizing
%   Moderate (natural convection):  h = 100  W/(m^2*K) — use for mass budgets
%   Aggressive (forced convection): h = 200  W/(m^2*K) — performance upper bound
P.h_convection_ox_wm2k   = 100;   % [W/(m^2*K)] Heat transfer coefficient for LOX tank (big effect)
P.h_convection_fuel_wm2k = 100;    % [W/(m^2*K)] Heat transfer coefficient for Fuel tank (not that big of an effect, there isn't as much of a temp. diff.)

% Nitrogen gas specific heat at constant pressure
% Approximately constant over the 90-294 K range at moderate pressures.
P.cp_n2_jkgk = 1040;              % [J/(kg*K)] Specific heat of N2 at constant pressure

% Timestep for the thermal simulation (forward Euler method)
% 0.01 s gives 900 steps over a 9 s burn, which is plenty for convergence.
P.dt_thermal_s = 0.01;            % [s] Timestep for the transient thermal model

% --- 1.8b Pressurant Blowdown Model Parameters ---
% The blowdown model tracks the pressurant tank gas temperature as it
% cools due to expansion during the burn.  This cooling means the gas
% entering the propellant tanks gets progressively colder, which increases
% the total pressurant mass requirement compared to the V1 thermal model.
%
% Toggle: true = coupled blowdown + thermal model (most accurate)
%         false = V1 thermal model only (constant inlet temp, optimistic)
P.use_blowdown_model = true; % [unitless] true - coupled blowdown model, false - V1 constant inlet temp

% Nitrogen gas properties (additional, needed for blowdown energy balance)
P.cv_n2_jkgk = 743;    % [J/(kg*K)]  Specific heat of N2 at constant volume
P.gamma_n2   = 1.4;    % [unitless]  Specific heat ratio for N2, gamma = cp/cv

% Internal convective heat transfer coefficient for the pressurant tank wall
% This controls how fast the aluminum tank wall can heat the expanding gas.

P.h_wall_pressurant_wm2k = 50;    % [W/(m^2*K)] Internal convective HTC for pressurant tank (conservative)

% Aluminum 6061 specific heat (for tank wall thermal mass)
P.c_wall_pressurant_jkgk = 896;   % [J/(kg*K)]  Specific heat of Al 6061

% --- 1.9 Vehicle Aerodynamic Properties ---
rasaero_data = readtable("4in_RASAero_CD_data.CSV");                                                                                    % %%% TEMPORARY VALUES %%% [~] Load in vehicle aerodynamic data from RASAero
P.rasaero_data_0_AoA = rasaero_data(rasaero_data.Alpha == 0,:); P.rasaero_data_0_AoA = P.rasaero_data_0_AoA(1:end-1,:);  % [~] Seperate out data for 0 degree angle of attack (also trim to be same size as 4 deg array)
P.rasaero_data_2_AoA = rasaero_data(rasaero_data.Alpha == 2,:); P.rasaero_data_2_AoA = P.rasaero_data_2_AoA(1:end-1,:);  % [~] Seperate out data for 2 degree angle of attack (also trim to be same size as 4 deg array)
P.rasaero_data_4_AoA = rasaero_data(rasaero_data.Alpha == 4,:);                                                            % [~] Seperate out data for 4 degree angle of attack

end