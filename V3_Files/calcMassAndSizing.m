function [m_dot_total_kgs, m_dot_ox_kgs, m_dot_fuel_kgs,                          ...
          m_ox_kg, m_fuel_kg,                                                       ...
          d_vehicle_outer_m, d_ox_tank_m, d_fuel_tank_m,                           ...
          t_cyl_ox_m, t_cyl_fuel_m,                                                 ...
          v_total_ox_tank_m3, v_total_fuel_tank_m3,                                 ...
          l_cyl_ox_m, l_cyl_fuel_m,                                                 ...
          m_empty_ox_tank_kg, m_empty_fuel_tank_kg,                                 ...
          m_pressurant_gas_kg, d_pressurant_tank_outer_m,                           ...
          r_pressurant_outer_m, r_pressurant_internal_m,                            ...
          t_pressurant_tank_m, l_cyl_pressurant_m,                                  ...
          m_empty_pressurant_tank_kg, m_misc_and_plumbing_kg,                       ...
          l_ox_tank_total_m, l_fuel_tank_total_m, l_pressurant_tank_total_m,       ...
          l_total_vehicle_m,                                                         ...
          m_total_kg, m_final_kg, twr_ratio, delta_v_ms,                            ...
          v_ox_m3, v_fuel_m3, v_pressurant_tank_m3,                                 ...
          m_pressurant_gas_isothermal_kg, v_storage_isothermal_L,                   ...
          m_pressurant_gas_cold_kg, v_storage_cold_L,                               ...
          m_pressurant_gas_warm_kg, v_storage_warm_L,                               ...
          m_pressurant_gas_thermal_kg, v_storage_thermal_L,                         ...
          T_ox_history_K, T_fuel_history_K, t_history_s] = calcMassAndSizing(P, C)
%{
---------------------------------------------------------------------------
calcMassAndSizing
  Single function that performs all mass and sizing calculations for the
  vehicle in the correct order. Internally organised into four stages:

    1. Propellant & Flow Rate Analysis
       Calculates mass flow rates and total propellant masses from engine
       performance inputs.

    2. Tank Sizing & Mass (LOX, Fuel, and Pressurant/COPV)
       Sizes all three tanks: wall thicknesses, internal volumes,
       cylindrical section lengths, and empty tank masses.

    3. Mass Assumption (Section 1.6 of getInputs.m)
       Computes m_misc_and_plumbing_kg once tank dry masses are known.

    4. Vehicle Geometry Stackup
       Computes total axial length of each tank and the simplified vehicle
       length. For 2:1 ellipsoidal caps, each cap height = r/2, so two
       caps together add d/2 to the cylinder length.

    5. Vehicle Mass Buildup & Static Performance
       Builds wet and dry masses and computes TWR and ideal delta-V.

Inputs:
  P - [struct] Design inputs from getInputs()
  C - [struct] Physical constants from getConstants()

Outputs:
  m_dot_total_kgs            - [kg/s]    Total propellant mass flow rate
  m_dot_ox_kgs               - [kg/s]    Oxidizer mass flow rate
  m_dot_fuel_kgs             - [kg/s]    Fuel mass flow rate
  m_ox_kg                    - [kg]      Total mass of oxidizer loaded into tank
  m_fuel_kg                  - [kg]      Total mass of fuel loaded into tank
  d_vehicle_outer_m          - [m]       Outer diameter of the vehicle airframe
  d_ox_tank_m                - [m]       Outer diameter of the LOX tank
  d_fuel_tank_m              - [m]       Outer diameter of the Fuel tank
  t_cyl_ox_m                 - [m]       LOX tank cylinder wall thickness
  t_cyl_fuel_m               - [m]       Fuel tank cylinder wall thickness
  v_total_ox_tank_m3         - [m^3]     Total internal volume of LOX tank (propellant + ullage)
  v_total_fuel_tank_m3       - [m^3]     Total internal volume of Fuel tank (propellant + ullage)
  l_cyl_ox_m                 - [m]       Length of LOX tank cylindrical section
  l_cyl_fuel_m               - [m]       Length of Fuel tank cylindrical section
  m_empty_ox_tank_kg         - [kg]      Mass of empty LOX tank
  m_empty_fuel_tank_kg       - [kg]      Mass of empty Fuel tank
  m_pressurant_gas_kg        - [kg]      Total mass of pressurant gas needed
  d_pressurant_tank_outer_m  - [m]       Final outer diameter of the pressurant tank
  r_pressurant_outer_m       - [m]       Outer radius of the pressurant tank
  r_pressurant_internal_m    - [m]       Internal radius of the pressurant tank (inside liner)
  t_pressurant_tank_m        - [m]       Pressurant tank wall thickness
  l_cyl_pressurant_m         - [m]       Length of pressurant tank cylindrical section
  m_empty_pressurant_tank_kg - [kg]      Total empty mass of the COPV (shell + liner)
  m_misc_and_plumbing_kg     - [kg]      Combined misc and plumbing mass (= total dry tank mass, see Section 1.6 of getInputs.m)
  l_ox_tank_total_m          - [m]       Total length of the LOX tank
  l_fuel_tank_total_m        - [m]       Total length of the Fuel tank
  l_pressurant_tank_total_m  - [m]       Total length of the Pressurant tank
  l_total_vehicle_m          - [m]       Simplified total vehicle length (tanks only). TODO: Add engine, avionics, recovery, etc.
  m_total_kg                 - [kg]      Total vehicle liftoff (wet) mass
  m_final_kg                 - [kg]      Vehicle mass after propellant is consumed (dry + residuals)
  twr_ratio                  - [unitless] Liftoff Thrust-to-Weight ratio
  delta_v_ms                 - [m/s]     Ideal Tsiolkovsky delta-V
  v_ox_m3                    - [m^3]     Volume of liquid oxidizer (propellant only, no ullage)
  v_fuel_m3                  - [m^3]     Volume of liquid fuel (propellant only, no ullage)
  v_pressurant_tank_m3       - [m^3]     Internal volume of the pressurant tank
  m_pressurant_gas_isothermal_kg - [kg]  Pressurant mass from the OLD isothermal model (for comparison)
  v_storage_isothermal_L     - [L]       Required storage volume from isothermal model (for comparison)
  m_pressurant_gas_cold_kg   - [kg]      COLD BOUND: pressurant mass if gas stays at LOX temp (90 K) — absolute max
  v_storage_cold_L           - [L]       COLD BOUND: required storage volume
  m_pressurant_gas_warm_kg   - [kg]      WARM BOUND: pressurant mass if gas stays at inlet temp (294 K) — absolute min
  v_storage_warm_L           - [L]       WARM BOUND: required storage volume
  m_pressurant_gas_thermal_kg - [kg]     Pressurant mass from the transient thermal model
  v_storage_thermal_L        - [L]       Required storage volume from thermal model
  T_ox_history_K             - [K]       LOX ullage gas temperature history (for plotting)
  T_fuel_history_K           - [K]       Fuel ullage gas temperature history (for plotting)
  t_history_s                - [s]       Time array for temperature histories (for plotting)
---------------------------------------------------------------------------
%}

%% --- 1. PROPELLANT & FLOW RATE ANALYSIS ---

% Calculate Prop Mass Flow Rates (Eq. 1,2,3)
m_dot_total_kgs = P.f_thrust_n/(P.i_sp_s*C.g_earth_ms2);              % [kg/s] Total propellant mass flow rate
m_dot_ox_kgs    = m_dot_total_kgs*(P.o_f_ratio/(1+P.o_f_ratio));       % [kg/s] Oxidizer mass flow rate
m_dot_fuel_kgs  = m_dot_total_kgs/(1+P.o_f_ratio);                     % [kg/s] Fuel mass flow rate

% Calculate Prop Masses (Eq. 4)
m_ox_kg   = m_dot_ox_kgs   * P.t_burn_s/(1-P.residual_fraction);       % [kg] Total mass of oxidizer
m_fuel_kg = m_dot_fuel_kgs * P.t_burn_s/(1-P.residual_fraction);       % [kg] Total mass of fuel

%% --- 2. TANK SIZING & MASS ---

% --- LOX & Fuel Tanks ---

% Calculate vehicle inner diameter and the available diameter for tanks
d_vehicle_outer_m = P.d_tank_outer_m + 2*P.wall_thickness_m + 2*P.inner_clearance_m; % [m] Outer diameter of the vehicle airframe

% Assign the calculated diameter to all tanks
d_ox_tank_m   = P.d_tank_outer_m;  % [m] Outer diameter of the LOX tank
d_fuel_tank_m = P.d_tank_outer_m;  % [m] Outer diameter of the Fuel tank

% The rest of your calculations that depend on these diameters would follow...

% Calculate Volume of Props (Eq. 5)
v_ox_m3   = m_ox_kg/P.ox_density_kgm3;     % [m^3] Volume of liquid oxidizer
v_fuel_m3 = m_fuel_kg/P.fuel_density_kgm3; % [m^3] Volume of liquid fuel

% Calculate Internal Volume of Prop Tanks (Volume of prop in addition to ullage) (Eq. 6)
v_total_ox_tank_m3   = v_ox_m3/(1-P.ullage_fraction_ox);     % [m^3] Total internal volume of LOX tank
v_total_fuel_tank_m3 = v_fuel_m3/(1-P.ullage_fraction_fuel); % [m^3] Total internal volume of Fuel tank

% Calculate tank outer radii (Basic Geometry)
r_ox_tank_m   = d_ox_tank_m/2;   % [m] Outer radius of cylindrical ox tank
r_fuel_tank_m = d_fuel_tank_m/2; % [m] Outer radius of cylindrical fuel tank

% Calculate tank design pressures (Eq. 10)
p_design_ox_pa   = P.p_op_ox_tank_pa   * P.safety_factor; % [Pa] Design pressure for LOX tank
p_design_fuel_pa = P.p_op_fuel_tank_pa * P.safety_factor; % [Pa] Design pressure for Fuel tank

% Calculate cylinder wall thicknesses (Eq. 11, ASME Hoop Stress)
if ~P.prescribe_t   % calculate optimal tank wall thickness only if tank thickness is not being prescribed
    t_cyl_ox_m   = (p_design_ox_pa   * d_ox_tank_m)   / (2 * P.material_allowable_stress_ox_pa   * P.joint_efficiency_ox_tank)   + P.corrosion_allowance_m; % [m] LOX tank cylinder thickness
    t_cyl_fuel_m = (p_design_fuel_pa * d_fuel_tank_m) / (2 * P.material_allowable_stress_fuel_pa * P.joint_efficiency_fuel_tank) + P.corrosion_allowance_m; % [m] Fuel tank cylinder thickness
else
    t_cyl_ox_m   = P.t_cyl_ox_prescribed_m;    % [m] LOX tank cylinder thickness (prescribed)
    t_cyl_fuel_m = P.t_cyl_fuel_prescribed_m;  % [m] Fuel tank cylinder thickness (prescribed)
end

% Calculate end cap volumes
v_caps_ox_m3   = (2/3)*pi*(r_ox_tank_m   - t_cyl_ox_m)^3;   % [m^3] Combined volume of the two 2:1 ellipsoidal end caps for the LOX tank
v_caps_fuel_m3 = (2/3)*pi*(r_fuel_tank_m - t_cyl_fuel_m)^3; % [m^3] Combined volume of the two 2:1 ellipsoidal end caps for the Fuel tank

% Calculate cylindrical section volumes (Eq. 8)
v_cyl_ox_m3   = v_total_ox_tank_m3   - v_caps_ox_m3;   % [m^3] Volume of the LOX tank's cylindrical section
v_cyl_fuel_m3 = v_total_fuel_tank_m3 - v_caps_fuel_m3; % [m^3] Volume of the Fuel tank's cylindrical section

% Calculate cylindrical section lengths (Eq. 9)
l_cyl_ox_m   = v_cyl_ox_m3   / (pi * (r_ox_tank_m   - t_cyl_ox_m)^2);   % [m] Required length of the LOX tank's cylindrical section
l_cyl_fuel_m = v_cyl_fuel_m3 / (pi * (r_fuel_tank_m - t_cyl_fuel_m)^2); % [m] Required length of the Fuel tank's cylindrical section

% Calculate empty tank masses
m_empty_ox_tank_kg   = P.material_density_ox_kgm3   * (pi * (r_ox_tank_m^2   - (r_ox_tank_m   - t_cyl_ox_m)^2)  *l_cyl_ox_m   + (2/3)*pi * (r_ox_tank_m^3   - (r_ox_tank_m   - t_cyl_ox_m)^3));   % [kg] Mass of the empty LOX tank
m_empty_fuel_tank_kg = P.material_density_fuel_kgm3 * (pi * (r_fuel_tank_m^2 - (r_fuel_tank_m - t_cyl_fuel_m)^2)*l_cyl_fuel_m + (2/3)*pi * (r_fuel_tank_m^3 - (r_fuel_tank_m - t_cyl_fuel_m)^3)); % [kg] Mass of the empty Fuel tank

% --- Pressurant Tank ---

% ======================================================================
%  ISOTHERMAL MODEL (OLD/BASELINE) — Conservative worst-case estimate
%  Assumes gas reaches thermal equilibrium with the propellant.
%  This is the original calculation, kept here for comparison.
% ======================================================================
n_ox_mol_isothermal   = (P.p_op_ox_tank_pa   * v_total_ox_tank_m3)   / (C.r_universal_jmolk * (P.ox_temp_k + P.pressurant_temp_k)/2);  % [mol] Moles of gas to displace oxidizer (avg temp of OX and pressurant)
n_fuel_mol_isothermal = (P.p_op_fuel_tank_pa * v_total_fuel_tank_m3) / (C.r_universal_jmolk * P.fuel_temp_k);                          % [mol] Moles of gas to displace fuel

n_total_mol_isothermal           = n_ox_mol_isothermal + n_fuel_mol_isothermal;                    % [mol] Total moles (isothermal)
m_pressurant_gas_isothermal_kg   = n_total_mol_isothermal * P.pressurant_molar_mass_kgmol;         % [kg]  Total pressurant mass (isothermal)
v_storage_isothermal_m3          = (n_total_mol_isothermal * C.r_universal_jmolk * P.pressurant_temp_k) / P.p_storage_pressurant_pa; % [m^3] Storage volume needed (isothermal)
v_storage_isothermal_L           = v_storage_isothermal_m3 * 1000;                                 % [L]   Storage volume in liters (isothermal)

% ======================================================================
%  ISOTHERMAL BOUNDS — Absolute min and max pressurant mass
%  These bracket the entire possible range using the ideal gas law.
%
%  IMPORTANT: The ullage does not fill the entire tank at end of burn
%  because a fraction of propellant remains as residuals (P.residual_fraction).
%  The actual final ullage volume is:
%    V_final = V_tank * [ullage_frac + (1 - ullage_frac)*(1 - residual_frac)]
%  This is the volume the gas must fill, and what the thermal model uses.
%
%  COLD BOUND (Maximum N2): Gas in the LOX tank is at LOX temperature
%    (90 K) and never warms up. This is the absolute worst case — cold
%    gas is dense, so you need the most mass to fill the volume.
%
%  WARM BOUND (Minimum N2): Gas in the LOX tank is at the inlet/storage
%    temperature (294 K) and never cools down. This is the absolute best
%    case — warm gas is less dense, so you need the least mass.
%
%  The fuel tank is at 294 K in both cases (fuel = ambient temp).
% ======================================================================

% Calculate actual final ullage volumes (accounting for residual propellant)
v_final_ullage_ox_m3   = v_total_ox_tank_m3   * (P.ullage_fraction_ox   + (1 - P.ullage_fraction_ox)   * (1 - P.residual_fraction)); % [m^3] LOX tank ullage volume at end of burn
v_final_ullage_fuel_m3 = v_total_fuel_tank_m3 * (P.ullage_fraction_fuel + (1 - P.ullage_fraction_fuel) * (1 - P.residual_fraction)); % [m^3] Fuel tank ullage volume at end of burn

% --- COLD BOUND (gas at LOX temp = 90 K) ---
n_ox_mol_cold   = (P.p_op_ox_tank_pa   * v_final_ullage_ox_m3)   / (C.r_universal_jmolk * P.ox_temp_k);   % [mol] LOX ullage gas at 90 K (worst case)
n_fuel_mol_cold = (P.p_op_fuel_tank_pa * v_final_ullage_fuel_m3) / (C.r_universal_jmolk * P.fuel_temp_k);  % [mol] Fuel ullage gas at 294 K (same as always)

n_total_mol_cold         = n_ox_mol_cold + n_fuel_mol_cold;                                               % [mol] Total moles (cold bound)
m_pressurant_gas_cold_kg = n_total_mol_cold * P.pressurant_molar_mass_kgmol;                              % [kg]  Total pressurant mass (cold bound)
v_storage_cold_m3        = (n_total_mol_cold * C.r_universal_jmolk * P.pressurant_temp_k) / P.p_storage_pressurant_pa; % [m^3] Storage volume (cold bound)
v_storage_cold_L         = v_storage_cold_m3 * 1000;                                                      % [L]   Storage volume in liters (cold bound)

% --- WARM BOUND (gas at inlet temp = 294 K) ---
n_ox_mol_warm   = (P.p_op_ox_tank_pa   * v_final_ullage_ox_m3)   / (C.r_universal_jmolk * P.pressurant_temp_k); % [mol] LOX ullage gas at 294 K (best case)
n_fuel_mol_warm = (P.p_op_fuel_tank_pa * v_final_ullage_fuel_m3) / (C.r_universal_jmolk * P.fuel_temp_k);       % [mol] Fuel ullage gas at 294 K

n_total_mol_warm         = n_ox_mol_warm + n_fuel_mol_warm;                                               % [mol] Total moles (warm bound)
m_pressurant_gas_warm_kg = n_total_mol_warm * P.pressurant_molar_mass_kgmol;                              % [kg]  Total pressurant mass (warm bound)
v_storage_warm_m3        = (n_total_mol_warm * C.r_universal_jmolk * P.pressurant_temp_k) / P.p_storage_pressurant_pa; % [m^3] Storage volume (warm bound)
v_storage_warm_L         = v_storage_warm_m3 * 1000;                                                      % [L]   Storage volume in liters (warm bound)

% ======================================================================
%  TRANSIENT THERMAL MODEL (NEW) — More physically accurate estimate
%  Tracks the actual ullage gas temperature over time using an energy
%  balance. Typically gives 30-60% less mass than the isothermal model
%  for the LOX tank. See calcPressurantThermal.m for full details.
% ======================================================================
if P.use_thermal_model
    [m_pressurant_gas_thermal_kg, v_storage_thermal_m3, v_storage_thermal_L, ...
     ~, ~, T_ox_history_K, T_fuel_history_K, t_history_s] = calcPressurantThermal( ...
        v_total_ox_tank_m3, v_total_fuel_tank_m3,                           ...
        m_dot_ox_kgs, m_dot_fuel_kgs,                                       ...
        t_cyl_ox_m, t_cyl_fuel_m,                                          ...
        d_ox_tank_m, d_fuel_tank_m,                                         ...
        P, C);
else
    % Thermal model disabled — fill outputs with isothermal values
    m_pressurant_gas_thermal_kg = m_pressurant_gas_isothermal_kg;
    v_storage_thermal_m3        = v_storage_isothermal_m3;
    v_storage_thermal_L         = v_storage_isothermal_L;
    T_ox_history_K              = [];
    T_fuel_history_K            = [];
    t_history_s                 = [];
end

% ======================================================================
%  SELECT ACTIVE PRESSURANT MASS
%  When the thermal model is enabled, use it as the design value.
%  When disabled, fall back to the isothermal (conservative) value.
% ======================================================================
if P.use_thermal_model
    m_pressurant_gas_kg = m_pressurant_gas_thermal_kg; % [kg] Use thermal model result
else
    m_pressurant_gas_kg = m_pressurant_gas_isothermal_kg; % [kg] Use isothermal model result (conservative)
end

if P.use_cots_tank
    % --- COTS Tank Path ---
    % Spec sheet values are used directly. COPV material math is skipped.
    % TODO: Replace placeholder values in getInputs.m Section 1.3b with
    %       actual values from the COTS tank spec sheet before using for design.
    m_empty_pressurant_tank_kg = P.cots_tank_mass_kg;           % [kg]  Empty tank mass from spec sheet
    d_pressurant_tank_outer_m  = P.cots_tank_outer_diameter_m;  % [m]   Outer diameter from spec sheet
    r_pressurant_outer_m       = P.cots_tank_outer_diameter_m / 2; % [m] Outer radius from spec sheet
    l_cyl_pressurant_m         = P.cots_tank_length_m;          % [m]   Total tank length used directly in geometry stackup
    % Wall thickness and internal radius are not available until spec sheet is received
    r_pressurant_internal_m    = NaN;   % %%% PLACEHOLDER %%% [m]   Internal radius — unknown until spec sheet received
    t_pressurant_tank_m        = NaN;   % %%% PLACEHOLDER %%% [m]   Wall thickness  — unknown until spec sheet received
    v_pressurant_tank_m3       = P.cots_tank_volume_m3; % [m^3] Internal volume of the COTS tank (from spec sheet)
else
    % --- Calculated COPV Path ---
    % Calculate internal volume of pressurant storage tank (Eq. 17)
    % Use the active pressurant mass (thermal or isothermal, depending on toggle)
    n_active_mol = m_pressurant_gas_kg / P.pressurant_molar_mass_kgmol; % [mol] Active total moles (from whichever model is selected)
    v_pressurant_tank_internal_m3 = (n_active_mol * C.r_universal_jmolk * P.pressurant_temp_k) / P.p_storage_pressurant_pa; % [m^3] Internal volume of the high-pressure storage tank

    % Calculate wall thickness of spherical pressurant tank (Eq. 19)
    p_design_pressurant_pa = P.p_storage_pressurant_pa * P.safety_factor; % [Pa] Design pressure for pressurant tank

    % Define inner and outer and liner radii based on the allowed diameter
    r_pressurant_outer_m = P.d_tank_outer_m / 2; % [m] Outer radius of the pressurant tank (based on allowed diameter)

    % Calculate wall thickness using hoop stress formula for a cylinder
    % We use the tank's outer radius for this calculation
    if ~P.prescribe_t   % calculate optimal tank wall thickness only if tank thickness is not being prescribed
        t_pressurant_tank_m = (p_design_pressurant_pa * r_pressurant_outer_m) / (P.material_allowable_stress_pressurant_pa * P.joint_efficiency_pressurant_tank) + P.corrosion_allowance_m; % [m] Pressurant tank wall thickness
    else
        t_pressurant_tank_m = P.t_pressurant_prescribed_m; % [m] Pressurant tank wall thickness (prescribed)
    end

    r_pressurant_internal_m    = r_pressurant_outer_m - t_pressurant_tank_m - P.t_liner_m; % [m] Final internal radius of the pressurant tank (inside the liner)
    r_pressurant_liner_outer_m = r_pressurant_internal_m + P.t_liner_m;                    % [m] Outer radius of the COPV liner

    % Calculate wall thickness using hoop stress formula for a cylinder
    % We use the liner's outer radius for this calculation
    t_pressurant_tank_m = (p_design_pressurant_pa * r_pressurant_outer_m) / (P.material_allowable_stress_pressurant_pa * P.joint_efficiency_pressurant_tank) + P.corrosion_allowance_m; % [m] Wall thickness

    % Calculate the internal volume of the two 2:1 ellipsoidal end caps
    v_caps_pressurant_m3 = (2/3) * pi * r_pressurant_internal_m^3; % [m^3] Combined internal volume of the two 2:1 ellipsoidal end caps

    % Calculate the required volume and length of the cylindrical section
    v_cyl_pressurant_m3 = v_pressurant_tank_internal_m3 - v_caps_pressurant_m3; % [m^3] Required internal volume of the pressurant tank's cylindrical section
    if v_cyl_pressurant_m3 < 0 % Check if the required volume fits within the end caps alone
        warning('Pressurant tank volume is too small; the end caps alone provide more volume than required. Consider a smaller diameter or spherical tank.');
        l_cyl_pressurant_m  = 0;  % Set length to zero if caps are sufficient
        v_cyl_pressurant_m3 = 0;
    else
        l_cyl_pressurant_m = v_cyl_pressurant_m3 / (pi * r_pressurant_internal_m^2); % [m] Length of cylindrical portion
    end

    % Calculate the mass of the empty pressurant tank (COPV)
    % Mass = (Volume of Carbon Fiber Shell + Volume of Liner) * Respective Densities
    m_shell_cyl_material  = pi * (r_pressurant_outer_m^2      - r_pressurant_liner_outer_m^2) * l_cyl_pressurant_m; % [m^3] Volume of the carbon fiber shell (cylindrical section)
    m_shell_caps_material = (2/3)*pi * (r_pressurant_outer_m^3 - r_pressurant_liner_outer_m^3);                      % [m^3] Volume of the carbon fiber shell (both 2:1 ellipsoidal end caps)
    m_liner_cyl_material  = pi * (r_pressurant_liner_outer_m^2 - r_pressurant_internal_m^2)   * l_cyl_pressurant_m; % [m^3] Volume of the liner (cylindrical section)
    m_liner_caps_material = (2/3)*pi * (r_pressurant_liner_outer_m^3 - r_pressurant_internal_m^3);                   % [m^3] Volume of the liner (both 2:1 ellipsoidal end caps)

    m_empty_pressurant_tank_kg = P.material_density_pressurant_kgm3 * (m_shell_cyl_material + m_shell_caps_material) + P.material_density_liner_kgm3 * (m_liner_cyl_material + m_liner_caps_material); % [kg] Total empty mass of COPV (Shell Mass + Liner Mass)

    % Calculate pressurant tank outer diameter for geometry stackup
    d_pressurant_tank_outer_m = 2 * r_pressurant_outer_m; % [m] Final outer diameter of the pressurant tank (for stackup)
    v_pressurant_tank_m3      = v_pressurant_tank_internal_m3; % [m^3] Calculated internal volume of the pressurant tank

end % end COTS/COPV toggle

%% --- 3. MASS ASSUMPTION (Section 1.6 of getInputs.m) ---
% The combined mass of plumbing and miscellaneous items is assumed equal to
% the total dry tank mass. Computed here because it depends on tank masses.
m_misc_and_plumbing_kg = m_empty_ox_tank_kg + m_empty_fuel_tank_kg + m_empty_pressurant_tank_kg; % [kg] Combined miscellaneous and plumbing mass (assumed equal to total dry tank mass per Section 1.6 of getInputs.m)

%% --- 4. VEHICLE GEOMETRY STACKUP ---
% For 2:1 ellipsoidal caps, each cap's height is r/2. Total added length from two caps is (r/2) + (r/2) = r = d/2.
% NOTE: For the COTS pressurant tank path, l_cyl_pressurant_m is set to the
%       full tank length from the spec sheet, so the + d/2 cap term below
%       will overestimate the total pressurant tank length. TODO: Update
%       once the COTS tank geometry is confirmed from the spec sheet.
l_ox_tank_total_m         = l_cyl_ox_m         + (d_ox_tank_m / 2);               % [m] Length of Lox cylinder + height of two 2:1 ellipsoidal caps
l_fuel_tank_total_m       = l_cyl_fuel_m       + (d_fuel_tank_m / 2);             % [m] Length of Fuel cylinder + height of two 2:1 ellipsoidal caps
l_pressurant_tank_total_m = l_cyl_pressurant_m + (d_pressurant_tank_outer_m / 2); % [m] Length of cylinder + combined height of two 2:1 ellipsoidal caps

% NOTE: This is a simplified length sum. A real stackup would be more complex.
l_total_vehicle_m = l_ox_tank_total_m + l_fuel_tank_total_m + l_pressurant_tank_total_m; % [m] Simplified total vehicle length (tanks only). TODO: Add engine, avionics, recovery, etc.

%% --- 5. VEHICLE MASS BUILDUP & STATIC PERFORMANCE ---

% Calculate total vehicle liftoff mass
m_total_kg = m_empty_ox_tank_kg + m_empty_fuel_tank_kg + m_ox_kg + m_fuel_kg + m_empty_pressurant_tank_kg + m_pressurant_gas_kg + m_misc_and_plumbing_kg; % [kg] Total vehicle liftoff mass

% Calculate final mass after propellant is consumed
m_final_kg = m_total_kg - (1 - P.residual_fraction) * (m_ox_kg + m_fuel_kg); % [kg] Final mass after propellant is consumed

% Calculate static performance metrics
twr_ratio  = (P.f_thrust_n / m_total_kg) / C.g_earth_ms2;              % [unitless] Thrust-to-Weight ratio
delta_v_ms = P.i_sp_s * C.g_earth_ms2 * log(m_total_kg / m_final_kg);  % [m/s] Ideal change in velocity

end