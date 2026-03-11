function [m_gas_total_thermal_kg, v_storage_thermal_m3, v_storage_thermal_L, ...
          m_gas_ox_thermal_kg, m_gas_fuel_thermal_kg,                       ...
          T_ox_history_K, T_fuel_history_K, t_history_s] = calcPressurantThermal( ...
    v_total_ox_tank_m3, v_total_fuel_tank_m3,                               ...
    m_dot_ox_kgs, m_dot_fuel_kgs,                                           ...
    t_cyl_ox_m, t_cyl_fuel_m,                                              ...
    d_ox_tank_m, d_fuel_tank_m,                                             ...
    P, C)
%{
---------------------------------------------------------------------------
calcPressurantThermal
  Calculates the total pressurant gas mass needed to push both propellants
  (LOX and Jet-A) out of their tanks using a transient ullage gas thermal
  model. This is more accurate than the simple isothermal model because it
  tracks the actual gas temperature over time using an energy balance.

  The isothermal model assumes the gas reaches thermal equilibrium with the
  propellant (worst case). This transient model accounts for:
    - Warm gas entering the ullage from the pressurant supply (~294 K)
    - Convective heat loss from the gas to the cold propellant surface
    - Growing ullage volume as propellant drains out

  The result is a lower (more realistic) pressurant mass requirement,
  especially for the LOX tank where the temperature difference is large.

  Reference: pressurant_thermal_analysis.pdf (Sections 4-8)

How It Works (Summary):
  For each propellant tank, the function steps through time from ignition
  to burnout. At each timestep it:
    1. Grows the ullage volume (propellant is leaving)
    2. Calculates how much new gas flows in to maintain pressure
    3. Updates the gas temperature using the energy balance (Eq. 11)
    4. Recalculates the gas mass from the ideal gas law (Eq. 15)
  The final gas mass at the end of the burn is the answer for that tank.
  Both tanks are summed to get the total pressurant requirement.

Inputs:
  v_total_ox_tank_m3   - [m^3]    Total internal volume of the LOX tank
  v_total_fuel_tank_m3 - [m^3]    Total internal volume of the Fuel tank
  m_dot_ox_kgs         - [kg/s]   Mass flow rate of LOX leaving the tank
  m_dot_fuel_kgs       - [kg/s]   Mass flow rate of Fuel leaving the tank
  t_cyl_ox_m           - [m]      Wall thickness of the LOX tank
  t_cyl_fuel_m         - [m]      Wall thickness of the Fuel tank
  d_ox_tank_m          - [m]      Outer diameter of the LOX tank
  d_fuel_tank_m        - [m]      Outer diameter of the Fuel tank
  P  - [struct] Design inputs from getInputs()
  C  - [struct] Physical constants from getConstants()

Outputs:
  m_gas_total_thermal_kg  - [kg]    Total pressurant gas mass (LOX + Fuel)
  v_storage_thermal_m3    - [m^3]   Required pressurant storage volume
  v_storage_thermal_L     - [L]     Required pressurant storage volume in liters
  m_gas_ox_thermal_kg     - [kg]    Pressurant gas mass for LOX tank only
  m_gas_fuel_thermal_kg   - [kg]    Pressurant gas mass for Fuel tank only
  T_ox_history_K          - [K]     LOX ullage gas temperature over time (for plotting)
  T_fuel_history_K        - [K]     Fuel ullage gas temperature over time (for plotting)
  t_history_s             - [s]     Time array corresponding to the temperature histories
---------------------------------------------------------------------------
%}

%% --- NITROGEN GAS PROPERTIES ---
% These are derived from the universal gas constant and N2 molar mass.
% See pressurant_thermal_analysis.pdf Section 6.

R_s_jkgk = C.r_universal_jmolk / P.pressurant_molar_mass_kgmol; % [J/(kg*K)] Specific gas constant for N2 (Eq. 17, ~296.8)
cp_n2_jkgk = P.cp_n2_jkgk;                                      % [J/(kg*K)] Specific heat at constant pressure for N2

%% --- TIME SETUP ---
dt_s = P.dt_thermal_s;                              % [s]     Timestep for the thermal simulation
N_steps = round(P.t_burn_s / dt_s);                 % [steps] Number of timesteps
t_history_s = (0:N_steps) * dt_s;                   % [s]     Time array (from 0 to burn time)

%% ====================================================================
%  LOX TANK - Transient Thermal Model
%  This is where the thermal model matters most. The gas enters warm
%  (~294 K) but the LOX is at ~90 K, so there is a large temperature
%  difference driving convective heat loss.
%  ====================================================================

% --- LOX Tank Parameters ---
p_ox_pa         = P.p_op_ox_tank_pa;        % [Pa]     LOX tank operating pressure
T_prop_ox_K     = P.ox_temp_k;              % [K]      LOX surface temperature
T_inlet_K       = P.pressurant_temp_k;      % [K]      Temperature of incoming N2 gas
h_ox_wm2k       = P.h_convection_ox_wm2k;  % [W/(m^2*K)] Heat transfer coefficient for LOX tank
rho_ox_kgm3     = P.ox_density_kgm3;       % [kg/m^3] Density of liquid oxygen

% Gas-liquid contact area (Eq. 16): cross-sectional area of the tank interior
d_ox_inner_m    = d_ox_tank_m - 2 * t_cyl_ox_m;          % [m]   Inner diameter of LOX tank
A_contact_ox_m2 = pi * (d_ox_inner_m / 2)^2;             % [m^2] Gas-liquid interface area

% Initial conditions (Section 4.7 of the PDF)
V_ullage_ox_m3  = P.ullage_fraction_ox * v_total_ox_tank_m3;            % [m^3] Initial ullage volume (Eq. 8)
T_gas_ox_K      = T_prop_ox_K;                                          % [K]   Initial gas temp = propellant temp (thermal equilibrium pre-launch)
m_gas_ox_kg     = p_ox_pa * V_ullage_ox_m3 / (R_s_jkgk * T_gas_ox_K); % [kg]  Initial gas mass from ideal gas law (Eq. 15)

% Volumetric flow rate of LOX leaving the tank (Eq. 6)
dVdt_ox_m3s = m_dot_ox_kgs / rho_ox_kgm3; % [m^3/s] Rate at which ullage volume grows

% Pre-allocate history arrays for plotting
T_ox_history_K = zeros(1, N_steps + 1); % [K] Temperature history
T_ox_history_K(1) = T_gas_ox_K;         % Store initial condition

% --- LOX Tank Time-Stepping Loop ---
for i = 1:N_steps

    % Step 5: Grow the ullage volume as LOX drains out (Eq. 7)
    dV_ox_m3 = dVdt_ox_m3s * dt_s;                                     % [m^3] Volume increment this timestep
    V_ullage_ox_m3 = V_ullage_ox_m3 + dV_ox_m3;                        % [m^3] Updated ullage volume
    V_ullage_ox_m3 = min(V_ullage_ox_m3, v_total_ox_tank_m3);          % [m^3] Clamp to tank volume (can't exceed tank)

    % Step 7: Calculate incoming gas mass to maintain pressure (Eq. 9)
    dm_in_ox_kg = p_ox_pa * dV_ox_m3 / (R_s_jkgk * T_gas_ox_K);      % [kg]    Mass of new gas entering this timestep
    mdot_in_ox_kgs = dm_in_ox_kg / dt_s;                               % [kg/s]  Instantaneous mass flow rate of pressurant

    % Step 8: Calculate temperature changes from the energy balance (Eq. 12, 13)
    % Term 1: Warm gas mixing — incoming warm gas heats up the ullage
    dT_mixing_K    = (mdot_in_ox_kgs / m_gas_ox_kg) * (T_inlet_K - T_gas_ox_K) * dt_s;           % [K] (Eq. 12)

    % Term 2: Convective heat loss — cold LOX surface cools the gas
    dT_heatloss_K  = -(h_ox_wm2k * A_contact_ox_m2 / (m_gas_ox_kg * cp_n2_jkgk)) * (T_gas_ox_K - T_prop_ox_K) * dt_s; % [K] (Eq. 13)

    % Step 9: Update the gas temperature (Eq. 14)
    T_gas_ox_K = T_gas_ox_K + dT_mixing_K + dT_heatloss_K;            % [K] New gas temperature

    % Step 10: Recalculate gas mass at the updated temperature and volume (Eq. 15)
    m_gas_ox_kg = p_ox_pa * V_ullage_ox_m3 / (R_s_jkgk * T_gas_ox_K); % [kg] Updated gas mass

    % Step 11: Store temperature for plotting
    T_ox_history_K(i + 1) = T_gas_ox_K;

end

% The final gas mass in the LOX ullage at the end of the burn
m_gas_ox_thermal_kg = m_gas_ox_kg; % [kg] Final LOX pressurant mass (thermal model)

%% ====================================================================
%  FUEL TANK - Transient Thermal Model
%  For Jet-A, the gas enters at ~294 K and the fuel is also at ~294 K,
%  so the temperature difference is essentially zero. The thermal model
%  gives the same result as the isothermal model for the fuel
%  tank, but we run it anyway for completeness and consistency.
%  ====================================================================

% --- Fuel Tank Parameters ---
p_fuel_pa        = P.p_op_fuel_tank_pa;        % [Pa]     Fuel tank operating pressure
T_prop_fuel_K    = P.fuel_temp_k;              % [K]      Fuel surface temperature
h_fuel_wm2k      = P.h_convection_fuel_wm2k;  % [W/(m^2*K)] Heat transfer coefficient for Fuel tank
rho_fuel_kgm3    = P.fuel_density_kgm3;        % [kg/m^3] Density of Jet-A

% Gas-liquid contact area (Eq. 16)
d_fuel_inner_m    = d_fuel_tank_m - 2 * t_cyl_fuel_m;        % [m]   Inner diameter of Fuel tank
A_contact_fuel_m2 = pi * (d_fuel_inner_m / 2)^2;             % [m^2] Gas-liquid interface area

% Initial conditions
V_ullage_fuel_m3  = P.ullage_fraction_fuel * v_total_fuel_tank_m3;                % [m^3] Initial ullage volume
T_gas_fuel_K      = T_prop_fuel_K;                                                % [K]   Initial gas temp = fuel temp
m_gas_fuel_kg     = p_fuel_pa * V_ullage_fuel_m3 / (R_s_jkgk * T_gas_fuel_K);   % [kg]  Initial gas mass

% Volumetric flow rate of fuel leaving the tank
dVdt_fuel_m3s = m_dot_fuel_kgs / rho_fuel_kgm3; % [m^3/s] Rate of ullage volume growth

% Pre-allocate
T_fuel_history_K = zeros(1, N_steps + 1);
T_fuel_history_K(1) = T_gas_fuel_K;

% --- Fuel Tank Time-Stepping Loop ---
for i = 1:N_steps

    % Grow ullage volume
    dV_fuel_m3 = dVdt_fuel_m3s * dt_s;
    V_ullage_fuel_m3 = V_ullage_fuel_m3 + dV_fuel_m3;
    V_ullage_fuel_m3 = min(V_ullage_fuel_m3, v_total_fuel_tank_m3);

    % Incoming gas mass
    dm_in_fuel_kg = p_fuel_pa * dV_fuel_m3 / (R_s_jkgk * T_gas_fuel_K);
    mdot_in_fuel_kgs = dm_in_fuel_kg / dt_s;

    % Temperature changes
    dT_mixing_K   = (mdot_in_fuel_kgs / m_gas_fuel_kg) * (T_inlet_K - T_gas_fuel_K) * dt_s;
    dT_heatloss_K = -(h_fuel_wm2k * A_contact_fuel_m2 / (m_gas_fuel_kg * cp_n2_jkgk)) * (T_gas_fuel_K - T_prop_fuel_K) * dt_s;

    % Update temperature
    T_gas_fuel_K = T_gas_fuel_K + dT_mixing_K + dT_heatloss_K;

    % Update gas mass
    m_gas_fuel_kg = p_fuel_pa * V_ullage_fuel_m3 / (R_s_jkgk * T_gas_fuel_K);

    % Store
    T_fuel_history_K(i + 1) = T_gas_fuel_K;

end

% Final fuel tank pressurant mass
m_gas_fuel_thermal_kg = m_gas_fuel_kg; % [kg] Final Fuel pressurant mass (thermal model)

%% ====================================================================
%  TOTAL PRESSURANT REQUIREMENT
%  Sum both tanks and calculate the required storage volume.
%  ====================================================================

% Total pressurant gas mass (Eq. 18)
m_gas_total_thermal_kg = m_gas_ox_thermal_kg + m_gas_fuel_thermal_kg; % [kg] Total pressurant gas mass

% Convert mass to moles for storage volume calculation
n_total_thermal_mol = m_gas_total_thermal_kg / P.pressurant_molar_mass_kgmol; % [mol] Total moles of pressurant

% Required pressurant storage volume at supply conditions (Eq. 5)
v_storage_thermal_m3 = (n_total_thermal_mol * C.r_universal_jmolk * P.pressurant_temp_k) / P.p_storage_pressurant_pa; % [m^3] Storage volume needed
v_storage_thermal_L  = v_storage_thermal_m3 * 1000; % [L] Storage volume in liters

end
