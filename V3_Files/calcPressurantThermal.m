function [m_gas_total_thermal_kg, v_storage_thermal_m3, v_storage_thermal_L, ...
          m_gas_ox_thermal_kg, m_gas_fuel_thermal_kg,                       ...
          T_ox_history_K, T_fuel_history_K, t_history_s,                    ...
          T_pt_history_K, T_wall_history_K,                                 ...
          p_pt_history_pa, m_pt_history_kg,                                 ...
          m_ox_ullage_history_kg, m_fuel_ullage_history_kg,                 ...
          m_gas_total_v1_kg, v_storage_v1_L] = calcPressurantThermal(       ...
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
  model.

  This function contains TWO models:
    V1 (Old/Default):  The original transient ullage thermal model.
       Assumes the pressurant gas enters at a CONSTANT temperature
       (P.pressurant_fill_temp_k = 294 K) for the entire burn.
       Two independent loops (one per tank).  Always runs for comparison.

    V2 (New/Blowdown):  Coupled pressurant supply blowdown + ullage model.
       Tracks the pressurant tank gas temperature as it cools due to
       expansion during the burn.
       The inlet temperature to the propellant tanks decreases over time,
       increasing the total pressurant mass requirement.  Also tracks the
       pressurant tank wall temperature (lumped capacitance model).
       Both propellant tanks are iterated in a SINGLE shared loop so the
       blowdown model can sum their combined mass demand each timestep.
       Only runs when P.use_blowdown_model = true.

  ENERGY BALANCE METHOD (used by both V1 and V2):
    Each propellant tank ullage is at constant pressure (regulated).
    At each timestep, the enthalpy balance is solved exactly:

      m_new * cp * T_new = m_old * cp * T_old + dm_in * cp * T_inlet + Q

    Combined with the ideal gas constraint at constant pressure:

      m_new * T_new = p * V_new / R_s

    This gives a closed-form solution for dm_in (mass entering from the
    pressurant tank) that guarantees mass conservation at every timestep:

      dm_in = (p*V_new/R_s - m_old*T_old - Q/cp) / T_inlet

    where Q = h * A * (T_prop - T_old) * dt is the heat transfer from the
    propellant surface to the gas (positive when propellant warms the gas,
    negative when gas loses heat to cold propellant like LOX).

    If dm_in < 0 (gas temperature rose enough that existing mass would
    overpressure at current volume), we set dm_in = 0.  Physically, the
    regulator does not allow backflow.  In this case, the temperature is
    updated from heat transfer only, and the ullage pressure may differ
    slightly from the regulated set point.  This only occurs during the
    first few timesteps when cold pre-pressed gas is rapidly heated.

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
  m_gas_total_thermal_kg  - [kg]    Total pressurant gas mass loaded into PT (active model: V1 or V2)
  v_storage_thermal_m3    - [m^3]   Required pressurant storage volume (active model)
  v_storage_thermal_L     - [L]     Required pressurant storage volume in liters (active model)
  m_gas_ox_thermal_kg     - [kg]    Pressurant gas mass in LOX ullage at end of burn (active model)
  m_gas_fuel_thermal_kg   - [kg]    Pressurant gas mass in Fuel ullage at end of burn (active model)
  T_ox_history_K          - [K]     LOX ullage gas temperature over time (active model)
  T_fuel_history_K        - [K]     Fuel ullage gas temperature over time (active model)
  t_history_s             - [s]     Time array for all history outputs
  T_pt_history_K          - [K]     Pressurant tank gas temperature over time (V2 only, empty if V1)
  T_wall_history_K        - [K]     Pressurant tank wall temperature over time (V2 only, empty if V1)
  p_pt_history_pa         - [Pa]    Pressurant tank pressure over time (V2 only, empty if V1)
  m_pt_history_kg         - [kg]    Pressurant tank remaining gas mass over time (V2 only, empty if V1)
  m_ox_ullage_history_kg  - [kg]    LOX ullage gas mass over time (active model)
  m_fuel_ullage_history_kg- [kg]    Fuel ullage gas mass over time (active model)
  m_gas_total_v1_kg       - [kg]    V1 (old) thermal model total pressurant mass (for comparison)
  v_storage_v1_L          - [L]     V1 (old) thermal model storage volume (for comparison)
---------------------------------------------------------------------------
%}

%% ====================================================================
%  COMMON SETUP - Used by both V1 and V2
%  ====================================================================

% --- Nitrogen Gas Properties ---
R_s_jkgk    = C.r_universal_jmolk / P.pressurant_molar_mass_kgmol; % [J/(kg*K)] Specific gas constant for N2 (~296.8)
cp_n2_jkgk  = P.cp_n2_jkgk;                                        % [J/(kg*K)] Specific heat at constant pressure for N2
cv_n2_jkgk  = P.cv_n2_jkgk;                                        % [J/(kg*K)] Specific heat at constant volume for N2

% --- Time Setup ---
dt_s    = P.dt_thermal_s;                % [s]     Timestep
N_steps = round(P.t_burn_s / dt_s);     % [steps] Number of timesteps
t_history_s = (0:N_steps) * dt_s;        % [s]     Time array (0 to burn time)

% --- LOX Tank Geometry ---
p_ox_pa          = P.p_op_ox_tank_pa;        % [Pa]     LOX tank operating pressure
T_prop_ox_K      = P.ox_temp_k;              % [K]      LOX surface temperature (90 K)
h_ox_wm2k        = P.h_convection_ox_wm2k;  % [W/(m^2*K)] HTC for LOX tank ullage
rho_ox_kgm3      = P.ox_density_kgm3;       % [kg/m^3] Density of liquid oxygen
d_ox_inner_m     = d_ox_tank_m - 2 * t_cyl_ox_m;        % [m]   Inner diameter of LOX tank
A_contact_ox_m2  = pi * (d_ox_inner_m / 2)^2;           % [m^2] Gas-liquid interface area (cross-section of tank)
dVdt_ox_m3s      = m_dot_ox_kgs / rho_ox_kgm3;          % [m^3/s] Rate of LOX ullage volume growth (propellant leaves)

% --- Fuel Tank Geometry ---
p_fuel_pa         = P.p_op_fuel_tank_pa;        % [Pa]     Fuel tank operating pressure
T_prop_fuel_K     = P.fuel_temp_k;              % [K]      Fuel surface temperature (294 K, ambient)
h_fuel_wm2k       = P.h_convection_fuel_wm2k;  % [W/(m^2*K)] HTC for Fuel tank ullage
rho_fuel_kgm3     = P.fuel_density_kgm3;        % [kg/m^3] Density of Jet-A
d_fuel_inner_m    = d_fuel_tank_m - 2 * t_cyl_fuel_m;    % [m]   Inner diameter of Fuel tank
A_contact_fuel_m2 = pi * (d_fuel_inner_m / 2)^2;         % [m^2] Gas-liquid interface area
dVdt_fuel_m3s     = m_dot_fuel_kgs / rho_fuel_kgm3;      % [m^3/s] Rate of Fuel ullage volume growth

% --- Fill / Inlet Temperature ---
T_fill_K = P.pressurant_fill_temp_k; % [K] Temperature at which the pressurant tank was filled (294 K)

% --- Maximum ullage volumes (accounting for residual propellant) ---
% At end of burn, some propellant remains as residuals.  The ullage cannot
% fill the entire tank.  Maximum ullage volume is:
%   V_ullage_max = V_tank * [ullage_frac + (1 - ullage_frac) * (1 - residual_frac)]
V_ullage_ox_max_m3   = v_total_ox_tank_m3   * (P.ullage_fraction_ox   + (1 - P.ullage_fraction_ox)   * (1 - P.residual_fraction)); % [m^3]
V_ullage_fuel_max_m3 = v_total_fuel_tank_m3 * (P.ullage_fraction_fuel + (1 - P.ullage_fraction_fuel) * (1 - P.residual_fraction)); % [m^3]

%% ====================================================================
%  V1: ORIGINAL THERMAL MODEL (constant inlet temperature)
%  This is the old model from the February 2026 analysis.
%  It always runs so results are available for comparison.
%  Two independent loops — one for LOX, one for Fuel.
%  Inlet temperature is constant at T_fill (294 K).
%
%  Energy balance method:
%    dm_in = (p*V_new/R - m_old*T_old - Q/cp) / T_inlet
%    See function header for derivation.
%  ====================================================================

% ---------- V1 LOX Tank ----------
% Initial conditions: ullage gas is pre-pressed and sits at LOX temperature
V_ullage_ox_v1_m3  = P.ullage_fraction_ox * v_total_ox_tank_m3;                    % [m^3] Initial LOX ullage volume
T_gas_ox_v1_K      = T_prop_ox_K;                                                  % [K]   Initial gas temp = LOX temp (90 K)
m_gas_ox_v1_kg     = p_ox_pa * V_ullage_ox_v1_m3 / (R_s_jkgk * T_gas_ox_v1_K);   % [kg]  Initial gas mass from ideal gas law
T_ox_v1_hist_K     = zeros(1, N_steps + 1);                                        % [K]   Pre-allocate temperature history
T_ox_v1_hist_K(1)  = T_gas_ox_v1_K;                                                %       Store initial condition

for i = 1:N_steps
    % Volume increment from propellant leaving the tank
    dV_ox_m3           = dVdt_ox_m3s * dt_s;                                        % [m^3] Volume freed this timestep
    V_ullage_ox_v1_m3  = min(V_ullage_ox_v1_m3 + dV_ox_m3, V_ullage_ox_max_m3);   % [m^3] Updated ullage volume (capped at max)

    % Heat transfer from propellant surface to gas
    % Q > 0 when T_prop > T_gas (propellant warms gas)
    % Q < 0 when T_prop < T_gas (gas loses heat to cold propellant, e.g. LOX)
    Q_ox_J = h_ox_wm2k * A_contact_ox_m2 * (T_prop_ox_K - T_gas_ox_v1_K) * dt_s;  % [J] Heat transfer this timestep

    % Closed-form mass inflow from enthalpy balance + ideal gas constraint:
    %   dm_in = (p*V_new/R - m_old*T_old - Q/cp) / T_inlet
    dm_in_ox_kg = (p_ox_pa * V_ullage_ox_v1_m3 / R_s_jkgk ...
                   - m_gas_ox_v1_kg * T_gas_ox_v1_K ...
                   - Q_ox_J / cp_n2_jkgk) / T_fill_K;                              % [kg] Mass entering from pressurant system

    if dm_in_ox_kg >= 0
        % Normal case: gas flows in from pressurant tank
        m_gas_ox_v1_kg  = m_gas_ox_v1_kg + dm_in_ox_kg;                            % [kg] Updated ullage gas mass
        T_gas_ox_v1_K   = p_ox_pa * V_ullage_ox_v1_m3 / (R_s_jkgk * m_gas_ox_v1_kg); % [K] Temperature from ideal gas law
    else
        % No inflow: regulator does not allow backflow.
        % Gas mass stays the same; only heat transfer changes temperature.
        dm_in_ox_kg     = 0;                                                        % [kg] No mass enters
        T_gas_ox_v1_K   = T_gas_ox_v1_K + Q_ox_J / (m_gas_ox_v1_kg * cp_n2_jkgk); % [K] Temperature from heat transfer only
    end

    T_ox_v1_hist_K(i+1) = T_gas_ox_v1_K;                                           % Store temperature for plotting
end
m_gas_ox_v1_final_kg = m_gas_ox_v1_kg; % [kg] V1 LOX tank final pressurant mass in ullage

% ---------- V1 Fuel Tank ----------
% Initial conditions: ullage gas is pre-pressed and sits at fuel temperature
V_ullage_fuel_v1_m3  = P.ullage_fraction_fuel * v_total_fuel_tank_m3;                     % [m^3] Initial Fuel ullage volume
T_gas_fuel_v1_K      = T_prop_fuel_K;                                                     % [K]   Initial gas temp = Fuel temp (294 K)
m_gas_fuel_v1_kg     = p_fuel_pa * V_ullage_fuel_v1_m3 / (R_s_jkgk * T_gas_fuel_v1_K);  % [kg]  Initial gas mass from ideal gas law
T_fuel_v1_hist_K     = zeros(1, N_steps + 1);                                             % [K]   Pre-allocate
T_fuel_v1_hist_K(1)  = T_gas_fuel_v1_K;

for i = 1:N_steps
    dV_fuel_m3           = dVdt_fuel_m3s * dt_s;                                           % [m^3] Volume freed this timestep
    V_ullage_fuel_v1_m3  = min(V_ullage_fuel_v1_m3 + dV_fuel_m3, V_ullage_fuel_max_m3);   % [m^3] Updated ullage volume

    % Heat transfer (for fuel tank at ambient temp, T_prop ~ T_gas, so Q ~ 0)
    Q_fuel_J = h_fuel_wm2k * A_contact_fuel_m2 * (T_prop_fuel_K - T_gas_fuel_v1_K) * dt_s; % [J]

    % Closed-form mass inflow
    dm_in_fuel_kg = (p_fuel_pa * V_ullage_fuel_v1_m3 / R_s_jkgk ...
                     - m_gas_fuel_v1_kg * T_gas_fuel_v1_K ...
                     - Q_fuel_J / cp_n2_jkgk) / T_fill_K;                                  % [kg]

    if dm_in_fuel_kg >= 0
        m_gas_fuel_v1_kg  = m_gas_fuel_v1_kg + dm_in_fuel_kg;                              % [kg] Updated mass
        T_gas_fuel_v1_K   = p_fuel_pa * V_ullage_fuel_v1_m3 / (R_s_jkgk * m_gas_fuel_v1_kg); % [K] Temp from ideal gas
    else
        dm_in_fuel_kg     = 0;                                                              % [kg] No mass enters
        T_gas_fuel_v1_K   = T_gas_fuel_v1_K + Q_fuel_J / (m_gas_fuel_v1_kg * cp_n2_jkgk); % [K] Heat transfer only
    end

    T_fuel_v1_hist_K(i+1) = T_gas_fuel_v1_K;
end
m_gas_fuel_v1_final_kg = m_gas_fuel_v1_kg; % [kg] V1 Fuel tank final pressurant mass in ullage

% --- V1 Totals ---
% Total pressurant mass = gas in LOX ullage + gas in Fuel ullage at end of burn.
% This is all the N2 that must be available: the pre-pressed gas plus the
% gas delivered during the burn.
m_gas_total_v1_kg    = m_gas_ox_v1_final_kg + m_gas_fuel_v1_final_kg;                                  % [kg]  V1 total pressurant mass
n_total_v1_mol       = m_gas_total_v1_kg / P.pressurant_molar_mass_kgmol;                              % [mol] V1 total moles
v_storage_v1_m3      = (n_total_v1_mol * C.r_universal_jmolk * T_fill_K) / P.p_storage_pressurant_pa;  % [m^3] V1 storage volume (at fill conditions)
v_storage_v1_L       = v_storage_v1_m3 * 1000;                                                          % [L]   V1 storage volume in liters

%% ====================================================================
%  V2: COUPLED BLOWDOWN + THERMAL MODEL (if enabled)
%  ====================================================================
%
%  The pressurant tank is a rigid, fixed-volume vessel.  As gas exits
%  through the regulator, the remaining gas expands and cools 
%  The aluminum tank wall provides some heat back to the gas.
%
%  The inlet temperature to the propellant tanks equals the current
%  pressurant tank gas temperature, which decreases over time.
%  This couples the blowdown model to the ullage thermal model.
%
%  For COTS tanks: the tank volume is fixed from the spec sheet.
%  For custom COPV: bisection is used to find the minimum tank volume
%  such that the final PT pressure stays above the propellant tank
%  operating pressure (otherwise the pressurant can't push propellant out).
%  ====================================================================

if P.use_blowdown_model

    % --- Blowdown model parameters ---
    h_wall_wm2k  = P.h_wall_pressurant_wm2k;   % [W/(m^2*K)] Internal HTC for pressurant tank wall
    c_wall_jkgk  = P.c_wall_pressurant_jkgk;   % [J/(kg*K)]  Specific heat of tank wall (Al 6061)

    if P.use_cots_tank
        % COTS tank: fixed geometry from spec sheet
        V_pt_m3     = P.cots_tank_volume_m3;        % [m^3] Internal volume
        m_wall_kg   = P.cots_tank_wall_mass_kg;      % [kg]  Tank wall mass (for thermal model)
        A_wall_m2   = P.cots_tank_internal_area_m2;  % [m^2] Internal surface area
        max_iters   = 1;                              % No bisection needed — fixed tank
    else
        % Custom COPV: bisect on volume so final PT pressure >= propellant tank pressure.
        % The PT must stay above operating pressure to push propellant out.
        p_final_target = max(p_ox_pa, p_fuel_pa);    % [Pa] Minimum acceptable final PT pressure
        V_low  = v_storage_v1_m3;                     % [m^3] Lower bound (V1 assumes full depletion)
        V_high = v_storage_v1_m3 * 20;                % [m^3] Upper bound (generous)
        max_iters = 30;                                % 30 bisection iterations give ~1e-9 relative accuracy
    end

    for iter = 1:max_iters

        % --- Set pressurant tank volume for this iteration ---
        if ~P.use_cots_tank
            if iter == 1
                % Initial guess: isothermal blowdown from storage to final target pressure
                V_pt_m3 = m_gas_total_v1_kg * R_s_jkgk * T_fill_K ...
                          / (P.p_storage_pressurant_pa - p_final_target);           % [m^3]
            else
                V_pt_m3 = (V_low + V_high) / 2;                                    % [m^3] Bisection midpoint
            end

            % Estimate wall properties assuming a spherical tank at this volume
            r_sph_m   = ((3 * V_pt_m3) / (4 * pi))^(1/3);                          % [m]   Equivalent sphere radius
            A_wall_m2 = 4 * pi * r_sph_m^2;                                         % [m^2] Sphere surface area

            % Wall thickness from hoop stress: t = p*r / (2*sigma)
            t_wall_m  = (P.p_storage_pressurant_pa * r_sph_m) ...
                        / (2 * P.material_allowable_stress_pressurant_pa);           % [m]
            t_wall_m  = max(t_wall_m, 0.003);                                       % [m]   Floor at 3 mm

            % Wall mass = surface area * thickness * density
            m_wall_kg = A_wall_m2 * t_wall_m * P.material_density_pressurant_kgm3;  % [kg]
            m_wall_kg = max(m_wall_kg, 1.0);                                        % [kg]  Floor at 1 kg
        end

        % ================================================================
        %  INITIAL CONDITIONS
        % ================================================================

        % --- Pressurant Tank ---
        T_pt_K   = T_fill_K;                                                        % [K]  Gas starts at fill temperature (at t=0)
        T_wall_K = T_fill_K;                                                        % [K]  Wall starts at fill temperature
        m_pt_kg  = P.p_storage_pressurant_pa * V_pt_m3 / (R_s_jkgk * T_fill_K);   % [kg] Initial gas mass in PT (ideal gas)

        % --- LOX Ullage ---
        % Pre-pressed gas is at LOX temperature (loaded before launch, sat at equilibrium)
        V_ullage_ox_m3   = P.ullage_fraction_ox * v_total_ox_tank_m3;               % [m^3] Initial LOX ullage volume
        T_gas_ox_K       = T_prop_ox_K;                                              % [K]   Initial gas temp = LOX temp (90 K)
        m_gas_ox_kg      = p_ox_pa * V_ullage_ox_m3 / (R_s_jkgk * T_gas_ox_K);    % [kg]  Initial gas mass

        % --- Fuel Ullage ---
        % Pre-pressed gas is at fuel temperature (ambient, ~294 K)
        V_ullage_fuel_m3 = P.ullage_fraction_fuel * v_total_fuel_tank_m3;            % [m^3] Initial Fuel ullage volume
        T_gas_fuel_K     = T_prop_fuel_K;                                            % [K]   Initial gas temp = Fuel temp (294 K)
        m_gas_fuel_kg    = p_fuel_pa * V_ullage_fuel_m3 / (R_s_jkgk * T_gas_fuel_K); % [kg] Initial gas mass

        % --- Total initial N2 in system (for conservation check) ---
        m_total_n2_initial_kg = m_pt_kg + m_gas_ox_kg + m_gas_fuel_kg;              % [kg]

        % --- Allocate history arrays only on the final iteration ---
        is_final_iter = (iter == max_iters);
        if is_final_iter
            T_ox_history_K           = zeros(1, N_steps + 1); T_ox_history_K(1)           = T_gas_ox_K;
            T_fuel_history_K         = zeros(1, N_steps + 1); T_fuel_history_K(1)         = T_gas_fuel_K;
            T_pt_history_K           = zeros(1, N_steps + 1); T_pt_history_K(1)           = T_pt_K;
            T_wall_history_K         = zeros(1, N_steps + 1); T_wall_history_K(1)         = T_wall_K;
            m_pt_history_kg          = zeros(1, N_steps + 1); m_pt_history_kg(1)          = m_pt_kg;
            p_pt_history_pa          = zeros(1, N_steps + 1); p_pt_history_pa(1)          = P.p_storage_pressurant_pa;
            m_ox_ullage_history_kg   = zeros(1, N_steps + 1); m_ox_ullage_history_kg(1)   = m_gas_ox_kg;
            m_fuel_ullage_history_kg = zeros(1, N_steps + 1); m_fuel_ullage_history_kg(1) = m_gas_fuel_kg;
        end

        % ================================================================
        %  COUPLED TIME-STEPPING LOOP
        %
        %  At each timestep:
        %    1. Set inlet temperature = current PT gas temperature
        %    2. Solve LOX ullage energy balance for dm_in_ox
        %    3. Solve Fuel ullage energy balance for dm_in_fuel
        %    4. Update PT state: mass, temperature,
        %       wall temperature
        %    5. Store histories
        % ================================================================
        for i = 1:N_steps

            % --- Step 1: Inlet temperature = PT gas temperature ---
            T_inlet_K = T_pt_K;                                                     % [K]

            % ==========================================================
            %  Step 2: LOX TANK ULLAGE UPDATE
            % ==========================================================

            % Volume increment from LOX leaving the tank
            dV_ox_m3       = dVdt_ox_m3s * dt_s;                                    % [m^3]
            V_ullage_ox_m3 = min(V_ullage_ox_m3 + dV_ox_m3, V_ullage_ox_max_m3);   % [m^3]

            % Heat transfer from LOX surface to ullage gas
            %   Q > 0 when propellant is warmer than gas (gas gains heat)
            %   Q < 0 when gas is warmer than propellant (gas loses heat to LOX)
            %   For LOX (90 K): initially Q ≈ 0 (gas also at 90 K), then Q < 0
            %   as hot pressurant enters and raises gas temp above 90 K.
            Q_ox_J = h_ox_wm2k * A_contact_ox_m2 * (T_prop_ox_K - T_gas_ox_K) * dt_s; % [J]

            % Closed-form mass inflow (see function header for derivation):
            %   dm_in = (p*V_new/R - m_old*T_old - Q/cp) / T_inlet
            dm_in_ox_kg = (p_ox_pa * V_ullage_ox_m3 / R_s_jkgk ...
                           - m_gas_ox_kg * T_gas_ox_K ...
                           - Q_ox_J / cp_n2_jkgk) / T_inlet_K;                     % [kg]

            if dm_in_ox_kg >= 0
                % Normal: gas flows from PT into LOX ullage
                m_gas_ox_kg = m_gas_ox_kg + dm_in_ox_kg;                            % [kg]
                T_gas_ox_K  = p_ox_pa * V_ullage_ox_m3 / (R_s_jkgk * m_gas_ox_kg);% [K]
            else
                % No inflow: regulator blocks backflow. Only heat transfer.
                dm_in_ox_kg = 0;                                                    % [kg]
                T_gas_ox_K  = T_gas_ox_K + Q_ox_J / (m_gas_ox_kg * cp_n2_jkgk);  % [K]
            end

            % ==========================================================
            %  Step 3: FUEL TANK ULLAGE UPDATE
            % ==========================================================

            dV_fuel_m3       = dVdt_fuel_m3s * dt_s;                                         % [m^3]
            V_ullage_fuel_m3 = min(V_ullage_fuel_m3 + dV_fuel_m3, V_ullage_fuel_max_m3);    % [m^3]

            % Heat transfer from fuel surface to ullage gas
            %   For fuel at 294 K: Q ≈ 0 when inlet is also ~294 K (early in burn).
            %   As PT cools (blowdown), inlet gas gets colder, so ullage cools and
            %   Q > 0 (fuel surface heats the gas). Handled automatically by sign.
            Q_fuel_J = h_fuel_wm2k * A_contact_fuel_m2 * (T_prop_fuel_K - T_gas_fuel_K) * dt_s; % [J]

            dm_in_fuel_kg = (p_fuel_pa * V_ullage_fuel_m3 / R_s_jkgk ...
                             - m_gas_fuel_kg * T_gas_fuel_K ...
                             - Q_fuel_J / cp_n2_jkgk) / T_inlet_K;                  % [kg]

            if dm_in_fuel_kg >= 0
                m_gas_fuel_kg = m_gas_fuel_kg + dm_in_fuel_kg;                       % [kg]
                T_gas_fuel_K  = p_fuel_pa * V_ullage_fuel_m3 / (R_s_jkgk * m_gas_fuel_kg); % [K]
            else
                dm_in_fuel_kg = 0;                                                   % [kg]
                T_gas_fuel_K  = T_gas_fuel_K + Q_fuel_J / (m_gas_fuel_kg * cp_n2_jkgk); % [K]
            end

            % ==========================================================
            %  Step 4: PRESSURANT TANK BLOWDOWN UPDATE
            % ==========================================================

            % Total mass leaving PT this timestep
            dm_out_kg    = dm_in_ox_kg + dm_in_fuel_kg;                             % [kg]
            mdot_out_kgs = dm_out_kg / dt_s;                                        % [kg/s]

            % Wall-to-gas heat transfer inside the pressurant tank
            % Positive when wall is warmer than gas (wall heats the expanding gas)
            Qdot_wall_W = h_wall_wm2k * A_wall_m2 * (T_wall_K - T_pt_K);          % [W]

            % Update pressurant gas temperature (forward Euler)
            %   m_pt * cv * dT_pt/dt = -mdot_out * Rs * T_pt + Qdot_wall
            %
            %   First term: cooling from mass outflow. Gas leaving carries
            %     enthalpy (cp*T), gas staying retains internal energy (cv*T).
            %     Net energy loss per unit mass leaving = (cp - cv)*T = Rs*T.
            %   Second term: heating from the aluminum tank wall.
            dT_pt_K = ((-mdot_out_kgs * R_s_jkgk * T_pt_K + Qdot_wall_W) ...
                       / (m_pt_kg * cv_n2_jkgk)) * dt_s;                           % [K]
            T_pt_K  = T_pt_K + dT_pt_K;                                            % [K]

            % Update wall temperature (lumped capacitance, forward Euler)
            %   Wall cools as it gives up heat to the expanding gas.
            %   Valid because Biot number << 1 for thin aluminum wall.
            dT_wall_K = (-Qdot_wall_W / (m_wall_kg * c_wall_jkgk)) * dt_s;        % [K]
            T_wall_K  = T_wall_K + dT_wall_K;                                      % [K]

            % Update pressurant tank mass
            m_pt_kg = m_pt_kg - dm_out_kg;                                          % [kg]

            % ==========================================================
            %  Step 5: STORE HISTORIES
            % ==========================================================
            if is_final_iter
                T_ox_history_K(i+1)           = T_gas_ox_K;
                T_fuel_history_K(i+1)         = T_gas_fuel_K;
                T_pt_history_K(i+1)           = T_pt_K;
                T_wall_history_K(i+1)         = T_wall_K;
                m_pt_history_kg(i+1)          = m_pt_kg;
                p_pt_history_pa(i+1)          = m_pt_kg * R_s_jkgk * T_pt_K / V_pt_m3; % [Pa]
                m_ox_ullage_history_kg(i+1)   = m_gas_ox_kg;
                m_fuel_ullage_history_kg(i+1) = m_gas_fuel_kg;
            end

        end % end time-stepping loop

        % --- Bisection update (custom COPV only) ---
        if ~P.use_cots_tank
            p_pt_final_pa = m_pt_kg * R_s_jkgk * T_pt_K / V_pt_m3;               % [Pa] Final PT pressure
            if p_pt_final_pa > p_final_target
                V_high = V_pt_m3;   % Surplus pressure -> tank is too big -> shrink upper bound
            else
                V_low = V_pt_m3;    % Pressure too low -> tank is too small -> raise lower bound
            end
        end

    end % end bisection loop

    % ================================================================
    %  V2 FINAL RESULTS
    % ================================================================
    m_gas_ox_thermal_kg    = m_gas_ox_kg;                                          % [kg] LOX ullage gas at end of burn
    m_gas_fuel_thermal_kg  = m_gas_fuel_kg;                                        % [kg] Fuel ullage gas at end of burn

    % Total mass loaded into the pressurant tank at fill time.
    % This is the gas in the PT at launch.  It includes:
    %   - Gas that will be delivered during the burn
    %   - Gas that remains in the PT at burnout (cannot be used)
    % It does NOT include the pre-pressed gas already in the ullages at launch.
    m_gas_total_thermal_kg = P.p_storage_pressurant_pa * V_pt_m3 / (R_s_jkgk * T_fill_K); % [kg]

    % ================================================================
    %  VERIFICATION CHECKS (printed to command window)
    % ================================================================

    % --- Mass Conservation Check ---
    % Total N2 in the system = PT mass + LOX ullage mass + Fuel ullage mass
    % This should remain constant throughout the burn (no gas leaves the system).
    m_total_n2_final_kg = m_pt_kg + m_gas_ox_kg + m_gas_fuel_kg;                   % [kg]
    mass_error_kg       = abs(m_total_n2_final_kg - m_total_n2_initial_kg);        % [kg]
    mass_error_pct      = mass_error_kg / m_total_n2_initial_kg * 100;             % [%]

    fprintf('\n');
    fprintf('  --- V2 Blowdown Verification Checks ---\n');
    fprintf('  [CHECK] Mass conservation error: %.2e kg (%.6f%%)\n', mass_error_kg, mass_error_pct);
    if mass_error_pct > 0.1
        warning('VAPOR:MassConservation', ...
            'N2 mass conservation error exceeds 0.1%%: %.4f kg (%.2f%%).', ...
            mass_error_kg, mass_error_pct);
    end

    % --- Monotonicity Check ---
    % PT pressure and temperature should decrease monotonically.
    if is_final_iter
        if any(diff(p_pt_history_pa) > 1e-6)   % allow tiny numerical noise
            warning('VAPOR:Monotonicity', ...
                'Pressurant tank pressure is not monotonically decreasing. Check for sign errors.');
            fprintf('  [CHECK] PT pressure monotonicity: FAILED\n');
        else
            fprintf('  [CHECK] PT pressure monotonicity: PASSED\n');
        end
    end

    % --- Minimum Pressure Check ---
    % For the COTS tank, verify that PT pressure stays above operating pressure
    % throughout the entire burn.  If it drops below, gas can no longer push
    % propellant out of the tanks.
    if P.use_cots_tank && is_final_iter
        p_pt_min_pa = min(p_pt_history_pa);                                        % [Pa] Minimum PT pressure during burn
        p_op_max_pa = max(p_ox_pa, p_fuel_pa);                                    % [Pa] Highest propellant tank operating pressure
        if p_pt_min_pa < p_op_max_pa
            idx_drop = find(p_pt_history_pa < p_op_max_pa, 1);                    % [index] First timestep where pressure is too low
            fprintf('  [CHECK] Minimum pressure: FAILED — PT drops to %.0f PSI at t=%.2f s (op pressure = %.0f PSI)\n', ...
                p_pt_min_pa * C.PA_TO_PSI, P.t_burn_s, p_op_max_pa * C.PA_TO_PSI);
            warning('VAPOR:PressureDrop', ...
                'PT pressure drops below operating pressure at t=%.2f s. COTS tank may be undersized.', ...
                t_history_s(idx_drop));
        else
            fprintf('  [CHECK] Minimum pressure: PASSED — Min PT = %.0f PSI > Op = %.0f PSI\n', ...
                p_pt_min_pa * C.PA_TO_PSI, p_op_max_pa * C.PA_TO_PSI);
        end
    end

    fprintf('  --- End V2 Checks ---\n\n');

else
    % --- Blowdown model disabled — use V1 results as primary output ---
    m_gas_ox_thermal_kg     = m_gas_ox_v1_final_kg;
    m_gas_fuel_thermal_kg   = m_gas_fuel_v1_final_kg;
    m_gas_total_thermal_kg  = m_gas_total_v1_kg;
    T_ox_history_K          = T_ox_v1_hist_K;
    T_fuel_history_K        = T_fuel_v1_hist_K;
    T_pt_history_K          = [];
    T_wall_history_K        = [];
    p_pt_history_pa         = [];
    m_pt_history_kg         = [];
    m_ox_ullage_history_kg  = [];
    m_fuel_ullage_history_kg= [];

end % end blowdown toggle

%% ====================================================================
%  STORAGE VOLUME CALCULATION
%
%  The required storage volume is always evaluated at FILL conditions
%  (294 K, 3000 PSI).  This answers the physical question: "Does the
%  required gas mass fit in the tank at the time of loading?"
%  See Section 6 of the blowdown analysis document.
%  ====================================================================

if P.use_blowdown_model
    % V2: the tank volume was either fixed (COTS) or found by bisection
    v_storage_thermal_m3 = V_pt_m3;                                                % [m^3]
else
    % V1 / Isothermal: compute from total mass at fill conditions
    n_total_thermal_mol  = m_gas_total_thermal_kg / P.pressurant_molar_mass_kgmol; % [mol]
    v_storage_thermal_m3 = (n_total_thermal_mol * C.r_universal_jmolk * T_fill_K) ...
                           / P.p_storage_pressurant_pa;                            % [m^3]
end

v_storage_thermal_L = v_storage_thermal_m3 * 1000;                                 % [L]

end