function [h_max_m, f_drag_max_n, maxq_pa, altitude_m, velocity_ms, acceleration_ms2, t_s, end_of_flight_event] = runFlightSimulation( ...
    m_total_kg, m_final_kg, m_dot_total_kgs, d_vehicle_outer_m, P, C)
%{
---------------------------------------------------------------------------
runFlightSimulation
  Runs a 1D time-stepping flight simulation using the ISA atmosphere model
  and RASAero drag data to estimate apogee and key aerodynamic loads.

Inputs:
  m_total_kg        - [kg]    Liftoff (wet) mass
  m_final_kg        - [kg]    Post-burnout (dry + residual) mass
  m_dot_total_kgs   - [kg/s]  Total propellant mass flow rate
  d_vehicle_outer_m - [m]     Vehicle outer diameter (for cross-sectional area)
  P  - [struct] Design inputs from getInputs()  (uses P.t_burn_s, P.f_thrust_n, P.rasaero_data_0_AoA)
  C  - [struct] Physical constants from getConstants()  (uses C.g_earth_ms2)

Outputs:
  h_max_m             - [m]      Estimated apogee altitude
  f_drag_max_n        - [N]      Maximum drag force experienced during flight
  maxq_pa             - [Pa]     Maximum dynamic pressure (Max-Q)
  altitude_m          - [m]      Altitude array over the simulation time
  velocity_ms         - [m/s]    Velocity array over the simulation time
  end_of_flight_event - [index]  Time index when the rocket returns to ground
---------------------------------------------------------------------------
%}

t_s                = 0:0.01:120;             % [s]     Time interval and step to analyze
altitude_m       = zeros(size(t_s));          % [m]     Pre-allocate altitude array
velocity_ms      = zeros(size(t_s));          % [m/s]   Pre-allocate velocity array
vehicle_mach     = zeros(size(t_s));          % [m/s]   Pre-allocate mach number array
acceleration_ms2 = zeros(size(t_s));          % [m/s^2] Pre-allocate acceleration array
f_drag_n         = zeros(size(t_s));          % [N]     Pre-allocate drag force array
thrust_n         = zeros(size(t_s)); thrust_n(t_s<=P.t_burn_s) = P.f_thrust_n; % [N] Make thrust array (should be zero after burnout)
m_total_sim_kg   = m_total_kg - m_dot_total_kgs .* t_s .* (thrust_n>0); m_total_sim_kg(t_s>P.t_burn_s) = m_final_kg; % [kg] Array of vehicle total mass over time
dt_s             = t_s(2) - t_s(1);            % [s]     time step
area_cross_m2    = pi * (d_vehicle_outer_m/2)^2; % [m^2] Vehicle cross-sectional area

% Step through time to calculate altitue, velocity, and acceleration using RK4
for n = 2:max(size(t_s))
    % Current state
    y0 = altitude_m(n-1);
    v0 = velocity_ms(n-1);
    
    % Mass and Thrust evaluated at current, mid, and next steps
    T_0   = thrust_n(n-1);
    m_0   = m_total_sim_kg(n-1);
    T_mid = (thrust_n(n-1) + thrust_n(n)) / 2;
    m_mid = (m_total_sim_kg(n-1) + m_total_sim_kg(n)) / 2;
    T_1   = thrust_n(n);
    m_1   = m_total_sim_kg(n);

    % --- Runge-Kutta 4th Order Steps ---
    % k1
    [a_k1, ~, ~, ~] = calc_accel(y0, v0, T_0, m_0, area_cross_m2, P, C);
    k1_y = v0;
    k1_v = a_k1;
    
    % k2
    [a_k2, ~, ~, ~] = calc_accel(y0 + k1_y * dt_s/2, v0 + k1_v * dt_s/2, T_mid, m_mid, area_cross_m2, P, C);
    k2_y = v0 + k1_v * dt_s/2;
    k2_v = a_k2;
    
    % k3
    [a_k3, ~, ~, ~] = calc_accel(y0 + k2_y * dt_s/2, v0 + k2_v * dt_s/2, T_mid, m_mid, area_cross_m2, P, C);
    k3_y = v0 + k2_v * dt_s/2;
    k3_v = a_k3;
    
    % k4
    [a_k4, ~, ~, ~] = calc_accel(y0 + k3_y * dt_s, v0 + k3_v * dt_s, T_1, m_1, area_cross_m2, P, C);
    k4_y = v0 + k3_v * dt_s;
    k4_v = a_k4;
    
    % --- Update State Arrays ---
    velocity_ms(n) = v0 + (dt_s / 6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);
    altitude_m(n)  = y0 + (dt_s / 6) * (k1_y + 2*k2_y + 2*k3_y + k4_y);
    
    % --- Record properties for the newly updated state ---
    [acceleration_ms2(n), f_drag_n(n), vehicle_mach(n), C_D] = calc_accel(altitude_m(n), velocity_ms(n), T_1, m_1, area_cross_m2, P, C);
    
    % end the simulation if the rocket reaches the ground
    if altitude_m(n) < 0
        break
    end
end

h_max_m      = max(altitude_m);                     % [m]  Estimated apogee
f_drag_max_n = max(f_drag_n);                        % [N]  Maximum drag force
maxq_pa      = f_drag_max_n/(C_D * area_cross_m2);  % [Pa] Maximum dynamic pressure (derived from F_drag = q * C_D * A)

end_of_flight_event = find(altitude_m(10:end)==0,1,'first'); % [index] Find time index when flight ends (rocket hits ground). Starts search after 10 steps to avoid t=0.

end

% =========================================================================
% Local Helper Function
% =========================================================================
function [a, f_drag, mach, cd] = calc_accel(y, v, thrust, mass, area_cross, P, C)
    % Bound y to >= 0 to prevent atmosisa from crashing if an RK4 intermediate step dips slightly negative
    [~, v_mach1, ~, rho, ~, ~] = atmosisa(max(0, y)); 
    mach = v / v_mach1;
    cd = interp1(P.rasaero_data_0_AoA.Mach, P.rasaero_data_0_AoA.CD, mach, "linear", "extrap");
    f_drag = -0.5 * rho * cd * area_cross * v * abs(v);
    a = (thrust + f_drag) / mass - C.g_earth_ms2;
end