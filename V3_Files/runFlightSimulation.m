function [h_max_m, f_drag_max_n, maxq_pa, altitude_m, velocity_ms, end_of_flight_event] = runFlightSimulation( ...
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

t                = 0:0.001:120;             % [s]     Time interval and step to analyze
altitude_m       = zeros(size(t));          % [m]     Pre-allocate altitude array
velocity_ms      = zeros(size(t));          % [m/s]   Pre-allocate velocity array
vehicle_mach     = zeros(size(t));          % [m/s]   Pre-allocate mach number array
acceleration_ms2 = zeros(size(t));          % [m/s^2] Pre-allocate acceleration array
f_drag_n         = zeros(size(t));          % [N]     Pre-allocate drag force array
thrust_n         = zeros(size(t)); thrust_n(t<=P.t_burn_s) = P.f_thrust_n; % [N] Make thrust array (should be zero after burnout)
m_total_sim_kg   = m_total_kg - m_dot_total_kgs .* t .* (thrust_n>0); m_total_sim_kg(t>P.t_burn_s) = m_final_kg; % [kg] Array of vehicle total mass over time
dt_s             = t(2) - t(1);            % [s]     time step
area_cross_m2    = pi * (d_vehicle_outer_m/2)^2; % [m^2] Vehicle cross-sectional area

% Step through time to calculate altitue, velocity, and acceleration through time
for n = 2:max(size(t))
    velocity_ms(n)      = velocity_ms(n-1) + acceleration_ms2(n-1) * dt_s;                                                                            % [m/s]   Velocity array
    altitude_m(n)       = altitude_m(n-1)  + velocity_ms(n-1) * dt_s;                                                                                 % [m]     Altitude array
    [~,v_mach1_ms,~,rho_kgm3,~,~] = atmosisa(altitude_m(n));                                                                                          % [kg/m^3] Get atmospheric density and speed of sound at current altitude using ISA model
    vehicle_mach(n)     = velocity_ms(n)./v_mach1_ms;                                                                                                  % [unitless] Calculate mach number
    C_D                 = interp1(P.rasaero_data_0_AoA.Mach, P.rasaero_data_0_AoA.CD, vehicle_mach(n),"linear","extrap");                             % [unitless] Calculate current C_D as a function of mach number and angle of attack
    f_drag_n(n)         = - 0.5 * rho_kgm3 * C_D * area_cross_m2 * velocity_ms(n) * abs(velocity_ms(n));                                             % [N]     Vehicle drag (using v*abs(v) to ensure drag always opposes velocity)
    acceleration_ms2(n) = (thrust_n(n) + f_drag_n(n)) / m_total_sim_kg(n) - C.g_earth_ms2;                                                           % [m/s^2] Vehicle acceleration
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
