%{
---------------------------------------------------------------------------
CU Sounding Rocket Lab - Liquid Propulsion Team
Tank Sizing and Vehicle Mass Model

Description:
This script performs a trade study for the design of the propellant and
pressurant tanks for the Taipan liquid rocket engine.

Authors: Hussain Almatruk, Jonathan Forte
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
o_f_ratio       = 2.5;         % [~] Oxidizer-to-Fuel mass ratio

% --- 1.2 Propellants & Pressurant ---
ox_density_kgm3   = 1141;                   % [kg/m^3] Density of Liquid Oxygen
fuel_density_kgm3 = 820;                    % [kg/m^3] Density of Jet-A
pressurant_molar_mass_kgmol = 0.0280134;    % [kg/mol] Molar mass of Nitrogen (N2)
p_op_ox_tank_pa = 850*6894.76;              % [Pa] Operating pressure of LOX tank TODO: Define this value
p_op_fuel_tank_pa = 850*6894.76;            % [Pa] Operating pressure of Fuel tank TODO: Define this value
p_storage_pressurant_pa = 3000*6894.76;     % [Pa] The pressure the Nitrogen is stored in dedicated tank (MEOP). TODO: Define this value
ox_temp_k         = 90;                     % [K] Temperature of LOX in tank TODO: Define this value
fuel_temp_k       = 294;                    % [K] Temperature of Jet-A in tank (Ambient probably)
pressurant_temp_k = 294;                    % [K] Temperature of pressurant gas in its tank
residual_fraction = 0.1;                    % [~] Percentage of fuel left over in tank after burnout

% --- 1.3 Vehicle Geometry & Materials ---
% Assumption: Both LOX and Fuel tanks are cylinders of the same diameter
D = 5*2.54/100;                     % [m] Diameter of the propellant tanks
d_ox_tank_m     = D;                % [m] Diameter of the LOX tank
d_fuel_tank_m   = D;                % [m] Diameter of the Fuel tank

% Material Properties for LOX Tank 
material_density_ox_kgm3         = 2840; % [kg/m^3] TODO: Define this value
material_allowable_stress_ox_pa  = 2.90*10^8; % [Pa] TODO: Define this value

% Material Properties for Fuel Tank 
material_density_fuel_kgm3       = 2840; % [kg/m^3] TODO: Define this value
material_allowable_stress_fuel_pa= 2.90*10^8; % [Pa] TODO: Define this value

% Material Properties for Pressurant Tank 
material_density_liner_kgm3             = 2840;     % [kg/m^3] TODO: Define this value
t_liner                                 = 0.003;    % [m] Thickness of COPV liner 
material_density_pressurant_kgm3        = 1800;     % [kg/m^3] TODO: Define this value
material_allowable_stress_pressurant_pa = 3.5*10^9; % [Pa] TODO: Define this value

% --- 1.4 Design Margins & Factors ---
safety_factor   = 1.5;        % [~] Safety factor for pressure vessels
% Joint efficiency can differ based on tank material and welding process
joint_efficiency_ox_tank        = 0.8; % [~] TODO: Define this value 
joint_efficiency_fuel_tank      = 0.8; % [~] TODO: Define this value 
joint_efficiency_pressurant_tank= 1.0; % [~] TODO: Define this value 

corrosion_allowance_m = 0.001;  % [m] Extra thickness for material degradation
ullage_fraction       = 0.1;    % [~] Percent of empty volume in tanks (e.g., 0.1 for 10%)

% --- 1.5 Estimated Masses (Non-Calculated) ---
m_misc_kg     = 30; % [kg] TODO: Estimate mass of payload, structure, fins, avionics, recovery
m_plumbing_kg = 10; % [kg] TODO: Estimate mass of valves and plumbing

% --- 1.6 Design Constraints ---

l_airframe_max_m = 3; %[m] Maximum allowable vehicle length TODO: Define this value 

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
v_total_ox_tank_m3 = v_ox_m3/(1-ullage_fraction);       % [m^3] Total internal volume of LOX tank
v_total_fuel_tank_m3 = v_fuel_m3/(1-ullage_fraction);   % [m^3] Total internal volume of Fuel tank

% Calculate tank design pressures (Eq. 10)
p_design_ox_pa = p_op_ox_tank_pa * safety_factor;     % [Pa] Design pressure for LOX tank
p_design_fuel_pa = p_op_fuel_tank_pa * safety_factor; % [Pa] Design pressure for Fuel tank

% Calculate cylinder wall thicknesses (Eq. 11, ASME Hoop Stress)
t_cyl_ox_m = (p_design_ox_pa * d_ox_tank_m) / (2 * material_allowable_stress_ox_pa * joint_efficiency_ox_tank) + corrosion_allowance_m; % [m] LOX tank cylinder thickness
t_cyl_fuel_m = (p_design_fuel_pa * d_fuel_tank_m) / (2 * material_allowable_stress_fuel_pa * joint_efficiency_fuel_tank) + corrosion_allowance_m; % [m] Fuel tank cylinder thickness

% Calculate tank internal radii (Basic Geometry)
r_ox_tank_m = d_ox_tank_m/2 - t_cyl_ox_m;       % [m] Internal radius of cylindrical ox tank
r_fuel_tank_m = d_fuel_tank_m/2 - t_cyl_fuel_m; % [m] Internal radius of cylindrical fuel tank

% Calculate end cap volumes (Eq. 7)
v_caps_ox_m3 = (4/3)*pi*r_ox_tank_m^3;     % [m^3] Combined volume of the two hemispherical end caps for the LOX tank
v_caps_fuel_m3 = (4/3)*pi*r_fuel_tank_m^3; % [m^3] Combined volume of the two hemispherical end caps for the Fuel tank

% Calculate cylindrical section volumes (Eq. 8)
v_cyl_ox_m3 = v_total_ox_tank_m3 - v_caps_ox_m3;     % [m^3] Volume of the LOX tank's cylindrical section
v_cyl_fuel_m3 = v_total_fuel_tank_m3 - v_caps_fuel_m3; % [m^3] Volume of the Fuel tank's cylindrical section

% Calculate cylindrical section lengths (Eq. 9)
l_cyl_ox_m = v_cyl_ox_m3 / (pi * r_ox_tank_m^2);     % [m] Required length of the LOX tank's cylindrical section
l_cyl_fuel_m = v_cyl_fuel_m3 / (pi * r_fuel_tank_m^2); % [m] Required length of the Fuel tank's cylindrical section

% Calculate end cap wall thicknesses (Eq. 12, Spherical Shell Stress)
% t_caps_ox_m = (p_design_ox_pa * d_ox_tank_m) / (4 * material_allowable_stress_ox_pa * joint_efficiency_ox_tank) + corrosion_allowance_m; % [m] LOX tank end cap thickness
% t_caps_fuel_m = (p_design_fuel_pa * d_fuel_tank_m) / (4 * material_allowable_stress_fuel_pa * joint_efficiency_fuel_tank) + corrosion_allowance_m; % [m] Fuel tank end cap thickness
t_caps_ox_m = t_cyl_ox_m;     % [m] LOX tank end cap thickness
t_caps_fuel_m = t_cyl_fuel_m; % [m] Fuel tank end cap thickness

% Calculate empty tank masses (Eq. 13)
m_empty_ox_tank_kg = material_density_ox_kgm3 * pi * ((D/2)^2-(D/2-t_cyl_ox_m)^2)*l_cyl_ox_m + 4/3*pi*((D/2)^3-(D/2-t_caps_ox_m)^3); % [kg] Mass of the empty LOX tank
m_empty_fuel_tank_kg = material_density_fuel_kgm3 * pi * ((D/2)^2-(D/2-t_cyl_fuel_m)^2)*l_cyl_fuel_m + 4/3 * pi * ((D/2)^3-(D/2-t_caps_fuel_m)^3); % [kg] Mass of the empty Fuel tank

% --- 2.3 - Pressurant System Sizing --- 

% Calculate moles of pressurant required to displace propellants (Eq. 14, Ideal Gas Law)
% Assumes pressurant gas cools to its storage temperature in the propellant tanks
n_ox_mol = (p_op_ox_tank_pa * v_ox_m3) / (r_universal_jmolk * ox_temp_k);       % [mol] Moles of gas to displace oxidizer
n_fuel_mol = (p_op_fuel_tank_pa * v_fuel_m3) / (r_universal_jmolk * fuel_temp_k); % [mol] Moles of gas to displace fuel

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
t_pressurant_tank_m = (p_design_pressurant_pa * (r_pressurant_tank_internal_m + t_liner)) / (2 * material_allowable_stress_pressurant_pa * joint_efficiency_pressurant_tank) + corrosion_allowance_m; % [m] Wall thickness of pressurant tank

% Calculate mass of empty pressurant tank (Eq. 20)
v_shell_pressurant_m3 = (4/3) * pi * ((r_pressurant_tank_internal_m + t_liner + t_pressurant_tank_m)^3 - (r_pressurant_tank_internal_m + t_liner)^3); % [m^3] Volume of the tank material
v_shell_liner_m3 = (4/3) * pi * ((r_pressurant_tank_internal_m + t_liner)^3 - (r_pressurant_tank_internal_m )^3); % [m^3] Volume of the tank liner material
m_empty_pressurant_tank_kg = material_density_pressurant_kgm3 * v_shell_pressurant_m3 + material_density_liner_kgm3 * v_shell_liner_m3; % [kg] Mass of the empty pressurant tank

% Calculate pressurant tank outer diameter
D_pressurant_tank_outer_m = 2 * (r_pressurant_tank_internal_m + t_liner + t_pressurant_tank_m);

if D_pressurant_tank_outer_m > D
    warning('Outer diameter of pressurant COPV exceeds maximum allowed diameter D. COPV converted to cylinder.');
end

if l_cyl_ox_m < 0 || l_cyl_fuel_m < 0
    warning('Propellant cylinder length is too short.');
end

% Change pressurant tank to cylinder if outer diameter too large
if D_pressurant_tank_outer_m > D
    t_pressurant_tank_m = (p_design_pressurant_pa * (r_pressurant_tank_internal_m + t_liner)) / (material_allowable_stress_pressurant_pa * joint_efficiency_pressurant_tank) + corrosion_allowance_m; % [m] Wall thickness of pressurant tank
    l_cyl_pressurant_m = (v_pressurant_tank_internal_m3 - 4*pi/3*(D/2 - t_pressurant_tank_m - t_liner)^3)/(pi * (D/2 - t_pressurant_tank_m - t_liner)^2); % [m] Length of cylindrical portion of pressurant tank
    m_empty_pressurant_tank_kg = 4*pi/3*(material_density_pressurant_kgm3 * ((D/2)^3 - (D/2 - t_pressurant_tank_m)^3) + material_density_liner_kgm3 * ((D/2 - t_pressurant_tank_m)^3 - (D/2 - t_pressurant_tank_m - t_liner)^3)) + pi*(material_density_pressurant_kgm3*((D/2)^2 - (D/2 - t_pressurant_tank_m)^2) + material_density_liner_kgm3*((D/2 - t_pressurant_tank_m)^2 - (D/2 - t_pressurant_tank_m - t_liner)^2)) * l_cyl_pressurant_m; % [kg] Mass of the empty pressurant tank 
end
% --- 2.4 - Vehicle Mass Buildup ---
% Placeholder for summing up all masses

m_total = m_empty_ox_tank_kg + m_empty_fuel_tank_kg + m_ox_kg + m_fuel_kg + m_empty_pressurant_tank_kg + m_pressurant_gas_kg + m_plumbing_kg + m_misc_kg;

% --- 2.5 - Performance Analysis ---
% Placeholder for TWR, Delta-V, and Apogee calculations

TWR = (f_thrust_n / m_total) / g_earth_ms2;
dV = i_sp_s * g_earth_ms2 * log(m_total/((m_total-(1-residual_fraction)*(m_ox_kg + m_fuel_kg))));

t = 0:0.001:120;                % [s] Time interval and step to analyze
altitude = zeros(size(t));      % [m] Pre-allocate altitude array
velocity = zeros(size(t));      % [m/s] Pre-allocate velocity array
acceleration = zeros(size(t));  % [m/s^2] Pre-allocate acceleration array
thrust = zeros(size(t)); thrust(t<=t_burn_s) = f_thrust_n;  % [N] Make thrust array (should be zero after burnout)
m_total_sim = m_total - m_dot_total_kgs .* t .* (thrust>0); m_total_sim(t>t_burn_s) = m_total - m_dot_total_kgs * t_burn_s; % [kg] Array of vehicle total mass over time
dt = t(2) - t(1);   % [s] time step
C_D = .5;           % [~] Vehicle drag coefficient
A = pi * (D/2)^2;   % [m^2] Vehicle cross-sectional area

% Step through time to calculate altitue, velocity, and acceleration through time
for i = 2:max(size(t))
velocity(i) = velocity(i-1) + acceleration(i-1) * dt;
altitude(i) = altitude(i-1) + velocity(i-1) * dt;
[~,~,~,rho,~,~] = atmosisa(altitude(i));
acceleration(i) = (thrust(i) - 0.5 * rho * C_D * A * velocity(i) * abs(velocity(i))) / m_total_sim(i) - g_earth_ms2;
% end the simulation if the rocket reaches the ground
if altitude(i) < 0
    break
end
end

endofflight = find(altitude(10:end)==0,1,'first'); % find time step when flight ends

%% 3.0 - OUTPUTS
% This section displays the final calculated values in a clean format.

% --- 3.1 - Vehicle Parameters ---

% --- 3.2 - Vehicle Performance ---
tiledlayout(2,2)

nexttile
plot(t(1:endofflight),altitude(1:endofflight)/0.3048);
grid on
ylabel('Altitude [ft]')
xlabel('Time [s]')
title('Altitude vs. Time')

nexttile
plot(t(1:endofflight),velocity(1:endofflight));
grid on
ylabel('Velocity [m/s]')
xlabel('Time [s]')
title('Velocity vs. Time')

nexttile
plot(t(1:endofflight),acceleration(1:endofflight)./g_earth_ms2);
grid on
ylabel('Acceleration [g]')
xlabel('Time [s]')
title('Acceleration vs. Time')

parameter = {'Vehicle Diameter [in]'; 'Burn Time [s]'; 'Maximum Altitude [ft]'; 'Maximum Velocity [m/s]'; 'Maximum Acceleration [g]'};
value = [D*100/2.54; t_burn_s; max(altitude)/0.3048; max(velocity); max(acceleration./g_earth_ms2)];
T = table(value, 'RowNames',parameter);

uitable("Data",T{:,:},'RowName',T.Properties.RowNames,'ColumnName', {'Value'},'Units', 'Normalized', 'Position',[0.6, 0.2, 0.3, 0.212]);

%% 4.0 - Test Section

% you can use this section temporarly to test that github works for you and the changes you make are actually working